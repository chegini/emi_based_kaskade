#ifndef EMI_BDDC_HH
#define EMI_BDDC_HH

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <boost/fusion/include/at_c.hpp>

#include "fem/assemble.hh"
#include "fem/gridBasics.hh"
#include "fem/variables.hh"
#include "mg/bddc.hpp"
#include "timestepping/semieuler.hh"
#include "timestepping/sdc.hh"
#include "utilities/threading.hh"

namespace EmiBddc
{
  struct CompressionReport
  {
    std::uint64_t originalSentBytes = 0;
    std::uint64_t originalReceivedBytes = 0;
    std::uint64_t sentBytes = 0;
    std::uint64_t receivedBytes = 0;
    double encodeTimeMs = 0.0;
    double decodeTimeMs = 0.0;
    std::uint64_t encodeCount = 0;
    std::uint64_t decodeCount = 0;
    std::uint64_t codebookBits = 0;

    template <class Traffic>
    void add(Traffic const& traffic)
    {
      originalSentBytes += traffic.originalSentBytes;
      originalReceivedBytes += traffic.originalReceivedBytes;
      sentBytes += traffic.sentBytes;
      receivedBytes += traffic.receivedBytes;
      encodeTimeMs += traffic.encodeTimeMs;
      decodeTimeMs += traffic.decodeTimeMs;
      encodeCount += traffic.encodeCount;
      decodeCount += traffic.decodeCount;
      codebookBits = std::max(codebookBits, traffic.codebookBits);
    }

    void addUncompressed(std::array<int,6> const& traffic)
    {
      std::uint64_t const bytes = static_cast<std::uint64_t>(traffic[3] + traffic[4]);
      originalSentBytes += bytes;
      originalReceivedBytes += bytes;
      sentBytes += bytes;
      receivedBytes += bytes;
    }

    void print(std::string const& label) const
    {
      std::uint64_t const raw = originalSentBytes;
      std::uint64_t const wire = sentBytes;
      std::int64_t const saved = static_cast<std::int64_t>(raw)
                               - static_cast<std::int64_t>(wire);
      double const ratio = wire > 0 ? static_cast<double>(raw)/wire : 1.0;
      double const percent = raw > 0 ? 100.0*static_cast<double>(saved)/raw : 0.0;
      std::cout << "Compression report (" << label << "):\n"
                << "  raw sent bytes: " << raw << "\n"
                << "  transmitted sent bytes: " << wire << "\n"
                << "  bits saved: " << saved*8 << "\n"
                << "  savings: " << percent << "%\n"
                << "  ratio: " << ratio << ":1\n"
                << "  shared codebook estimate: " << codebookBits << " bits\n"
                << "  raw received bytes: " << originalReceivedBytes << "\n"
                << "  received bytes: " << receivedBytes << "\n";
    }
  };

  struct ProfileTimes
  {
    using Clock = std::chrono::steady_clock;

    double residualAssembly = 0.0;
    double localCollocationSystem = 0.0;
    double interfaceSetup = 0.0;
    double subdomainConstruction = 0.0;
    double transferSetup = 0.0;
    double solverConstruction = 0.0;
    double rhsSetup = 0.0;
    double bddcSolve = 0.0;
    double sdcUpdate = 0.0;
    double adaptivity = 0.0;
    double localOperatorExtraction = 0.0;
    size_t bddcSystems = 0;
    size_t bddcSubdomainConstructions = 0;
    size_t bddcSolveCalls = 0;

    static double seconds(Clock::time_point start)
    {
      return std::chrono::duration<double>(Clock::now()-start).count();
    }

    void print(char const* label) const
    {
      std::cout << "Profile (wall seconds, " << label << "):\n"
                << "  local operator extraction: " << localOperatorExtraction << "\n"
                << "  residual assembly: " << residualAssembly << "\n"
                << "  local collocation matrix/RHS: " << localCollocationSystem << "\n"
                << "  interface constraint setup: " << interfaceSetup << "\n"
                << "  subdomain construction/factorization: " << subdomainConstruction << "\n"
                << "  transfer configuration: " << transferSetup << "\n"
                << "  BDDC solver/coarse setup: " << solverConstruction << "\n"
                << "  BDDC RHS setup: " << rhsSetup << "\n"
                << "  BDDC subdomain constructions: " << bddcSubdomainConstructions << "\n"
                << "  BDDC solve calls: " << bddcSolveCalls
                << " across " << bddcSystems << " interval systems (wall seconds: "
                << bddcSolve << ", average calls/system: "
                << (bddcSystems > 0 ? static_cast<double>(bddcSolveCalls)/bddcSystems : 0.0)
                << ")\n"
                << "  SDC update/bookkeeping: " << sdcUpdate << "\n"
                << "  algebraic adaptivity: " << adaptivity << "\n";
    }
  };

  template <class Matrix, class Vector>
  struct BddcData
  {
    std::vector<int> tags;  // Material tag identifying each subdomain.
    std::vector<std::vector<size_t>> localDofs;  // Global DOF indices, in local order, per subdomain.
    std::vector<std::vector<Kaskade::BDDC::LocalDof>> sharedDofs;  // Local copies of each shared global DOF.
    std::vector<std::vector<int>> dofSubdomains;  // Inverse map: global DOF to owning subdomains.
    std::vector<std::vector<size_t>> dofCells;  // Inverse map: global DOF to incident leaf-cell indices.
    std::vector<int> subdomainSizes;  // Number of local DOFs per subdomain.
    // Combined local LHS for non-SDC BDDC; SDC+BDDC instead restricts mass and stiffness separately.
    std::vector<Matrix> localMatrices;
    // Weights for global-to-local residuals and local-to-global corrections; copies of a DOF sum to one.
    std::vector<Vector> weights;
  };

  // Build partition metadata and overlap weights independently of the operator matrix.
  template <class Grid, class Space, class Material, class Matrix, class Vector>
  BddcData<Matrix,Vector> buildBddcData(Grid const& grid,
                                        Space const& space,
                                        Material const& material,
                                        size_t nDofs,
                                        int nTasks)
  {
    using Kaskade::BDDC::LocalDof;

    // Cells with one material tag form a subdomain; interface DOFs consequently have multiple local copies.
    std::map<int,std::set<size_t>> dofsByTag;
    Dune::FieldVector<double,Grid::dimension> zero(0.0);

    for (auto const& cell : Dune::elements(grid.leafGridView()))
    {
      int const tag = material.value(cell,zero);
      auto const dofs = space.mapper().globalIndices(cell);
      auto& tagDofs = dofsByTag[tag];
      tagDofs.insert(dofs.begin(),dofs.end());
    }

    if (dofsByTag.empty())
      throw std::runtime_error("BDDC setup found no material subdomains.");

    BddcData<Matrix,Vector> data;
    data.dofCells.resize(nDofs);
    for (auto const& cell : Dune::elements(grid.leafGridView()))
    {
      size_t const cellIndex = space.indexSet().index(cell);
      for (size_t dof : space.mapper().globalIndices(cell))
        if (dof < nDofs)
          data.dofCells[dof].push_back(cellIndex);
    }
    for (auto& cells : data.dofCells)
    {
      std::sort(cells.begin(),cells.end());
      cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    }

    data.tags.reserve(dofsByTag.size());
    data.localDofs.reserve(dofsByTag.size());
    data.subdomainSizes.reserve(dofsByTag.size());

    std::vector<std::vector<LocalDof>> dofOccurrences(nDofs);

    for (auto const& [tag,dofs] : dofsByTag)
    {
      int const subdomain = static_cast<int>(data.tags.size());
      data.tags.push_back(tag);
      data.localDofs.emplace_back(dofs.begin(),dofs.end());
      if (data.localDofs.back().size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("BDDC local subdomain is too large for int local indices.");
      data.subdomainSizes.push_back(static_cast<int>(data.localDofs.back().size()));

      for (size_t local = 0; local < data.localDofs.back().size(); ++local)
        dofOccurrences[data.localDofs.back()[local]].push_back({subdomain,static_cast<int>(local)});
    }

    for (auto const& occurrence : dofOccurrences)
      if (occurrence.size() > 1)
        data.sharedDofs.push_back(occurrence);

    data.dofSubdomains.resize(nDofs);
    for (size_t global = 0; global < dofOccurrences.size(); ++global)
      for (auto const& occurrence : dofOccurrences[global])
        data.dofSubdomains[global].push_back(occurrence.s);

    if (data.sharedDofs.empty() && data.tags.size() > 1)
      throw std::runtime_error("BDDC setup found multiple subdomains but no shared dofs.");

    data.weights.resize(data.localDofs.size());
    size_t const taskLimit = nTasks > 0 ? static_cast<size_t>(nTasks)
                                        : std::numeric_limits<size_t>::max();
    auto buildSubdomainData = [&](size_t subdomain)
    {
      Vector weight(data.localDofs[subdomain].size());
      weight = 1.0;
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        size_t const ownerCount = std::max<size_t>(1,dofOccurrences[global].size());
        weight[local][0] = 1.0/static_cast<double>(ownerCount);
      }
      data.weights[subdomain] = std::move(weight);
    };
    if (taskLimit > 1)
      Kaskade::parallelFor(0,data.localDofs.size(),buildSubdomainData,taskLimit);
    else
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        buildSubdomainData(subdomain);

    std::cout << "BDDC subdomains: " << data.tags.size()
              << ", shared dofs: " << data.sharedDofs.size() << "\n";
    return data;
  }

  template <class Matrix, class Vector>
  std::vector<Matrix> extractLocalMatrices(BddcData<Matrix,Vector> const& data,
                                           Matrix const& globalMatrix,
                                           int nTasks)
  {
    std::vector<Matrix> localMatrices(data.localDofs.size());
    size_t const taskLimit = nTasks > 0 ? static_cast<size_t>(nTasks)
                                        : std::numeric_limits<size_t>::max();
    auto extractSubdomainMatrix = [&](size_t subdomain)
    {
      localMatrices[subdomain] = Matrix(data.localDofs[subdomain],globalMatrix);
    };
    // Each subdomain writes to a separate output slot; the source operator is read-only.
    if (taskLimit > 1)
      Kaskade::parallelFor(0,data.localDofs.size(),extractSubdomainMatrix,taskLimit);
    else
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        extractSubdomainMatrix(subdomain);
    return localMatrices;
  }

  template <class Matrix, class Vector>
  std::vector<Vector> scatterToSubdomains(BddcData<Matrix,Vector> const& data,
                                          Vector const& global)
  {
    // Copy global entries to local vectors, scaling shared DOFs so their local copies sum to the global value.
    std::vector<Vector> local(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
    {
      local[subdomain] = Vector(data.localDofs[subdomain].size());
      local[subdomain] = 0.0;
      for (size_t i = 0; i < data.localDofs[subdomain].size(); ++i)
      {
        size_t const globalDof = data.localDofs[subdomain][i];
        local[subdomain][i] = data.weights[subdomain][i][0] * global[globalDof];
      }
    }
    return local;
  }

  template <class Matrix, class Vector>
  Vector combineSubdomainVector(BddcData<Matrix,Vector> const& data,
                                std::vector<Vector> const& local,
                                size_t nDofs)
  {
    // Weighted assembly is the reverse operation to scatterToSubdomains.
    Vector global(nDofs);
    global = 0.0;
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      for (size_t i = 0; i < data.localDofs[subdomain].size(); ++i)
      {
        size_t const globalDof = data.localDofs[subdomain][i];
        global[globalDof] += data.weights[subdomain][i][0] * local[subdomain][i];
      }
    return global;
  }

  template <class Matrix>
  void addMatrixNeighborhood(Matrix const& A, std::vector<std::vector<size_t>>& dofNeighborhood)
  {
    for (size_t row = 0; row < A.N(); ++row)
    {
      dofNeighborhood[row].push_back(row);
      for (auto col = A[row].begin(); col != A[row].end(); ++col)
      {
        dofNeighborhood[row].push_back(col.index());
        if (col.index() < dofNeighborhood.size())
          dofNeighborhood[col.index()].push_back(row);
      }
    }
  }

  template <class Matrix>
  std::vector<std::vector<size_t>> buildMatrixNeighborhood(Matrix const& mass, Matrix const& stiffness)
  {
    // Algebraic adaptivity must retain couplings from either operator, not just the mass graph.
    std::vector<std::vector<size_t>> dofNeighborhood(mass.N());
    addMatrixNeighborhood(mass,dofNeighborhood);
    addMatrixNeighborhood(stiffness,dofNeighborhood);

    for (auto& neighbors : dofNeighborhood)
    {
      std::sort(neighbors.begin(),neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());
    }

    return dofNeighborhood;
  }

  inline std::vector<int> allSubdomains(size_t nSubdomains)
  {
    std::vector<int> activeIds(nSubdomains);
    std::iota(activeIds.begin(),activeIds.end(),0);
    return activeIds;
  }

  // Configure compression only for transfer types that provide the optional
  // compression API. The ordinary SpaceTransfer remains unchanged.
  template <class Transfer, class Options>
  void configureTransfer(Transfer& transfer, Options const& options)
  {
    if constexpr (requires(Transfer& candidate) {
                    candidate.setQuantizationBits(0);
                    candidate.setRestrictEncoding(true);
                    candidate.setProlongateEncoding(true);
                    candidate.setRestrictTransform(true);
                    candidate.setProlongateTransform(true);
                    candidate.setRestrictBitlengthEncoding(true);
                    candidate.setProlongateBitlengthEncoding(true);
                  })
    {
      transfer.setQuantizationBits(options.bddcCompressionBits);
      transfer.enableTransform(options.bddcGraphLifting
                                 ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                                 : Kaskade::BDDC::TransformType::NONE);
      transfer.setRestrictEncoding(options.bddcHuffman != 0);
      transfer.setProlongateEncoding(options.bddcHuffman != 0);
      transfer.setRestrictTransform(options.bddcGraphLifting != 0);
      transfer.setProlongateTransform(options.bddcGraphLifting != 0);
      transfer.setRestrictBitlengthEncoding(options.bddcBitlength != 0);
      transfer.setProlongateBitlengthEncoding(options.bddcBitlength != 0);
    }
  }

  // Construct independent BDDC subdomains concurrently, then move them into the
  // ordered vector required by BDDCSolver. Each worker writes to a distinct slot.
  template <class Subdomain, class Matrix, class Interfaces>
  std::vector<Subdomain> constructSubdomains(std::vector<Matrix> const& localOperators,
                                             Interfaces const& interfaces,
                                             int nThreads)
  {
    std::vector<std::optional<Subdomain>> slots(localOperators.size());
    std::exception_ptr constructionError;
    std::mutex constructionErrorMutex;
    auto construct = [&](size_t subdomain)
    {
      try
      {
        slots[subdomain].emplace(static_cast<int>(subdomain),localOperators[subdomain],interfaces);
        if constexpr (requires(Subdomain& subdomainObject, Matrix const& matrix) {
                        subdomainObject.transfer().setGraphLiftingFromMatrix(matrix);
                      })
          slots[subdomain]->transfer().setGraphLiftingFromMatrix(localOperators[subdomain]);
      }
      catch (...)
      {
        std::lock_guard<std::mutex> lock(constructionErrorMutex);
        if (!constructionError)
          constructionError = std::current_exception();
      }
    };

    size_t const taskLimit = nThreads > 0 ? static_cast<size_t>(nThreads)
                                         : std::numeric_limits<size_t>::max();
    if (taskLimit > 1)
      Kaskade::parallelFor(0,localOperators.size(),construct,taskLimit);
    else
      for (size_t subdomain = 0; subdomain < localOperators.size(); ++subdomain)
        construct(subdomain);

    if (constructionError)
      std::rethrow_exception(constructionError);

    std::vector<Subdomain> subdomains;
    subdomains.reserve(localOperators.size());
    for (auto& slot : slots)
    {
      if (!slot)
        throw std::runtime_error("BDDC subdomain construction did not populate every slot.");
      subdomains.emplace_back(std::move(*slot));
    }
    return subdomains;
  }

  template <class Matrix, class Vector, class Options>
  std::vector<int> selectActiveSubdomains(BddcData<Matrix,Vector> const& data,
                                          std::vector<std::vector<Vector>> const& corrections,
                                          std::vector<int> const& activeIds,
                                          std::vector<std::vector<size_t>> const& dofNeighborhood,
                                          double sdcContraction,
                                          Options const& options)
  {
    if (!options.algebraicAdaptivity || options.algebraicAdaptivityTolerance <= 0.0)
      return allSubdomains(data.localDofs.size());

    // Apply the correction-based AA estimate, then expand touched DOFs to complete owner subdomains.
    std::set<int> nextActiveIds;
    for (int subdomain : activeIds)
    {
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        double duMax = 0.0;
        for (auto const& pointCorrections : corrections)
          duMax = std::max(duMax,std::abs(pointCorrections[subdomain][local][0]));

        bool const selected = sdcContraction >= 1.0
                           || sdcContraction*duMax/(1.0-sdcContraction) > options.algebraicAdaptivityTolerance;
        if (!selected)
          continue;

        // BDDC subproblems are assembled and constrained subdomain-wise.
        // Therefore AA may detect activity at one dof, but the next solve must
        // activate every local dof of each touched owner subdomain.
        size_t const globalDof = data.localDofs[subdomain][local];
        for (size_t neighborDof : dofNeighborhood[globalDof])
          for (int owner : data.dofSubdomains[neighborDof])
            nextActiveIds.insert(owner);
      }
    }

    return std::vector<int>(nextActiveIds.begin(),nextActiveIds.end());
  }

// Evaluate node residuals in local subdomain ordering. Under AA, assemble cells incident to active DOFs.
template <class Functional, class State, class StateU, class TimeGrid, class Assembler, class Matrix, class Vector, class Options, class IndexSet>
void computeBddcSdcResiduals(Functional& F,
                               Kaskade::SemiImplicitEulerStep<Functional>& equation,
                               State const& stateAtStart,
                               std::vector<StateU> const& collocationStates,
                               TimeGrid const& grid,
                               BddcData<Matrix,Vector> const& data,
                               size_t nDofs,
                               int sweep,
                               double t,
                               double dt,
                               Assembler& assembler,
                               Options const& options,
                               std::vector<int> const& activeIds,
                               IndexSet const& indexSet,
                               std::vector<std::vector<Vector>>& residuals)
  {
    using SemiLinearization = Kaskade::SemiLinearizationAtInner<Kaskade::SemiImplicitEulerStep<Functional>>;

    State stateTmp(stateAtStart);
    State dstateTmp(stateAtStart);
    dstateTmp *= 0.0;
    equation.setTau(dt);
    F.Mass_stiff(0);
    Vector fullResidual(nDofs);

    bool const allSubdomainsActive = activeIds.size() == data.localDofs.size();
    std::vector<char> activeCells;
    if (!allSubdomainsActive)
    {
      activeCells.assign(indexSet.size(0),0);
      // Include every cell touching an active local DOF so its assembled residual row is complete.
      for (int subdomain : activeIds)
        for (size_t globalDof : data.localDofs[subdomain])
          for (size_t cell : data.dofCells[globalDof])
            activeCells[cell] = 1;
    }

    auto const& points = grid.points();
    for (int i = 0; i < points.N(); ++i)
    {
      // This autonomous EMI model has identical first-sweep guesses at every node, so the residual can be reused.
      if (sweep == 0 && i > 0)
      {
        residuals[i] = residuals[i-1];
        continue;
      }

      F.time(t + points[i] - points[0]);
      boost::fusion::at_c<0>(stateTmp.data) = collocationStates[i];
      if (allSubdomainsActive)
      {
        assembler.assemble(SemiLinearization(equation,stateTmp,stateTmp,dstateTmp),
                           Assembler::RHS,
                           options.assemblyThreads);
      }
      else
      {
        auto cellFilter = [&activeCells,&indexSet](auto const& cell)
        {
          return activeCells[indexSet.index(cell)] != 0;
        };
        assembler.template assemble<Kaskade::AssemblyDetail::TakeAllBlocks>(
                           SemiLinearization(equation,stateTmp,stateTmp,dstateTmp),
                           cellFilter,
                           Assembler::RHS,
                           options.assemblyThreads);
      }

      auto rhs = assembler.rhs();
      rhs.write(fullResidual.begin());
      // Remove the SemiImplicitEuler dt factor so the SDC residual has the unscaled spatial-operator convention.
      fullResidual *= (1.0/dt);
      residuals[i] = scatterToSubdomains(data,fullResidual);
    }
  }

  template <class Matrix, class StateU, class Vector>
  void computeBddcMassDifferences(BddcData<Matrix,Vector> const& data,
                                  std::vector<Matrix> const& localMassMatrices,
                                  std::vector<StateU> const& collocationStates,
                                  std::vector<std::vector<Vector>>& massDifferences)
  {
    // Compute the local restriction of M (u_i - u_{i+1}) for each collocation interval.
    int const intervals = static_cast<int>(collocationStates.size())-1;

    for (int i = 0; i < intervals; ++i)
      for (auto& localDifference : massDifferences[i])
        localDifference = 0.0;

    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      for (size_t localRow = 0; localRow < data.localDofs[subdomain].size(); ++localRow)
        for (auto col = localMassMatrices[subdomain][localRow].begin();
             col != localMassMatrices[subdomain][localRow].end(); ++col)
        {
          size_t const globalColumn = data.localDofs[subdomain][col.index()];
          for (int i = 0; i < intervals; ++i)
          {
            auto const diff = collocationStates[i].coefficients()[globalColumn]
                            - collocationStates[i+1].coefficients()[globalColumn];
            massDifferences[i][subdomain][localRow] += *col * diff;
          }
        }
  }

  // Form the SDC collocation LHS J = M - collocationWeight * A for one subdomain.
  template <class Matrix>
  void combineSdcMatrix(Matrix& J,
                        Matrix const& mass,
                        Matrix const& stiffness,
                        double stiffnessWeight)
  {
    J = mass;
    for (size_t row = 0; row < J.N(); ++row)
    {
      auto colJ = J[row].begin();
      auto const endJ = J[row].end();
      auto colA = stiffness[row].begin();
      auto const endA = stiffness[row].end();

      while (colJ != endJ && colA != endA)
      {
        if (colJ.index() == colA.index())
        {
          *colJ -= stiffnessWeight * *colA;
          ++colJ;
          ++colA;
        }
        else if (colJ.index() < colA.index())
          ++colJ;
        else
          ++colA;
      }
    }
  }

  template <class Transfer, class Matrix, class Vector, class Options>
  typename Matrix::field_type sdcIterationStepBddc(Kaskade::SDCTimeGrid const& grid,
                                                   Kaskade::SDCTimeGrid::RealMatrix const& Shat,
                                                   BddcData<Matrix,Vector> const& data,
                                                   std::vector<Matrix> const& localMassMatrices,
                                                   std::vector<Matrix> const& localStiffnessMatrices,
                                                   std::vector<std::vector<Vector>> const& residuals,
                                                   std::vector<std::vector<Vector>> const& massDifferences,
                                                   std::vector<std::vector<Vector>>& corrections,
                                                   std::vector<int> const& activeIds,
                                                   Options const& options,
                                                   ProfileTimes* profile = nullptr,
                                                   CompressionReport* compressionReport = nullptr)
  {
    using BddcSubdomain = Kaskade::BDDC::Subdomain<1,double,double,Transfer>;

    auto const& points = grid.points();
    int const intervals = points.size()-1;
    auto const& S = grid.integrationMatrix();
    typename Matrix::field_type norm = 0.0;

    for (auto& localCorrection : corrections[0])
      localCorrection = 0.0;

    auto const interfaceSetupStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
    Kaskade::BDDC::InterfaceAverages<1,int> interfaceAverages(data.sharedDofs,
                                                              data.subdomainSizes,
                                                              options.bddcInterfaceTypes);
    if (profile)
      profile->interfaceSetup += ProfileTimes::seconds(interfaceSetupStart);

    std::vector<int> solverActiveIds(activeIds);
    bool useCgSolver = options.bddcUseCg != 0;
    bool verbose = options.bddcVerbose != 0;

    for (int i = 1; i <= intervals; ++i)
    {
      if (profile)
        ++profile->bddcSystems;
      auto const systemStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      // The diagonal SDC coefficient varies by interval; restricted operators are reused across intervals.
      std::vector<Matrix> localJ;
      localJ.reserve(data.localDofs.size());
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      {
        localJ.push_back(localMassMatrices[subdomain]);
        combineSdcMatrix(localJ.back(),
                         localMassMatrices[subdomain],
                         localStiffnessMatrices[subdomain],
                         Shat[i-1][i]);
      }

      std::vector<Vector> rhs(data.localDofs.size());
      std::vector<Vector> tmp(data.localDofs.size());
      // Build the interval RHS from the mass jump, residual quadrature, and earlier sweep corrections.
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      {
        rhs[subdomain] = massDifferences[i-1][subdomain];
        localMassMatrices[subdomain].umv(corrections[i-1][subdomain],rhs[subdomain]);

        for (int j = 0; j <= intervals; ++j)
          rhs[subdomain].axpy(S[i-1][j],residuals[j][subdomain]);

        tmp[subdomain] = Vector(data.localDofs[subdomain].size());
        tmp[subdomain] = 0.0;
        for (int j = 0; j < i; ++j)
          tmp[subdomain].axpy(Shat[i-1][j],corrections[j][subdomain]);
        localStiffnessMatrices[subdomain].umv(tmp[subdomain],rhs[subdomain]);
      }
      if (profile)
      {
        profile->localCollocationSystem += ProfileTimes::seconds(systemStart);
      }

      auto const subdomainStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      if (profile)
        profile->bddcSubdomainConstructions += data.localDofs.size();

      // Each construction factors one interval matrix; these independent factorizations
      // are parallelized, but the resulting mutable subdomains are not reused across systems.
      auto subdomains = constructSubdomains<BddcSubdomain>(localJ,interfaceAverages,
                                                            options.assemblyThreads);
      if (profile)
        profile->subdomainConstruction += ProfileTimes::seconds(subdomainStart);

      auto const transferStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      for (auto& subdomain : subdomains)
        configureTransfer(subdomain.transfer(),options);
      if (profile)
        profile->transferSetup += ProfileTimes::seconds(transferStart);

      auto const solverStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      Kaskade::BDDC::BDDCSolver<BddcSubdomain> solver(subdomains,
                                                      interfaceAverages.coarseConstraints(),
                                                      solverActiveIds,
                                                      useCgSolver,
                                                      verbose);
      if (profile)
        profile->solverConstruction += ProfileTimes::seconds(solverStart);

      auto const rhsSetupStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      solver.setRhs(rhs);
      if (profile)
        profile->rhsSetup += ProfileTimes::seconds(rhsSetupStart);

      // Solve the collocation interval system to the requested BDDC residual tolerance.
      auto const solveStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      double residual = std::numeric_limits<double>::infinity();
      int iteration = 0;
      for (; iteration < options.bddcIterations; ++iteration)
      {
        residual = solver.solve();
        if (profile)
          ++profile->bddcSolveCalls;
        if (residual < options.bddcTolerance)
          break;
      }
      if (compressionReport)
      {
        if (options.bddcCompression)
          compressionReport->add(solver.compressionTraffic());
        else
          compressionReport->addUncompressed(solver.traffic());
      }
      if (profile)
        profile->bddcSolve += ProfileTimes::seconds(solveStart);

      for (int subdomain : activeIds)
        corrections[i][subdomain] = subdomains[subdomain].getSolution();

      if (verbose)
        std::cout << "    BDDC SDC interval " << i << "/" << intervals
                  << ": iterations=" << std::min(iteration+1,options.bddcIterations)
                  << ", residual=" << residual << "\n";

      for (int subdomain : activeIds)
        norm += (points[i]-points[i-1]) * (corrections[i][subdomain] * rhs[subdomain]);
    }

    return std::sqrt(std::max<typename Matrix::field_type>(0.0,norm/(points[intervals]-points[0])));
  }

template <class Transfer, class Functional, class Spaces, class State, class Element, class Matrix, class Vector, class Options, class IndexSet>
State runBddcSdc(Functional& F,
                   Spaces const& spaces,
                   State state,
                   Element&,
                   BddcData<Matrix,Vector> const& data,
                   Matrix const& mass,
                   Matrix const& stiffness,
                   size_t nDofs,
                   int steps,
                   Options const& options,
                   IndexSet const& indexSet)
  {
    using SemiLinearization = Kaskade::SemiLinearizationAtInner<Kaskade::SemiImplicitEulerStep<Functional>>;
    using Assembler = Kaskade::VariationalFunctionalAssembler<SemiLinearization>;
    using StateU = typename boost::fusion::result_of::value_at_c<typename State::Sequence,0>::type;

    ProfileTimes profile;
    auto const extractionStart = options.profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
    // The operators are fixed for this model, so extract their local blocks once and reuse them.
    auto const localMassMatrices = extractLocalMatrices(data,mass,options.assemblyThreads);
    auto const localStiffnessMatrices = extractLocalMatrices(data,stiffness,options.assemblyThreads);
    auto const dofNeighborhood = buildMatrixNeighborhood(mass,stiffness);
    if (options.profile)
      profile.localOperatorExtraction = ProfileTimes::seconds(extractionStart);

    Assembler assembler(spaces);
    Kaskade::SemiImplicitEulerStep<Functional> equation(&F,options.dt);
    CompressionReport compressionReport;

    for (int step = 0; step < steps; ++step)
    {
      double const t = step * options.dt;
      double const stepEnd = std::min(options.finalTime,t+options.dt);
      double const stepDt = stepEnd - t;
      F.time(t);

      Kaskade::RadauTimeGrid grid(options.sdcStartCollocationPoints,t,stepEnd);
      std::vector<StateU> collocationStates(grid.points().N(),Kaskade::component<0>(state));
      double sdcContraction = options.sdcInitialContraction;
      std::vector<double> sweepNorms;
      // BDDC updates complete subdomains; AA may reduce this set after the first sweep.
      std::vector<int> activeIds = allSubdomains(data.localDofs.size());
      int sweep = 0;

      std::cout << "BDDC+SDC step " << step+1 << "/" << steps
                << ", time [" << t << ", " << stepEnd << "]\n";

      for (; sweep < options.maximumSdcSweeps; ++sweep)
      {
        std::vector<std::vector<Vector>> residuals(grid.points().N());
        auto const residualStart = options.profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
        computeBddcSdcResiduals(F,equation,state,collocationStates,grid,data,nDofs,
                                sweep,t,stepDt,assembler,options,activeIds,indexSet,residuals);
        if (options.profile)
          profile.residualAssembly += ProfileTimes::seconds(residualStart);

        auto const updateStart = options.profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
        Kaskade::SDCTimeGrid::RealMatrix Shat;
        if (options.sdcSweepType == 0)
          Kaskade::eulerIntegrationMatrix(grid,Shat);
        else if (options.sdcSweepType == 1)
          Kaskade::luIntegrationMatrix(grid,Shat);
        else
          throw std::runtime_error("sdcSweepType must be 0 (Euler) or 1 (LU)");

        std::vector<std::vector<Vector>> massDifferences(grid.points().N()-1);
        for (auto& intervalData : massDifferences)
          for (auto const& localDofs : data.localDofs)
            intervalData.emplace_back(localDofs.size());
        computeBddcMassDifferences(data,localMassMatrices,collocationStates,massDifferences);
        if (options.profile)
          profile.sdcUpdate += ProfileTimes::seconds(updateStart);

        std::vector<std::vector<Vector>> corrections(grid.points().N());
        for (auto& pointData : corrections)
          for (auto const& localDofs : data.localDofs)
          {
            pointData.emplace_back(localDofs.size());
            pointData.back() = 0.0;
          }

        sweepNorms.push_back(sdcIterationStepBddc<Transfer>(grid,Shat,data,localMassMatrices,localStiffnessMatrices,
                                                            residuals,massDifferences,corrections,activeIds,options,
                                                            options.profile ? &profile : nullptr,
                                                            options.bddcCompressionReport ? &compressionReport : nullptr));

        auto const correctionUpdateStart = options.profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
        for (int i = 1; i < grid.points().N(); ++i)
        {
          Vector globalCorrection = combineSubdomainVector(data,corrections[i],nDofs);
          for (size_t j = 0; j < nDofs; ++j)
            collocationStates[i].coefficients()[j] += globalCorrection[j];
        }
        if (options.profile)
          profile.sdcUpdate += ProfileTimes::seconds(correctionUpdateStart);

        if (sweepNorms.size() > 1)
        {
          double const c = sweepNorms.back()/sweepNorms[sweepNorms.size()-2];
          sdcContraction = std::sqrt(c*sdcContraction);
        }

        std::cout << "  sweep " << sweep+1
                  << ": ||du||=" << sweepNorms.back()
                  << ", contraction=" << sdcContraction
                  << ", active subdomains=" << activeIds.size() << "\n";

        bool const reachedMinimumSweeps = sweep+1 >= options.minimumSdcSweeps;
        bool const smallCorrection = sweepNorms.back() < options.sdcTolerance;
        bool const reliableContraction = sdcContraction < 1.0;
        bool const estimatedSmall = reliableContraction
                                 && sweepNorms.back()*sdcContraction/(1.0-sdcContraction) <= options.sdcAbsoluteTolerance;
        if (reachedMinimumSweeps && (smallCorrection || estimatedSmall))
          break;

        auto const adaptivityStart = options.profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
        std::vector<int> nextActiveIds = selectActiveSubdomains(data,corrections,activeIds,
                                                                dofNeighborhood,sdcContraction,options);
        if (options.profile)
          profile.adaptivity += ProfileTimes::seconds(adaptivityStart);
        if (options.algebraicAdaptivity && options.algebraicAdaptivityTolerance > 0.0)
        {
          if (nextActiveIds.empty())
          {
            std::cout << "  AA selected no active subdomains for the next sweep\n";
            break;
          }
          std::cout << "  AA selected active subdomains for next sweep: "
                    << nextActiveIds.size() << "\n";
        }
        activeIds.swap(nextActiveIds);

        if (grid.points().N() < options.sdcCollocationPoints+1 && sweep+1 < options.maximumSdcSweeps)
        {
          Kaskade::SDCTimeGrid::RealMatrix prolongation;
          grid.refine(prolongation);
          std::vector<StateU> refinedStates(grid.points().N(),collocationStates[0]);
          for (int i = 0; i < static_cast<int>(refinedStates.size()); ++i)
          {
            refinedStates[i] = 0.0;
            for (int j = 0; j < static_cast<int>(collocationStates.size()); ++j)
              refinedStates[i].axpy(prolongation[i][j],collocationStates[j]);
          }
          collocationStates.swap(refinedStates);
        }
      }

      Kaskade::component<0>(state) = collocationStates.back();
      F.time(stepEnd);
    }

    if (options.profile)
      profile.print("SDC + BDDC");
    if (options.bddcCompressionReport)
      compressionReport.print("SDC + BDDC");
    return state;
  }

  template <class Transfer, class Functional, class Spaces, class State, class Element, class Matrix, class Vector, class Options>
  State runBddc(Functional& F,
                Spaces const& spaces,
                State state,
                Element&,
                BddcData<Matrix,Vector> const& data,
                size_t nDofs,
                int steps,
                Options const& options)
  {
    using SemiLinearization = Kaskade::SemiLinearizationAtInner<Kaskade::SemiImplicitEulerStep<Functional>>;
    using Assembler = Kaskade::VariationalFunctionalAssembler<SemiLinearization>;
    using TransmissionScalar = double;
    using BddcSubdomain = Kaskade::BDDC::Subdomain<1,double,double,Transfer>;

    Kaskade::BDDC::InterfaceAverages<1,int> interfaceAverages(data.sharedDofs,
                                                              data.subdomainSizes,
                                                              options.bddcInterfaceTypes);
    auto const interfaceCount = interfaceAverages.interfaceCount();
    std::cout << "BDDC coarse interfaces: corners=" << interfaceCount[0]
              << ", edges=" << interfaceCount[1]
              << ", faces=" << interfaceCount[2] << "\n";

    Assembler assembler(spaces);
    auto zeroState(state);
    auto stepState(state);
    CompressionReport compressionReport;
    zeroState *= 0.0;
    stepState *= 0.0;

    for (int step = 0; step < steps; ++step)
    {
      double const time = step*options.dt;
      F.time(time);
      F.Mass_stiff(10);

      Kaskade::SemiImplicitEulerStep<Functional> equation(&F,options.dt);
      equation.setTau(options.dt);
      assembler.assemble(SemiLinearization(equation,state,state,zeroState),
                         Assembler::RHS,
                         options.assemblyThreads);

      Vector rhs(nDofs);
      assembler.rhs().write(rhs.begin());

      std::vector<Vector> subdomainRhs(data.localDofs.size());
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      {
        subdomainRhs[subdomain] = Vector(data.localDofs[subdomain].size());
        subdomainRhs[subdomain] = 0.0;
        for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        {
          size_t const global = data.localDofs[subdomain][local];
          subdomainRhs[subdomain][local] = data.weights[subdomain][local][0] * rhs[global];
        }
      }

      auto subdomains = constructSubdomains<BddcSubdomain>(data.localMatrices,interfaceAverages,
                                                            options.assemblyThreads);

      if constexpr (requires(Transfer& transfer) {
                      transfer.setQuantizationBits(0);
                      transfer.setRestrictEncoding(true);
                      transfer.setProlongateEncoding(true);
                      transfer.setRestrictTransform(true);
                      transfer.setProlongateTransform(true);
                      transfer.setRestrictBitlengthEncoding(true);
                      transfer.setProlongateBitlengthEncoding(true);
                    })
      {
        for (auto& subdomain : subdomains)
        {
          subdomain.transfer().setQuantizationBits(options.bddcCompressionBits);
          subdomain.transfer().enableTransform(options.bddcGraphLifting
                                                 ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                                                 : Kaskade::BDDC::TransformType::NONE);
          subdomain.transfer().setRestrictEncoding(options.bddcHuffman != 0);
          subdomain.transfer().setProlongateEncoding(options.bddcHuffman != 0);
          subdomain.transfer().setRestrictTransform(options.bddcGraphLifting != 0);
          subdomain.transfer().setProlongateTransform(options.bddcGraphLifting != 0);
          subdomain.transfer().setRestrictBitlengthEncoding(options.bddcBitlength != 0);
          subdomain.transfer().setProlongateBitlengthEncoding(options.bddcBitlength != 0);
        }
      }

      std::vector<int> activeIds(subdomains.size());
      std::iota(activeIds.begin(),activeIds.end(),0);
      bool useCgSolver = options.bddcUseCg != 0;
      bool verbose = options.bddcVerbose != 0;

      Kaskade::BDDC::BDDCSolver<BddcSubdomain> solver(subdomains,
                                                      interfaceAverages.coarseConstraints(),
                                                      activeIds,
                                                      useCgSolver,
                                                      verbose);
      solver.setRhs(subdomainRhs);

      double residual = std::numeric_limits<double>::infinity();
      int iteration = 0;
      for (; iteration < options.bddcIterations; ++iteration)
      {
        residual = solver.solve();
        if (verbose)
          std::cout << "  BDDC iteration " << iteration+1 << ": residual=" << residual << "\n";
        if (residual < options.bddcTolerance)
          break;
      }
      if (options.bddcCompressionReport)
      {
        if (options.bddcCompression)
          compressionReport.add(solver.compressionTraffic());
        else
          compressionReport.addUncompressed(solver.traffic());
      }

      Vector stepVector(nDofs);
      stepVector = 0.0;
      for (size_t subdomain = 0; subdomain < subdomains.size(); ++subdomain)
      {
        auto localSolution = subdomains[subdomain].getSolution();
        for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        {
          size_t const global = data.localDofs[subdomain][local];
          stepVector[global] += data.weights[subdomain][local][0] * localSolution[local];
        }
      }

      for (size_t i = 0; i < nDofs; ++i)
        boost::fusion::at_c<0>(stepState.data).coefficients()[i] = stepVector[i];

      Kaskade::component<0>(state) += Kaskade::component<0>(stepState);
      if (verbose)
        std::cout << "BDDC step " << step+1 << "/" << steps
                  << ": iterations=" << std::min(iteration+1,options.bddcIterations)
                  << ", residual=" << residual << "\n";
    }

    if (options.bddcCompressionReport)
      compressionReport.print("BDDC");
    return state;
  }
}

#endif
