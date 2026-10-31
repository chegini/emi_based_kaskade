#ifndef EMI_BDDC_HH
#define EMI_BDDC_HH

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

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

#include "EMI_bddc_stopping.hh"
#include "EMI_bddc_condition.hh"
#include "EMI_bddc_convergence.hh"
#include "EMI_bddc_precision.hh"

namespace EmiBddc
{
  inline bool mpiRankZero(bool mpiEnabled)
  {
#ifdef KASKADE_HAVE_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (mpiEnabled && initialized)
    {
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      return rank == 0;
    }
#else
    (void)mpiEnabled;
#endif
    return true;
  }

  struct CompressionReport
  {
    std::uint64_t originalSentBytes = 0;
    std::uint64_t originalReceivedBytes = 0;
    std::uint64_t sentBytes = 0;
    std::uint64_t receivedBytes = 0;
    std::uint64_t streamBytes = 0;
    std::uint64_t streamPaddingBytes = 0;
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
      if constexpr (requires {
                      traffic.streamBytes;
                      traffic.streamPaddingBytes;
                    })
      {
        streamBytes += traffic.streamBytes;
        streamPaddingBytes += traffic.streamPaddingBytes;
      }
      encodeTimeMs += traffic.encodeTimeMs;
      decodeTimeMs += traffic.decodeTimeMs;
      encodeCount += traffic.encodeCount;
      decodeCount += traffic.decodeCount;
      codebookBits = std::max(codebookBits, traffic.codebookBits);
    }

    void addUncompressed(std::array<int,6> const& traffic, bool mpiEnabled = false)
    {
      // Domain::traffic() already holds MPI-global logical transfer counts.
      // Contribute them once; printGlobal() reduces rank-local report fields.
      if (!mpiRankZero(mpiEnabled))
        return;
      std::uint64_t const bytes = static_cast<std::uint64_t>(traffic[3])
                                + static_cast<std::uint64_t>(traffic[4]);
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
                << "  logical stream bytes: " << streamBytes << "\n"
                << "  stream padding bytes: " << streamPaddingBytes << "\n"
                << "  raw received bytes: " << originalReceivedBytes << "\n"
                << "  received bytes: " << receivedBytes << "\n";
    }

    // Reduce communication counters across MPI ranks before printing. The
    // shared codebook is common metadata, so report its maximum size once.
    void printGlobal(std::string const& label, bool mpiEnabled, bool show = true) const
    {
#ifdef KASKADE_HAVE_MPI
      int initialized = 0;
      MPI_Initialized(&initialized);
      if (mpiEnabled && initialized)
      {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        unsigned long long localBytes[6] = {
          static_cast<unsigned long long>(originalSentBytes),
          static_cast<unsigned long long>(originalReceivedBytes),
          static_cast<unsigned long long>(sentBytes),
          static_cast<unsigned long long>(receivedBytes),
          static_cast<unsigned long long>(streamBytes),
          static_cast<unsigned long long>(streamPaddingBytes)};
        unsigned long long globalBytes[6] = {0,0,0,0,0,0};
        MPI_Allreduce(localBytes,globalBytes,6,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

        double localTimes[2] = {encodeTimeMs,decodeTimeMs};
        double globalTimes[2] = {0.0,0.0};
        MPI_Allreduce(localTimes,globalTimes,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        unsigned long long localCounts[2] = {
          static_cast<unsigned long long>(encodeCount),
          static_cast<unsigned long long>(decodeCount)};
        unsigned long long globalCounts[2] = {0,0};
        MPI_Allreduce(localCounts,globalCounts,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

        unsigned long long localCodebook = static_cast<unsigned long long>(codebookBits);
        unsigned long long globalCodebook = 0;
        MPI_Allreduce(&localCodebook,&globalCodebook,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);

        if (rank != 0)
          return;

        CompressionReport global;
        global.originalSentBytes = globalBytes[0];
        global.originalReceivedBytes = globalBytes[1];
        global.sentBytes = globalBytes[2];
        global.receivedBytes = globalBytes[3];
        global.streamBytes = globalBytes[4];
        global.streamPaddingBytes = globalBytes[5];
        global.encodeTimeMs = globalTimes[0];
        global.decodeTimeMs = globalTimes[1];
        global.encodeCount = globalCounts[0];
        global.decodeCount = globalCounts[1];
        global.codebookBits = globalCodebook;
        if (show)
          global.print(label);
        return;
      }
#else
      (void)mpiEnabled;
#endif
      if (show)
        print(label);
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
    double correctionAssembly = 0.0;
    double trueResidualChecks = 0.0;
    double sdcUpdate = 0.0;
    double adaptivity = 0.0;
    double localOperatorExtraction = 0.0;
    size_t bddcSystems = 0;
    size_t bddcSolverConstructions = 0;
    size_t bddcSubdomainConstructions = 0;
    size_t bddcSolveCalls = 0;
    size_t trueResidualCheckCalls = 0;

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
                << "  BDDC solver/coarse constructions: " << bddcSolverConstructions << "\n"
                << "  BDDC subdomain constructions: " << bddcSubdomainConstructions << "\n"
                << "  BDDC solve calls: " << bddcSolveCalls
                << " across " << bddcSystems << " systems (wall seconds: "
                << bddcSolve << ", average calls/system: "
                << (bddcSystems > 0 ? static_cast<double>(bddcSolveCalls)/bddcSystems : 0.0)
                << ")\n"
                << "  correction assembly/reduction: " << correctionAssembly << "\n"
                << "  true global residual checks: " << trueResidualCheckCalls
                << " (wall seconds, excluding correction assembly: " << trueResidualChecks << ")\n"
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
    // Optional labels per sharedDofs entry; empty retains owner-only coarse groups.
    std::vector<int> interfaceComponents;
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
                                        int nTasks,
                                        bool mpiEnabled, bool detailedOutput = true)
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

    if (detailedOutput && mpiRankZero(mpiEnabled))
      std::cout << "BDDC subdomains: " << data.tags.size()
                << ", shared dofs: " << data.sharedDofs.size() << "\n";
    return data;
  }

  // Match the legacy EMI BDDC local spaces: keep every DOF supported by the
  // cells of a material tag and add the coordinate-matched DOFs on adjacent
  // cells of a different tag.  The legacy mapper has distinct global IDs on
  // the two sides of a material interface, so matrix-support membership alone
  // does not recover these interface-neighbor copies.
  template <class Grid, class Space, class Material, class Matrix, class Vector>
  void expandSubdomainsToLegacyInterfaceDofs(BddcData<Matrix,Vector>& data,
                                              Grid const& grid,
                                              Space& space,
                                              Material const& material,
                                              bool mpiEnabled, bool detailedOutput = true)
  {
    using LocalDof = Kaskade::BDDC::LocalDof;
    Dune::FieldVector<double,Grid::dimension> zero(0.0);
    std::map<int,size_t> tagToSubdomain;
    for (size_t subdomain = 0; subdomain < data.tags.size(); ++subdomain)
      tagToSubdomain[data.tags[subdomain]] = subdomain;

    auto localNodes = [&space](auto const& cell)
    {
      typename Space::Evaluator evaluator(space);
      evaluator.moveTo(cell);
      auto const& interpolationNodes = evaluator.shapeFunctions().interpolationNodes();
      auto const globalDofs = space.mapper().globalIndices(cell);
      if (interpolationNodes.size() != globalDofs.size())
        throw std::runtime_error("BDDC interface mapper and interpolation nodes have different sizes.");

      using Coordinate = Dune::FieldVector<double,Grid::dimension>;
      std::vector<std::pair<size_t,Coordinate>> nodes;
      nodes.reserve(globalDofs.size());
      for (size_t i = 0; i < globalDofs.size(); ++i)
      {
        auto const x = cell.geometry().global(interpolationNodes[i]);
        Coordinate coordinate;
        for (int d = 0; d < Grid::dimension; ++d)
          coordinate[d] = x[d];
        nodes.emplace_back(globalDofs[i],coordinate);
      }
      return nodes;
    };

    std::vector<std::set<size_t>> expanded(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      expanded[subdomain].insert(data.localDofs[subdomain].begin(),
                                 data.localDofs[subdomain].end());

    for (auto const& cell : Dune::elements(grid.leafGridView()))
    {
      int const insideTag = material.value(cell,zero);
      auto inside = tagToSubdomain.find(insideTag);
      if (inside == tagToSubdomain.end())
        continue;
      auto const insideNodes = localNodes(cell);
      for (auto const& intersection : Dune::intersections(grid.leafGridView(),cell))
      {
        if (!intersection.neighbor())
          continue;
        auto const& outsideCell = intersection.outside();
        int const outsideTag = material.value(outsideCell,zero);
        if (insideTag == outsideTag)
          continue;
        auto outside = tagToSubdomain.find(outsideTag);
        if (outside == tagToSubdomain.end())
          continue;
        auto const outsideNodes = localNodes(outsideCell);
        for (auto const& insideNode : insideNodes)
          for (auto const& outsideNode : outsideNodes)
          {
            bool sameCoordinate = true;
            for (int d = 0; d < Grid::dimension; ++d)
              sameCoordinate = sameCoordinate &&
                insideNode.second[d] == outsideNode.second[d];
            if (!sameCoordinate)
              continue;
            expanded[inside->second].insert(outsideNode.first);
            expanded[outside->second].insert(insideNode.first);
          }
      }
    }

    data.sharedDofs.clear();
    data.dofSubdomains.assign(data.dofCells.size(),{});
    data.subdomainSizes.clear();
    std::vector<std::vector<LocalDof>> occurrences(data.dofCells.size());
    for (size_t subdomain = 0; subdomain < expanded.size(); ++subdomain)
    {
      data.localDofs[subdomain].assign(expanded[subdomain].begin(),
                                       expanded[subdomain].end());
      if (data.localDofs[subdomain].size() >
          static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("BDDC local subdomain is too large for int local indices.");
      data.subdomainSizes.push_back(static_cast<int>(data.localDofs[subdomain].size()));
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        if (global >= occurrences.size())
          throw std::runtime_error("BDDC interface DOF is outside the global vector.");
        occurrences[global].push_back(
          {static_cast<int>(subdomain),static_cast<int>(local)});
      }
    }

    for (size_t global = 0; global < occurrences.size(); ++global)
    {
      if (occurrences[global].size() > 1)
        data.sharedDofs.push_back(occurrences[global]);
      for (auto const& occurrence : occurrences[global])
        data.dofSubdomains[global].push_back(occurrence.s);
    }

    if (data.sharedDofs.empty() && data.tags.size() > 1)
      throw std::runtime_error(
        "BDDC interface expansion found multiple subdomains but no shared dofs.");

    if (detailedOutput && mpiRankZero(mpiEnabled))
    {
      size_t totalDofs = 0;
      for (auto const& localDofs : data.localDofs)
        totalDofs += localDofs.size();
      std::cout << "BDDC legacy interface-neighbor local spaces: total local dofs="
                << totalDofs << ", shared dofs=" << data.sharedDofs.size() << "\n";
    }
  }

  // Add one algebraic overlap layer to the material-tag subdomains. The
  // original tag partition is useful as a starting point, but a local
  // principal matrix must contain both endpoints of every nonzero coupling
  // that it is expected to represent. This is essential for membrane and
  // gap-junction entries, whose two DOFs can belong to different tags.
  template <class Matrix, class Vector>
  void expandSubdomainsToMatrixCouplings(BddcData<Matrix,Vector>& data,
                                         Matrix const& globalMatrix,
                                         int nTasks,
                                         bool mpiEnabled)
  {
    using Kaskade::BDDC::LocalDof;

    std::vector<std::set<size_t>> expanded(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      expanded[subdomain].insert(data.localDofs[subdomain].begin(),
                                 data.localDofs[subdomain].end());

    for (auto row = globalMatrix.begin(); row != globalMatrix.end(); ++row)
    {
      size_t const rowIndex = row.index();
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        if ((*col)[0][0] == 0.0)
          continue;
        size_t const colIndex = col.index();
        if (rowIndex >= data.dofSubdomains.size() ||
            colIndex >= data.dofSubdomains.size())
          continue;

        for (int subdomain : data.dofSubdomains[rowIndex])
          expanded[static_cast<size_t>(subdomain)].insert(colIndex);
        for (int subdomain : data.dofSubdomains[colIndex])
          expanded[static_cast<size_t>(subdomain)].insert(rowIndex);
      }
    }

    data.sharedDofs.clear();
    data.dofSubdomains.assign(data.dofCells.size(),{});
    data.subdomainSizes.clear();
    std::vector<std::vector<LocalDof>> occurrences(data.dofCells.size());

    for (size_t subdomain = 0; subdomain < expanded.size(); ++subdomain)
    {
      data.localDofs[subdomain].assign(expanded[subdomain].begin(),
                                       expanded[subdomain].end());
      data.subdomainSizes.push_back(static_cast<int>(data.localDofs[subdomain].size()));
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        occurrences[data.localDofs[subdomain][local]].push_back(
          {static_cast<int>(subdomain),static_cast<int>(local)});
    }

    for (size_t global = 0; global < occurrences.size(); ++global)
    {
      if (occurrences[global].size() > 1)
        data.sharedDofs.push_back(occurrences[global]);
      for (auto const& occurrence : occurrences[global])
        data.dofSubdomains[global].push_back(occurrence.s);
    }

    data.weights.resize(data.localDofs.size());
    size_t const taskLimit = nTasks > 0 ? static_cast<size_t>(nTasks)
                                        : std::numeric_limits<size_t>::max();
    auto buildWeights = [&](size_t subdomain)
    {
      Vector weight(data.localDofs[subdomain].size());
      weight = 1.0;
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        size_t const ownerCount = std::max<size_t>(1,occurrences[global].size());
        weight[local][0] = 1.0/static_cast<double>(ownerCount);
      }
      data.weights[subdomain] = std::move(weight);
    };
    if (taskLimit > 1)
      Kaskade::parallelFor(0,data.localDofs.size(),buildWeights,taskLimit);
    else
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        buildWeights(subdomain);

    if (mpiRankZero(mpiEnabled))
    {
      size_t totalDofs = 0;
      for (auto const& localDofs : data.localDofs)
        totalDofs += localDofs.size();
      std::cout << "BDDC algebraic overlap: total local dofs=" << totalDofs
                << ", shared dofs=" << data.sharedDofs.size() << "\n";
    }
  }

  // Add only neighboring-cell DOFs that actually couple across a material
  // interface. This retains membrane/gap couplings without adding unrelated
  // cell-interior DOFs to the local spaces.
  template <class Grid, class Space, class Material, class Matrix, class Vector>
  void expandSubdomainsToInterfaceNeighbors(BddcData<Matrix,Vector>& data,
                                            Grid const& grid,
                                            Space const& space,
                                            Material const& material,
                                            Matrix const& globalMatrix,
                                            int nTasks,
                                            bool mpiEnabled, bool detailedOutput = true)
  {
    using LocalDof = Kaskade::BDDC::LocalDof;
    Dune::FieldVector<double,Grid::dimension> zero(0.0);
    std::map<int,size_t> tagToSubdomain;
    for (size_t subdomain = 0; subdomain < data.tags.size(); ++subdomain)
      tagToSubdomain[data.tags[subdomain]] = subdomain;

    std::vector<std::set<size_t>> expanded(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      expanded[subdomain].insert(data.localDofs[subdomain].begin(),
                                 data.localDofs[subdomain].end());

    for (auto const& cell : Dune::elements(grid.leafGridView()))
    {
      int const insideTag = material.value(cell,zero);
      auto inside = tagToSubdomain.find(insideTag);
      if (inside == tagToSubdomain.end())
        continue;
      for (auto const& intersection : Dune::intersections(grid.leafGridView(),cell))
      {
        if (!intersection.neighbor())
          continue;
        auto const& outsideCell = intersection.outside();
        int const outsideTag = material.value(outsideCell,zero);
        if (insideTag == outsideTag)
          continue;
        auto outside = tagToSubdomain.find(outsideTag);
        if (outside == tagToSubdomain.end())
          continue;
        auto const insideDofs = space.mapper().globalIndices(cell);
        auto const outsideDofs = space.mapper().globalIndices(outsideCell);
        auto hasCoupling = [&globalMatrix](size_t rowIndex, size_t columnIndex)
        {
          auto row = globalMatrix[rowIndex];
          for (auto col = row.begin(); col != row.end(); ++col)
            if (col.index() == columnIndex)
              return (*col)[0][0] != 0.0;
          return false;
        };
        for (size_t insideDof : insideDofs)
          for (size_t outsideDof : outsideDofs)
            if (hasCoupling(insideDof,outsideDof) || hasCoupling(outsideDof,insideDof))
            {
              expanded[inside->second].insert(outsideDof);
              expanded[outside->second].insert(insideDof);
            }
      }
    }

    data.sharedDofs.clear();
    data.dofSubdomains.assign(data.dofCells.size(),{});
    data.subdomainSizes.clear();
    std::vector<std::vector<LocalDof>> occurrences(data.dofCells.size());
    for (size_t subdomain = 0; subdomain < expanded.size(); ++subdomain)
    {
      data.localDofs[subdomain].assign(expanded[subdomain].begin(),
                                       expanded[subdomain].end());
      data.subdomainSizes.push_back(static_cast<int>(data.localDofs[subdomain].size()));
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        occurrences[data.localDofs[subdomain][local]].push_back(
          {static_cast<int>(subdomain),static_cast<int>(local)});
    }
    for (size_t global = 0; global < occurrences.size(); ++global)
    {
      if (occurrences[global].size() > 1)
        data.sharedDofs.push_back(occurrences[global]);
      for (auto const& occurrence : occurrences[global])
        data.dofSubdomains[global].push_back(occurrence.s);
    }

    data.weights.resize(data.localDofs.size());
    size_t const taskLimit = nTasks > 0 ? static_cast<size_t>(nTasks)
                                        : std::numeric_limits<size_t>::max();
    auto buildWeights = [&](size_t subdomain)
    {
      Vector weight(data.localDofs[subdomain].size());
      weight = 1.0;
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        size_t const count = std::max<size_t>(1,occurrences[global].size());
        weight[local][0] = 1.0/static_cast<double>(count);
      }
      data.weights[subdomain] = std::move(weight);
    };
    if (taskLimit > 1)
      Kaskade::parallelFor(0,data.localDofs.size(),buildWeights,taskLimit);
    else
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        buildWeights(subdomain);

    if (detailedOutput && mpiRankZero(mpiEnabled))
    {
      size_t totalDofs = 0;
      for (auto const& localDofs : data.localDofs)
        totalDofs += localDofs.size();
      std::cout << "BDDC interface-neighbor overlap: total local dofs="
                << totalDofs << ", shared dofs=" << data.sharedDofs.size() << "\n";
    }
  }

  template <class Matrix, class Vector>
  void buildDiagonalWeights(BddcData<Matrix,Vector>& data,
                            std::vector<Matrix> const& localMatrices,
                            Matrix const& globalMatrix,
                            bool mpiEnabled, bool detailedOutput = true)
  {
    bool const outputRankZero = mpiRankZero(mpiEnabled);
    for (size_t subdomain = 0; subdomain < localMatrices.size(); ++subdomain)
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        double const globalDiagonal = globalMatrix[global][global][0][0];
        double const localDiagonal = localMatrices[subdomain][local][local][0][0];
        if (std::abs(globalDiagonal) > 1e-14)
          data.weights[subdomain][local][0] = localDiagonal/globalDiagonal;
      }
    if (detailedOutput && outputRankZero)
      std::cout << "BDDC weights: diagonal-based local/global contributions\n";
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

  // Build additive local operators from disjoint material contributions.
  // Their embedded sum must reproduce the global matrix exactly.
  template <class Matrix, class Vector>
  void setLocalMatricesFromContributions(BddcData<Matrix,Vector>& data,
                                         std::vector<Matrix> const& contributions,
                                         bool mpiEnabled, bool detailedOutput = true)
  {
    using Kaskade::BDDC::LocalDof;
    if (contributions.size() != data.tags.size())
      throw std::runtime_error("BDDC local contribution count does not match material subdomains.");
    if (data.localDofs.size() != contributions.size())
      throw std::runtime_error("BDDC local-space count does not match material subdomains.");
    std::vector<Matrix> localMatrices(contributions.size());
    for (size_t subdomain = 0; subdomain < contributions.size(); ++subdomain)
    {
      if (data.localDofs[subdomain].empty())
        throw std::runtime_error("BDDC material subdomain has an empty local space.");
      std::set<size_t> localDofSet(data.localDofs[subdomain].begin(),
                                   data.localDofs[subdomain].end());
      auto const& contribution = contributions[subdomain];
      for (auto row = contribution.begin(); row != contribution.end(); ++row)
        for (auto col = row.begin(); col != row.end(); ++col)
          if ((*col)[0][0] != 0.0)
          {
            if (!localDofSet.contains(row.index()) ||
                !localDofSet.contains(col.index()))
              throw std::runtime_error("legacy BDDC local space omits a nonzero material coupling endpoint.");
          }
      bool hasNonzero = false;
      for (auto row = contribution.begin(); row != contribution.end() && !hasNonzero; ++row)
        for (auto col = row.begin(); col != row.end(); ++col)
          if ((*col)[0][0] != 0.0)
          {
            hasNonzero = true;
            break;
          }
      if (!hasNonzero)
        throw std::runtime_error("BDDC material subdomain has an empty assembled operator contribution.");
      localMatrices[subdomain] = Matrix(data.localDofs[subdomain],contribution);
    }

    std::vector<std::vector<LocalDof>> occurrences(data.dofSubdomains.size());
    data.subdomainSizes.resize(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
    {
      if (data.localDofs[subdomain].size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("BDDC local subdomain is too large for int local indices.");
      data.subdomainSizes[subdomain] = static_cast<int>(data.localDofs[subdomain].size());
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        occurrences[data.localDofs[subdomain][local]].push_back(
          {static_cast<int>(subdomain),static_cast<int>(local)});
    }

    data.sharedDofs.clear();
    data.dofSubdomains.assign(occurrences.size(),{});
    for (size_t global = 0; global < occurrences.size(); ++global)
    {
      if (occurrences[global].size() > 1)
        data.sharedDofs.push_back(occurrences[global]);
      for (auto const& occurrence : occurrences[global])
        data.dofSubdomains[global].push_back(occurrence.s);
    }

    data.localMatrices = std::move(localMatrices);
    data.weights.resize(data.localDofs.size());
    size_t totalLocalDofs = 0;
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
    {
      totalLocalDofs += data.localDofs[subdomain].size();
      Vector weight(data.localDofs[subdomain].size());
      weight = 1.0;
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        size_t const copies = std::max<size_t>(1,occurrences[global].size());
        weight[local][0] = 1.0/static_cast<double>(copies);
      }
      data.weights[subdomain] = std::move(weight);
    }

    if (detailedOutput && mpiRankZero(mpiEnabled))
      std::cout << "BDDC legacy-interface local spaces: total local dofs="
                << totalLocalDofs << ", shared global dofs="
                << data.sharedDofs.size() << "\n";
  }

  template <class Matrix, class Vector>
  void dumpLocalBddcData(BddcData<Matrix,Vector> const& data,
                         std::string const& directory,
                         Vector const& globalRhs,
                         int interfaceTypes)
  {
    std::filesystem::create_directories(directory);
    std::ofstream matrices(directory + "/local_matrices.csv");
    std::ofstream dofsAndWeights(directory + "/local_dofs_weights.csv");
    std::ofstream localRhs(directory + "/local_rhs.csv");
    std::ofstream sharedDofs(directory + "/shared_dofs.csv");
    std::ofstream coarseConstraints(directory + "/coarse_constraints.csv");
    if (!matrices || !dofsAndWeights || !localRhs || !sharedDofs || !coarseConstraints)
      throw std::runtime_error("cannot open BDDC local-data diagnostic files in " + directory);

    matrices << "subdomain,tag,row_global,col_global,value\n" << std::setprecision(17);
    dofsAndWeights << "subdomain,tag,local_dof,global_dof,weight\n" << std::setprecision(17);
    localRhs << "subdomain,tag,local_dof,global_dof,weight,value\n" << std::setprecision(17);
    sharedDofs << "global_dof,subdomain,tag,local_dof,weight\n" << std::setprecision(17);
    coarseConstraints << "constraint_id,type,subdomain,tag,global_dof,local_dof\n";
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
    {
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        dofsAndWeights << subdomain << ',' << data.tags[subdomain] << ',' << local << ','
                       << data.localDofs[subdomain][local] << ','
                       << data.weights[subdomain][local][0] << '\n';
        localRhs << subdomain << ',' << data.tags[subdomain] << ',' << local << ','
                 << data.localDofs[subdomain][local] << ','
                 << data.weights[subdomain][local][0] << ','
                 << data.weights[subdomain][local][0]
                    * globalRhs[data.localDofs[subdomain][local]][0] << '\n';
      }

      auto const& matrix = data.localMatrices[subdomain];
      for (auto row = matrix.begin(); row != matrix.end(); ++row)
        for (auto col = row.begin(); col != row.end(); ++col)
          if ((*col)[0][0] != 0.0)
            matrices << subdomain << ',' << data.tags[subdomain] << ','
                     << data.localDofs[subdomain][row.index()] << ','
                     << data.localDofs[subdomain][col.index()] << ','
                     << (*col)[0][0] << '\n';
    }

    for (auto const& occurrences : data.sharedDofs)
    {
      if (occurrences.empty())
        continue;
      size_t const global = data.localDofs[static_cast<size_t>(occurrences.front().s)]
                                           [static_cast<size_t>(occurrences.front().i)];
      for (auto const& occurrence : occurrences)
      {
        size_t const subdomain = static_cast<size_t>(occurrence.s);
        size_t const local = static_cast<size_t>(occurrence.i);
        sharedDofs << global << ',' << occurrence.s << ',' << data.tags[subdomain] << ','
                   << occurrence.i << ',' << data.weights[subdomain][local][0] << '\n';
      }
    }

    Kaskade::BDDC::InterfaceAverages<1,int> interfaceAverages(
      data.sharedDofs,data.subdomainSizes,interfaceTypes,data.interfaceComponents);
    auto const& actualConstraints = interfaceAverages.coarseConstraints();
    auto const actualCounts = interfaceAverages.interfaceCount();
    // Read the actual C rows: several constraints can share the same owners
    // when disconnected components have separate averages.
    std::vector<std::vector<size_t>> rowsBySubdomain(data.tags.size());
    for (size_t id = 0; id < actualConstraints.size(); ++id)
      for (int subdomain : actualConstraints[id].first)
        rowsBySubdomain.at(subdomain).push_back(id);
    for (size_t subdomain = 0; subdomain < data.tags.size(); ++subdomain)
    {
      auto const C = interfaceAverages.template coarseConstraint<double>(subdomain);
      if (C.N() != rowsBySubdomain[subdomain].size())
        throw std::runtime_error("BDDC dump: constraint row mapping mismatch");
      for (auto row = C.begin(); row != C.end(); ++row)
      {
        size_t const id = rowsBySubdomain[subdomain][row.index()];
        size_t count = 0;
        for (auto col = row->begin(); col != row->end(); ++col)
          ++count;
        char const* type = actualConstraints[id].first.size() == 2 ? "face"
                         : count == 1 ? "corner" : "edge";
        for (auto col = row->begin(); col != row->end(); ++col)
          coarseConstraints << id << ',' << type << ',' << subdomain << ','
                            << data.tags[subdomain] << ','
                            << data.localDofs[subdomain][col.index()] << ','
                            << col.index() << '\n';
      }
    }

    std::ofstream summary(directory + "/coarse_summary.txt");
    if (!summary)
      throw std::runtime_error("cannot open BDDC coarse-space summary in " + directory);
    summary << "shared_global_dofs," << data.sharedDofs.size() << '\n'
            << "enabled_constraints," << actualConstraints.size() << '\n'
            << "corners," << actualCounts[0] << '\n'
            << "edges," << actualCounts[1] << '\n'
            << "faces," << actualCounts[2] << '\n'
            << "interface_type_mask," << interfaceTypes << '\n';
  }

  // Assign every global nonzero matrix entry to exactly one subdomain that
  // contains both endpoints. Principal extraction duplicates entries when
  // endpoints occur in several overlapping subdomains; this removes those
  // duplicate copies so the embedded local matrices add up to the global
  // matrix exactly. The local operators still require a solver test.
  template <class Matrix, class Vector>
  void assignLocalMatrixEntryOwners(BddcData<Matrix,Vector> const& data,
                                    std::vector<Matrix>& localMatrices,
                                    Matrix const& globalMatrix,
                                    bool mpiEnabled)
  {
    bool const outputRankZero = mpiRankZero(mpiEnabled);

    using Key = std::pair<size_t,size_t>;
    std::map<Key,size_t> owners;
    size_t unowned = 0;
    for (auto row = globalMatrix.begin(); row != globalMatrix.end(); ++row)
    {
      size_t const globalRow = row.index();
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        if ((*col)[0][0] == 0.0)
          continue;
        size_t const globalCol = col.index();
        size_t owner = data.localDofs.size();
        if (globalRow < data.dofSubdomains.size() &&
            globalCol < data.dofSubdomains.size())
          for (int candidate : data.dofSubdomains[globalRow])
            if (std::find(data.dofSubdomains[globalCol].begin(),
                          data.dofSubdomains[globalCol].end(),candidate) !=
                data.dofSubdomains[globalCol].end())
              owner = std::min(owner,static_cast<size_t>(candidate));
        if (owner == data.localDofs.size())
          ++unowned;
        else
          owners[{globalRow,globalCol}] = owner;
      }
    }

    size_t removedDuplicates = 0;
    for (size_t subdomain = 0; subdomain < localMatrices.size(); ++subdomain)
    {
      auto& localMatrix = localMatrices[subdomain];
      for (auto row = localMatrix.begin(); row != localMatrix.end(); ++row)
      {
        size_t const localRow = row.index();
        size_t const globalRow = data.localDofs[subdomain][localRow];
        for (auto col = row.begin(); col != row.end(); ++col)
        {
          if ((*col)[0][0] == 0.0)
            continue;
          size_t const localCol = col.index();
          size_t const globalCol = data.localDofs[subdomain][localCol];
          auto const owner = owners.find({globalRow,globalCol});
          if (owner == owners.end() || owner->second != subdomain)
          {
            (*col)[0][0] = 0.0;
            ++removedDuplicates;
          }
        }
      }
    }

    if (outputRankZero)
      std::cout << "BDDC matrix entry ownership: owned=" << owners.size()
                << " unowned=" << unowned
                << " removedDuplicateCopies=" << removedDuplicates << "\n";
  }

  // Reconstruct the sum of local matrices in global coordinates. This is
  // intended for small debugging meshes; avoid the sparse-map overhead for
  // production-size systems.
  template <class Matrix, class Vector>
  void printLocalMatrixReconstruction(BddcData<Matrix,Vector> const& data,
                                      std::vector<Matrix> const& localMatrices,
                                      Matrix const& globalMatrix,
                                      size_t nDofs,
                                      bool mpiEnabled)
  {
    if (!mpiRankZero(mpiEnabled))
      return;
    if (nDofs > 5000)
    {
      std::cout << "BDDC local matrix reconstruction: skipped for N=" << nDofs
                << " (diagnostic limit is 5000)\n";
      return;
    }

    using Key = std::pair<size_t,size_t>;
    struct Entry { double value = 0.0; size_t count = 0; };
    std::map<Key,Entry> reconstructed;
    size_t localNonzeros = 0;
    size_t maxMultiplicity = 0;

    for (size_t subdomain = 0; subdomain < localMatrices.size(); ++subdomain)
    {
      auto const& localMatrix = localMatrices[subdomain];
      for (auto row = localMatrix.begin(); row != localMatrix.end(); ++row)
      {
        size_t const localRow = row.index();
        if (localRow >= data.localDofs[subdomain].size())
          continue;
        size_t const globalRow = data.localDofs[subdomain][localRow];
        for (auto col = row.begin(); col != row.end(); ++col)
        {
          double const value = (*col)[0][0];
          if (value == 0.0)
            continue;
          size_t const localCol = col.index();
          if (localCol >= data.localDofs[subdomain].size())
            continue;
          size_t const globalCol = data.localDofs[subdomain][localCol];
          auto& entry = reconstructed[{globalRow,globalCol}];
          entry.value += value;
          ++entry.count;
          ++localNonzeros;
          maxMultiplicity = std::max(maxMultiplicity,entry.count);
        }
      }
    }

    size_t globalNonzeros = 0;
    size_t missing = 0;
    size_t duplicated = 0;
    std::map<Key,double> globalEntries;
    double maxDifference = 0.0;
    double squaredDifference = 0.0;
    for (auto row = globalMatrix.begin(); row != globalMatrix.end(); ++row)
    {
      size_t const globalRow = row.index();
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const reference = (*col)[0][0];
        if (reference == 0.0)
          continue;
        ++globalNonzeros;
        size_t const globalCol = col.index();
        globalEntries[{globalRow,globalCol}] = reference;
        auto const it = reconstructed.find({globalRow,globalCol});
        double const reconstructedValue = it == reconstructed.end()
                                            ? 0.0 : it->second.value;
        if (it == reconstructed.end())
          ++missing;
        if (it != reconstructed.end() && it->second.count > 1)
          ++duplicated;
        double const difference = reconstructedValue-reference;
        maxDifference = std::max(maxDifference,std::abs(difference));
        squaredDifference += difference*difference;
      }
    }

    size_t extra = 0;
    for (auto const& [key,entry] : reconstructed)
      if (globalEntries.find(key) == globalEntries.end())
        ++extra;

    std::cout << "BDDC local matrix reconstruction: "
              << "globalNonzeros=" << globalNonzeros
              << " localNonzeros=" << localNonzeros
              << " reconstructedNonzeros=" << reconstructed.size()
              << " missing=" << missing
              << " duplicated=" << duplicated
              << " extra=" << extra
              << " maxMultiplicity=" << maxMultiplicity
              << " maxAbsDifference=" << maxDifference
              << " l2Difference=" << std::sqrt(squaredDifference)
              << "\n";
  }

  template <class Matrix, class Vector>
  void printLocalRhsReconstruction(BddcData<Matrix,Vector> const& data,
                                   std::vector<Vector> const& subdomainRhs,
                                   Vector const& rhs,
                                   bool mpiEnabled)
  {
    if (!mpiRankZero(mpiEnabled))
      return;
    std::vector<double> reconstructed(rhs.size(),0.0);
    for (size_t subdomain = 0; subdomain < subdomainRhs.size(); ++subdomain)
      for (size_t local = 0; local < subdomainRhs[subdomain].size(); ++local)
        reconstructed[data.localDofs[subdomain][local]] += subdomainRhs[subdomain][local];

    double maxDifference = 0.0;
    double squaredDifference = 0.0;
    for (size_t global = 0; global < rhs.size(); ++global)
    {
      double const difference = reconstructed[global]-rhs[global];
      maxDifference = std::max(maxDifference,std::abs(difference));
      squaredDifference += difference*difference;
    }
    std::cout << "BDDC local RHS reconstruction: maxAbsDifference="
              << maxDifference << " l2Difference="
              << std::sqrt(squaredDifference) << "\n";
  }

  // Check whether every nonzero global matrix coupling is represented in at
  // least one local BDDC matrix. A principal submatrix extraction preserves
  // A_ij only when row i and column j occur in a common subdomain.
  template <class Matrix, class Vector>
  void printUncoveredMatrixCouplings(BddcData<Matrix,Vector> const& data,
                                     Matrix const& globalMatrix,
                                     char const* label,
                                     bool mpiEnabled, bool detailedOutput = true)
  {
    if (!mpiRankZero(mpiEnabled))
      return;

    size_t nonzeros = 0;
    size_t covered = 0;
    size_t uncovered = 0;
    double maxUncoveredAbs = 0.0;
    size_t firstRow = 0;
    size_t firstCol = 0;
    bool haveFirst = false;

    for (auto row = globalMatrix.begin(); row != globalMatrix.end(); ++row)
    {
      size_t const rowIndex = row.index();
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const value = (*col)[0][0];
        if (value == 0.0)
          continue;
        ++nonzeros;

        size_t const colIndex = col.index();
        bool commonSubdomain = false;
        if (rowIndex < data.dofSubdomains.size() &&
            colIndex < data.dofSubdomains.size())
          for (int rowSubdomain : data.dofSubdomains[rowIndex])
            if (std::find(data.dofSubdomains[colIndex].begin(),
                          data.dofSubdomains[colIndex].end(),
                          rowSubdomain) != data.dofSubdomains[colIndex].end())
            {
              commonSubdomain = true;
              break;
            }

        if (commonSubdomain)
          ++covered;
        else
        {
          ++uncovered;
          if (std::abs(value) > maxUncoveredAbs)
          {
            maxUncoveredAbs = std::abs(value);
            if (!haveFirst)
            {
              firstRow = rowIndex;
              firstCol = colIndex;
              haveFirst = true;
            }
          }
        }
      }
    }

    // Always report missing couplings, even in quiet mode.
    if (!detailedOutput && uncovered == 0)
      return;
    std::cout << "BDDC coupling coverage (" << label << "): "
              << "nonzeros=" << nonzeros
              << " covered=" << covered
              << " uncovered=" << uncovered
              << " maxUncoveredAbs=" << maxUncoveredAbs;
    if (haveFirst)
      std::cout << " firstUncovered=(" << firstRow << "," << firstCol << ")";
    std::cout << "\n";
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

  // AA decisions are made from rank-local corrections. In owner-only MPI mode,
  // combine them before the next sweep so every rank enters the same solver
  // collectives with the same global active-subdomain set.
  inline std::vector<int> synchronizeActiveSubdomains(std::vector<int> activeIds,
                                                       size_t nSubdomains,
                                                       bool mpiEnabled)
  {
#ifdef KASKADE_HAVE_MPI
    if (mpiEnabled)
    {
      std::vector<int> localMask(nSubdomains,0);
      std::vector<int> globalMask(nSubdomains,0);
      for (int id : activeIds)
        if (id >= 0 && static_cast<size_t>(id) < nSubdomains)
          localMask[static_cast<size_t>(id)] = 1;

      MPI_Allreduce(localMask.data(),globalMask.data(),static_cast<int>(nSubdomains),
                    MPI_INT,MPI_MAX,MPI_COMM_WORLD);

      activeIds.clear();
      for (size_t id = 0; id < nSubdomains; ++id)
        if (globalMask[id] != 0)
          activeIds.push_back(static_cast<int>(id));
    }
#else
    (void)nSubdomains;
    (void)mpiEnabled;
#endif
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
                    candidate.setRestrictZigzagEncoding(true);
                    candidate.setProlongateZigzagEncoding(true);
                  })
    {
      configureTransferPrecision(transfer,options);
      if constexpr (requires(Transfer& candidate) { candidate.setCompressionMinInterfaceSize(0); }) {
        transfer.setCompressionMinInterfaceSize(static_cast<std::size_t>(std::max(0, options.bddcCompressionMinInterfaceSize)));
      }
      if constexpr (requires(Transfer& candidate) { candidate.setCompressionThresholdInterfaceTypes(0); }) {
        transfer.setCompressionThresholdInterfaceTypes(options.bddcCompressionThresholdInterfaceTypes);
      }
      transfer.enableTransform((options.bddcGraphLifting || options.bddcGraphLiftingAfterQuantization)
                                 ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                                 : Kaskade::BDDC::TransformType::NONE);
      transfer.setRestrictEncoding(options.bddcHuffman != 0);
      transfer.setProlongateEncoding(options.bddcHuffman != 0);
      transfer.setRestrictTransform(options.bddcGraphLifting != 0);
      transfer.setProlongateTransform(options.bddcGraphLifting != 0);
      transfer.setRestrictTransformAfterQuantization(options.bddcGraphLiftingAfterQuantization != 0);
      transfer.setProlongateTransformAfterQuantization(options.bddcGraphLiftingAfterQuantization != 0);
      transfer.setRestrictBitlengthEncoding(options.bddcBitlength != 0);
      transfer.setProlongateBitlengthEncoding(options.bddcBitlength != 0);
      transfer.setRestrictZigzagEncoding(options.bddcZigzag != 0);
      transfer.setProlongateZigzagEncoding(options.bddcZigzag != 0);
    }
  }

  // Construct only the subdomains owned by this rank. In serial mode the
  // ownership layout contains every subdomain, so the same path also covers
  // non-MPI BDDC without maintaining a replicated MPI implementation.
  template <class Subdomain, class Matrix, class Interfaces>
  Kaskade::BDDC::OwnedSubdomainStorage<Subdomain>
  constructOwnedSubdomains(std::vector<Matrix> const& localOperators,
                           Interfaces const& interfaces,
                           int nThreads,
                           bool mpiEnabled,
                           bool masterWorker=false)
  {
    Kaskade::BDDC::MpiDomainLayout layout;
    layout.configure(mpiEnabled, localOperators.size(), masterWorker);
    Kaskade::BDDC::OwnedSubdomainStorage<Subdomain> subdomains(localOperators.size());
    std::exception_ptr constructionError;
    std::mutex constructionErrorMutex;
    auto construct = [&](size_t subdomain)
    {
      if (!layout.owns(static_cast<int>(subdomain)))
        return;
      try
      {
        subdomains.emplace(subdomain, static_cast<int>(subdomain),
                           localOperators[subdomain], interfaces);
        if constexpr (requires(Subdomain& subdomainObject, Matrix const& matrix) {
                        subdomainObject.transfer().setGraphLiftingFromMatrix(matrix);
                      })
          subdomains[subdomain].transfer().setGraphLiftingFromMatrix(localOperators[subdomain]);
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
                                                              options.bddcInterfaceTypes,
                                                              data.interfaceComponents);
    if (profile)
      profile->interfaceSetup += ProfileTimes::seconds(interfaceSetupStart);

    std::vector<int> solverActiveIds(activeIds);
    bool useCgSolver = options.bddcUseCg != 0;
    // Explicit iteration verbosity also applies in quiet mode.
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
      bool const mpiEnabled = options.mpi != 0;
      auto subdomains = constructOwnedSubdomains<BddcSubdomain>(localJ,interfaceAverages,
                                                                 options.assemblyThreads,
                                                                 mpiEnabled,
                                                                 options.masterWorker != 0);
      if (profile)
        profile->subdomainConstruction += ProfileTimes::seconds(subdomainStart);

      auto const transferStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      for (size_t id = 0; id < subdomains.size(); ++id)
        if (subdomains.has(id))
        {
          configureTransfer(subdomains[id].transfer(),options);
          // Corner DOFs are coarse-space representatives.  When requested,
          // keep them out of the compressed payload while retaining the
          // normal BDDC transfer and coarse-space construction.
          if (options.bddcExcludeCornersFromCompression != 0)
            if constexpr (requires(Transfer& t, std::vector<int> const& d) {
                            t.setExcludedLocalDofs(d);
                          })
              subdomains[id].transfer().setExcludedLocalDofs(interfaceAverages.cornerDofs(id));
        }
      if (profile)
        profile->transferSetup += ProfileTimes::seconds(transferStart);

      auto const solverStart = profile ? ProfileTimes::Clock::now() : ProfileTimes::Clock::time_point{};
      Kaskade::BDDC::BDDCSolver<BddcSubdomain,decltype(subdomains)> solver(subdomains,
                                                                            interfaceAverages.coarseConstraints(),
                                                                            solverActiveIds,
                                                                            useCgSolver,
                                                                            verbose,
                                                                            mpiEnabled,
                                                                            options.masterWorker != 0);
      if (options.bddcHistogram)
        solver.setHistogramOutput(options.outputDir + "/histograms");
      if (profile)
      {
        profile->solverConstruction += ProfileTimes::seconds(solverStart);
        ++profile->bddcSolverConstructions;
      }

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
        if (solver.stopped())
        {
          if (verbose || (options.quiet && mpiRankZero(options.mpi != 0)))
            std::cout << "    BDDC stopped before applying the current correction: "
                      << (solver.stagnated() ? "stagnation" : "breakdown") << "\n";
          break;
        }
        if (residual < options.bddcTolerance)
          break;
      }
      if (compressionReport)
      {
        if (options.bddcCompression)
          compressionReport->add(solver.compressionTraffic());
        else
          compressionReport->addUncompressed(solver.traffic(),options.mpi != 0);
      }
      if (profile)
        profile->bddcSolve += ProfileTimes::seconds(solveStart);

      for (int subdomain : activeIds)
        if (subdomains.has(subdomain))
          corrections[i][subdomain] = subdomains[subdomain].getSolution();

      if (verbose)
        std::cout << "    BDDC SDC interval " << i << "/" << intervals
                  << ": iterations=" << std::min(iteration+1,options.bddcIterations)
                  << ", residual=" << residual << "\n";

      for (int subdomain : activeIds)
        norm += (points[i]-points[i-1]) * (corrections[i][subdomain] * rhs[subdomain]);
    }

    double globalNorm = norm;
#ifdef KASKADE_HAVE_MPI
    if (options.mpi != 0)
      MPI_Allreduce(&norm,&globalNorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    return std::sqrt(std::max<typename Matrix::field_type>(0.0,globalNorm/(points[intervals]-points[0])));
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
#ifdef KASKADE_HAVE_MPI
          if (options.mpi != 0)
          {
            Vector reducedCorrection(nDofs);
            MPI_Allreduce(globalCorrection.data(),reducedCorrection.data(),
                          static_cast<int>(nDofs),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            globalCorrection = std::move(reducedCorrection);
          }
#endif
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
        nextActiveIds = synchronizeActiveSubdomains(std::move(nextActiveIds),
                                                    data.localDofs.size(),
                                                    options.mpi != 0);
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

    if (options.profile && !options.quiet)
      profile.print("SDC + BDDC");
    if (options.bddcCompressionReport)
      compressionReport.printGlobal("SDC + BDDC",options.mpi != 0,!options.quiet);
    return state;
  }

  template <class Transfer, class Functional, class Spaces, class State, class Element, class Matrix, class Vector, class Options>
  State runBddc(Functional& F,
                Spaces const& spaces,
                State state,
                Element&,
                BddcData<Matrix,Vector> const& data,
                Matrix const& globalMatrix,
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
                                                              options.bddcInterfaceTypes,
                                                              data.interfaceComponents);
    auto const interfaceCount = interfaceAverages.interfaceCount();
    if (mpiRankZero(options.mpi != 0))
      std::cout << "BDDC coarse interfaces: corners=" << interfaceCount[0]
                << ", edges=" << interfaceCount[1]
                << ", faces=" << interfaceCount[2] << "\n";

    Assembler assembler(spaces);
    auto zeroState(state);
    auto stepState(state);
    CompressionReport compressionReport;
    ProfileTimes profile;
    zeroState *= 0.0;
    stepState *= 0.0;

    // Fixed local operators: retain the two local factorizations across timesteps.
    auto const subdomainStart = options.profile
                              ? ProfileTimes::Clock::now()
                              : ProfileTimes::Clock::time_point{};
    auto subdomains = constructOwnedSubdomains<BddcSubdomain>(data.localMatrices,
                                                               interfaceAverages,
                                                               options.assemblyThreads,
                                                               options.mpi != 0,
                                                               options.masterWorker != 0);
    if (options.profile)
    {
      profile.subdomainConstruction += ProfileTimes::seconds(subdomainStart);
      profile.bddcSubdomainConstructions += subdomains.size();
    }

    auto const transferStart = options.profile
                             ? ProfileTimes::Clock::now()
                             : ProfileTimes::Clock::time_point{};
    if constexpr (requires(Transfer& transfer) {
                    transfer.setQuantizationBits(0);
                    transfer.setRestrictEncoding(true);
                    transfer.setProlongateEncoding(true);
                    transfer.setRestrictTransform(true);
                    transfer.setProlongateTransform(true);
                    transfer.setRestrictBitlengthEncoding(true);
                    transfer.setProlongateBitlengthEncoding(true);
                    transfer.setRestrictZigzagEncoding(true);
                    transfer.setProlongateZigzagEncoding(true);
                  })
    {
      for (size_t id = 0; id < subdomains.size(); ++id)
      {
        if (!subdomains.has(id))
          continue;
        auto& subdomain = subdomains[id];
        configureTransferPrecision(subdomain.transfer(),options);
        subdomain.transfer().setCompressionMinInterfaceSize(
            static_cast<std::size_t>(
                std::max(0, options.bddcCompressionMinInterfaceSize)));
        subdomain.transfer().setCompressionThresholdInterfaceTypes(
            options.bddcCompressionThresholdInterfaceTypes);
        subdomain.transfer().enableTransform((options.bddcGraphLifting || options.bddcGraphLiftingAfterQuantization)
                                               ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                                               : Kaskade::BDDC::TransformType::NONE);
        subdomain.transfer().setRestrictEncoding(options.bddcHuffman != 0);
        subdomain.transfer().setProlongateEncoding(options.bddcHuffman != 0);
        subdomain.transfer().setRestrictTransform(options.bddcGraphLifting != 0);
        subdomain.transfer().setProlongateTransform(options.bddcGraphLifting != 0);
        subdomain.transfer().setRestrictTransformAfterQuantization(options.bddcGraphLiftingAfterQuantization != 0);
        subdomain.transfer().setProlongateTransformAfterQuantization(options.bddcGraphLiftingAfterQuantization != 0);
        subdomain.transfer().setRestrictBitlengthEncoding(options.bddcBitlength != 0);
        subdomain.transfer().setProlongateBitlengthEncoding(options.bddcBitlength != 0);
        subdomain.transfer().setRestrictZigzagEncoding(options.bddcZigzag != 0);
        subdomain.transfer().setProlongateZigzagEncoding(options.bddcZigzag != 0);
        // Apply the same corner policy in the non-SDC BDDC path.  The
        // transfer object performs the filtering consistently for both
        // restriction and prolongation.
        if (options.bddcExcludeCornersFromCompression != 0)
          if constexpr (requires(Transfer& t, std::vector<int> const& d) {
                          t.setExcludedLocalDofs(d);
                        })
            subdomain.transfer().setExcludedLocalDofs(interfaceAverages.cornerDofs(id));
      }
    }
    if (options.profile)
      profile.transferSetup += ProfileTimes::seconds(transferStart);

    if (!options.quiet && mpiRankZero(options.mpi != 0))
      std::cout << "BDDC local and coarse setup reused for " << steps << " time steps\n";

    // The fixed coarse operator is assembled and factored once as well.
    std::vector<int> activeIds(subdomains.size());
    std::iota(activeIds.begin(),activeIds.end(),0);
    bool useCgSolver = options.bddcUseCg != 0;
    // Explicit iteration verbosity also applies in quiet mode.
    // Select one table; the report uses the same accepted updates and stopping checks.
    bool verbose = options.bddcVerbose != 0 && !options.bddcConvergenceReport;

    auto const solverStart = options.profile
                           ? ProfileTimes::Clock::now()
                           : ProfileTimes::Clock::time_point{};
    Kaskade::BDDC::BDDCSolver<BddcSubdomain,decltype(subdomains)> solver(subdomains,
                                                    interfaceAverages.coarseConstraints(),
                                                    activeIds,
                                                    useCgSolver,
                                                    verbose,
                                                    options.mpi != 0,
                                                    options.masterWorker != 0);
    if constexpr (requires { options.bddcRestrictionPasses; })
      solver.setRestrictionPasses(options.bddcRestrictionPasses);
    if (options.profile)
    {
      profile.solverConstruction += ProfileTimes::seconds(solverStart);
      ++profile.bddcSolverConstructions;
    }
    if (options.bddcHistogram)
      solver.setHistogramOutput(options.outputDir + "/histograms");

    // These work vectors have fixed sizes. Each RHS and solution is fully
    // overwritten below; the accumulated local correction is cleared each step.
    Vector rhs(nDofs);
    Vector stepVector(nDofs);
    std::vector<Vector> subdomainRhs;
    subdomainRhs.reserve(data.localDofs.size());
    for (auto const& dofs : data.localDofs)
      subdomainRhs.emplace_back(dofs.size());
    std::vector<typename BddcSubdomain::XVector> subdomainSolutions(subdomains.size());
    for (size_t id = 0; id < subdomains.size(); ++id)
      if (subdomains.has(id))
        subdomainSolutions[id].resize(data.localDofs[id].size());
    std::vector<double> localStep(nDofs, 0.0);
    std::vector<double> globalStep(nDofs, 0.0);
    bool const rankZero = mpiRankZero(options.mpi != 0);
    bool const relativeStopping = options.bddcRelativeTolerance > 0.0;
    auto checkConvergenceOutput = [&](bool success)
    {
      if (!options.bddcConvergenceReport)
        return;
      int ready = success ? 1 : 0;
#ifdef KASKADE_HAVE_MPI
      if (options.mpi != 0)
        MPI_Bcast(&ready,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      if (!ready)
        throw std::runtime_error("Cannot write BDDC convergence CSV in " + options.outputDir);
    };
    bool const computeDiagnostics = options.diagnostics && rankZero;
    Vector diagnosticResidual((computeDiagnostics || (relativeStopping && rankZero)) ? nDofs : 0);

    // Assemble the actual accepted correction, not the CG search direction.
    // This is needed after every accepted update in true-residual mode.
    auto assembleCurrentStep = [&]()
    {
      auto const start = options.profile ? ProfileTimes::Clock::now()
                                         : ProfileTimes::Clock::time_point{};
      std::fill(localStep.begin(),localStep.end(),0.0);
      for (size_t subdomain = 0; subdomain < subdomains.size(); ++subdomain)
      {
        if (!subdomains.has(subdomain))
          continue;
        auto& localSolution = subdomainSolutions[subdomain];
        subdomains[subdomain].getSolution(localSolution);
        for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        {
          size_t const global = data.localDofs[subdomain][local];
          localStep[global] += data.weights[subdomain][local][0] * localSolution[local];
        }
      }
#ifdef KASKADE_HAVE_MPI
      if (options.mpi != 0)
        MPI_Allreduce(localStep.data(),globalStep.data(),static_cast<int>(nDofs),
                      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      else
#endif
        globalStep = localStep;
      for (size_t i = 0; i < nDofs; ++i)
        stepVector[i] = globalStep[i];
      if (options.profile)
        profile.correctionAssembly += ProfileTimes::seconds(start);
    };
    if (relativeStopping && rankZero)
      std::cout << "BDDC stopping rule: ||b-A*du||_2 / ||b||_2 <= "
                << options.bddcRelativeTolerance
                << " (zero RHS: absolute residual; legacy internal tolerance is inactive)\n";
    if (verbose && rankZero)
      std::cout << "BDDC table legend: step length=alpha; ||du||_A=energy norm of the unscaled search direction.\n"
                << "  p[0]=internal residual dot BDDC-preconditioned residual; not ||b||_2.\n"
                << "  Table values are evaluated before each update; du in the stopping rule is the accumulated correction.\n"
                << "  rel uses the assembled global residual and fixed timestep RHS norm (absolute residual for zero RHS).\n";
    if (options.bddcConditionEstimate && rankZero && !options.quiet)
      std::cout << "BDDC condition diagnostic: RHS-dependent CG/Lanczos Ritz estimate "
                << "after interior elimination; fixed SPD preconditioner assumed, "
                << "not the exact condition number of A\n";

    for (int step = 0; step < steps; ++step)
    {
      double const time = step*options.dt;
      F.time(time);
      F.Mass_stiff(10);

      Kaskade::SemiImplicitEulerStep<Functional> equation(&F,options.dt);
      equation.setTau(options.dt);
      auto const residualAssemblyStart = options.profile
                                       ? ProfileTimes::Clock::now()
                                       : ProfileTimes::Clock::time_point{};
      assembler.assemble(SemiLinearization(equation,state,state,zeroState),
                         Assembler::RHS,
                         options.assemblyThreads);
      if (options.profile)
        profile.residualAssembly += ProfileTimes::seconds(residualAssemblyStart);

      assembler.rhs().write(rhs.begin());

      auto const rhsSetupStart = options.profile
                               ? ProfileTimes::Clock::now()
                               : ProfileTimes::Clock::time_point{};
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
      {
        for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
        {
          size_t const global = data.localDofs[subdomain][local];
          subdomainRhs[subdomain][local] = data.weights[subdomain][local][0] * rhs[global];
        }
      }
      if (options.profile)
        profile.rhsSetup += ProfileTimes::seconds(rhsSetupStart);
      if (options.diagnostics && step == 0)
        printLocalRhsReconstruction(data,subdomainRhs,rhs,options.mpi != 0);

      auto const solverRhsStart = options.profile
                                ? ProfileTimes::Clock::now()
                                : ProfileTimes::Clock::time_point{};
      solver.restartRhs(subdomainRhs);
      if (options.profile)
      {
        profile.rhsSetup += ProfileTimes::seconds(solverRhsStart);
        ++profile.bddcSystems;
      }

      double const rhsNorm = rankZero && (computeDiagnostics || relativeStopping)
                           ? std::sqrt(rhs.dot(rhs)) : 0.0;
      auto measureCurrentResidual = [&]()
      {
        auto const start = options.profile ? ProfileTimes::Clock::now()
                                           : ProfileTimes::Clock::time_point{};
        TrueResidual measured;
        if (rankZero)
        {
          diagnosticResidual = rhs;
          diagnosticResidual *= -1.0;
          globalMatrix.umv(stepVector,diagnosticResidual);
          measured = {std::sqrt(diagnosticResidual.dot(diagnosticResidual)),rhsNorm};
        }
#ifdef KASKADE_HAVE_MPI
        // The global matrix, RHS and correction are replicated. Evaluate their
        // norm once on root; summing replicated norms would overcount them.
        // Broadcast the measurement so every rank makes the same stop decision.
        if (relativeStopping && options.mpi != 0)
        {
          double norms[2] = {measured.absolute,measured.rhsNorm};
          MPI_Bcast(norms,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
          measured = {norms[0],norms[1]};
        }
#endif
        if (options.profile)
        {
          profile.trueResidualChecks += ProfileTimes::seconds(start);
          ++profile.trueResidualCheckCalls;
        }
        return measured;
      };

      ConvergenceReport convergenceReport;
      std::string const convergencePath = options.bddcConvergenceReport
        ? options.outputDir + "/bddc-convergence-step" + std::to_string(step+1) + ".csv"
        : std::string{};
      if (options.bddcConvergenceReport)
        checkConvergenceOutput(!rankZero || convergenceReport.open(
          convergencePath,step+1,options.bddcRelativeTolerance,std::cout));
      int reportedUpdates = 0;
      std::optional<double> reportedAlpha;

      double residual = std::numeric_limits<double>::infinity();
      int solveCalls = 0;
      CgConditionEstimator conditionEstimator;
      auto recordCgStep = [&]()
      {
        if (options.bddcConditionEstimate && rankZero)
        {
          auto const& cgStep = solver.lastCgStep();
          if (cgStep.applied)
            conditionEstimator.append(cgStep.alpha,cgStep.beta);
        }
      };
      TrueResidualResult trueResult;
      if (relativeStopping)
      {
        trueResult = solveToTrueResidual(options.bddcIterations,options.bddcRelativeTolerance,
          [&]()
          {
            auto const start = options.profile ? ProfileTimes::Clock::now()
                                               : ProfileTimes::Clock::time_point{};
            double const internalNorm = solver.solve();
            recordCgStep();
            if (options.bddcConvergenceReport)
            {
              auto const& cgStep = solver.lastCgStep();
              if (cgStep.applied)
              {
                ++reportedUpdates;
                reportedAlpha = cgStep.alpha;
              }
            }
            if (options.profile)
            {
              profile.bddcSolve += ProfileTimes::seconds(start);
              ++profile.bddcSolveCalls;
            }
            return internalNorm;
          },
          [&]() { return solver.stopped(); },
          [&]()
          {
            assembleCurrentStep();
            TrueResidual const measured = measureCurrentResidual();
            if (options.bddcConvergenceReport && rankZero)
              convergenceReport.record(reportedUpdates,reportedAlpha,measured);
            return measured;
          });
        solveCalls = trueResult.calls;
        residual = trueResult.internalNorm;
      }
      else
      {
        auto const solveStart = options.profile
                              ? ProfileTimes::Clock::now()
                              : ProfileTimes::Clock::time_point{};
        for (int iteration = 0; iteration < options.bddcIterations; ++iteration)
        {
          if (options.profile)
            ++profile.bddcSolveCalls;
          ++solveCalls;
          residual = solver.solve();
          recordCgStep();
          if (solver.stopped())
          {
            if (verbose || (options.quiet && rankZero))
              std::cout << "  BDDC stopped before applying the current correction: "
                        << (solver.stagnated() ? "stagnation" : "breakdown") << "\n";
            break;
          }
          if (residual < options.bddcTolerance)
            break;
        }
        if (options.profile)
          profile.bddcSolve += ProfileTimes::seconds(solveStart);
        assembleCurrentStep();
      }
      if (options.bddcConvergenceReport)
      {
        checkConvergenceOutput(!rankZero || convergenceReport.close());
        if (rankZero)
          std::cout << "BDDC convergence data: " << convergencePath << "\n";
      }
      if (relativeStopping || computeDiagnostics)
      {
        TrueResidual const measured = relativeStopping ? trueResult.measured
                                                        : measureCurrentResidual();
        if (rankZero)
        {
          if (options.quiet)
            std::cout << "BDDC step " << step+1 << "/" << steps
                      << ": iterations=" << solveCalls << " rel=" << measured.relative();
          else
            std::cout << "BDDC global residual step " << step+1
                    << ": abs=" << measured.absolute << " rel=" << measured.relative();
          if (relativeStopping && options.quiet)
            std::cout << " status=" << (trueResult.converged() ? "converged" : "not-converged");
          if (relativeStopping && !options.quiet)
            std::cout << " target=" << options.bddcRelativeTolerance
                      << " status=" << (trueResult.converged() ? "converged" : "not-converged");
          if (!options.quiet)
            std::cout << " iterations=" << solveCalls;
          if (!options.quiet || !options.bddcConditionEstimate)
            std::cout << "\n";
        }
      }
      if (options.quiet && rankZero && !relativeStopping && !computeDiagnostics)
      {
        std::cout << "BDDC step " << step+1 << "/" << steps
                  << ": iterations=" << solveCalls << " internal_norm=" << residual;
        if (!options.bddcConditionEstimate)
          std::cout << "\n";
      }
      if (options.bddcConditionEstimate && rankZero)
      {
        auto const estimate = conditionEstimator.estimate();
        if (options.quiet)
        {
          if (std::string(estimate.status) == "estimated")
            std::cout << " condition_estimate=" << estimate.condition << "\n";
          else
            std::cout << " condition_estimate=unavailable(" << estimate.status << ")\n";
        }
        else
        {
          std::cout << "BDDC condition estimate step " << step+1;
          if (std::string(estimate.status) == "estimated")
            std::cout << ": lambda_min=" << estimate.lambdaMin
                      << " lambda_max=" << estimate.lambdaMax
                      << " condition_estimate=" << estimate.condition;
          else
            std::cout << ": condition_estimate=unavailable";
          std::cout << " krylov_dim=" << estimate.dimension
                    << " status=" << estimate.status << "\n";
        }
      }
      if (options.quiet && rankZero)
        std::cout << std::flush;
      if (relativeStopping && !trueResult.converged())
      {
        char const* reason = trueResult.status == TrueResidualStatus::iterationLimit
                           ? "iteration limit"
                           : trueResult.status == TrueResidualStatus::nonfiniteResidual
                           ? "nonfinite true residual or RHS norm"
                           : solver.stagnated() ? "solver stagnation" : "solver breakdown";
        if (rankZero)
          std::cout << "BDDC step " << step+1 << " failed after " << solveCalls
                    << " solver calls: " << reason << ". State was not advanced.\n" << std::flush;
        throw std::runtime_error(std::string("BDDC true relative residual target not reached: ")
                                 + reason);
      }

      for (size_t i = 0; i < nDofs; ++i)
        boost::fusion::at_c<0>(stepState.data).coefficients()[i] = stepVector[i];

      Kaskade::component<0>(state) += Kaskade::component<0>(stepState);
      if (verbose && rankZero && !options.quiet)
      {
        std::cout << "BDDC step " << step+1 << "/" << steps
                  << ": iterations=" << solveCalls;
        if (relativeStopping)
          std::cout << ", true relative residual=" << trueResult.measured.relative() << "\n";
        else
          std::cout << ", residual=" << residual << "\n";
      }
    }

    // The retained solver reports cumulative traffic. Add it once at the end.
    if (options.bddcCompressionReport)
    {
      if (options.bddcCompression)
        compressionReport.add(solver.compressionTraffic());
      else
        compressionReport.addUncompressed(solver.traffic(),options.mpi != 0);
    }

    if (options.profile && !options.quiet)
      profile.print("BDDC");
    if (options.bddcCompressionReport)
      compressionReport.printGlobal("BDDC",options.mpi != 0,!options.quiet);
    return state;
  }
}

#endif
