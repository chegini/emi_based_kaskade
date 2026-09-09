#ifndef EMI_BDDC_HH
#define EMI_BDDC_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
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

namespace EmiBddc
{
  template <class Matrix, class Vector>
  struct BddcData
  {
    std::vector<int> tags;
    std::vector<std::vector<size_t>> localDofs;
    std::vector<std::vector<Kaskade::BDDC::LocalDof>> sharedDofs;
    std::vector<int> subdomainSizes;
    std::vector<Matrix> localMatrices;
    std::vector<Vector> weights;
  };

  template <class Grid, class Space, class Material, class Matrix, class Vector>
  BddcData<Matrix,Vector> buildBddcData(Grid const& grid,
                                        Space const& space,
                                        Material const& material,
                                        Matrix const& globalMatrix,
                                        size_t nDofs)
  {
    using Kaskade::BDDC::LocalDof;

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

    if (data.sharedDofs.empty() && data.tags.size() > 1)
      throw std::runtime_error("BDDC setup found multiple subdomains but no shared dofs.");

    data.localMatrices.reserve(data.localDofs.size());
    data.weights.reserve(data.localDofs.size());
    for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
    {
      data.localMatrices.emplace_back(data.localDofs[subdomain],globalMatrix);

      Vector weight(data.localDofs[subdomain].size());
      weight = 1.0;
      for (size_t local = 0; local < data.localDofs[subdomain].size(); ++local)
      {
        size_t const global = data.localDofs[subdomain][local];
        size_t const ownerCount = std::max<size_t>(1,dofOccurrences[global].size());
        weight[local][0] = 1.0/static_cast<double>(ownerCount);
      }
      data.weights.push_back(weight);
    }

    std::cout << "BDDC subdomains: " << data.tags.size()
              << ", shared dofs: " << data.sharedDofs.size() << "\n";
    return data;
  }

  template <class Matrix, class Vector>
  std::vector<Matrix> extractLocalMatrices(BddcData<Matrix,Vector> const& data,
                                           Matrix const& globalMatrix)
  {
    std::vector<Matrix> localMatrices;
    localMatrices.reserve(data.localDofs.size());
    for (auto const& localDofs : data.localDofs)
      localMatrices.emplace_back(localDofs,globalMatrix);
    return localMatrices;
  }

  template <class Matrix, class Vector>
  std::vector<Vector> scatterToSubdomains(BddcData<Matrix,Vector> const& data,
                                          Vector const& global)
  {
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

  template <class Functional, class State, class StateU, class TimeGrid, class Assembler, class Matrix, class Vector, class Options>
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
                               std::vector<std::vector<Vector>>& residuals)
  {
    using SemiLinearization = Kaskade::SemiLinearizationAtInner<Kaskade::SemiImplicitEulerStep<Functional>>;

    State stateTmp(stateAtStart);
    State dstateTmp(stateAtStart);
    dstateTmp *= 0.0;
    equation.setTau(dt);
    F.Mass_stiff(0);
    Vector fullResidual(nDofs);

    auto const& points = grid.points();
    for (int i = 0; i < points.N(); ++i)
    {
      if (sweep == 0 && i > 0)
      {
        residuals[i] = residuals[i-1];
        continue;
      }

      F.time(t + points[i] - points[0]);
      boost::fusion::at_c<0>(stateTmp.data) = collocationStates[i];
      assembler.assemble(SemiLinearization(equation,stateTmp,stateTmp,dstateTmp),
                         Assembler::RHS,
                         options.assemblyThreads);

      auto rhs = assembler.rhs();
      rhs.write(fullResidual.begin());
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

  template <class Matrix, class Vector, class Options>
  typename Matrix::field_type sdcIterationStepBddc(Kaskade::SDCTimeGrid const& grid,
                                                   Kaskade::SDCTimeGrid::RealMatrix const& Shat,
                                                   BddcData<Matrix,Vector> const& data,
                                                   std::vector<Matrix> const& localMassMatrices,
                                                   std::vector<Matrix> const& localStiffnessMatrices,
                                                   std::vector<std::vector<Vector>> const& residuals,
                                                   std::vector<std::vector<Vector>> const& massDifferences,
                                                   std::vector<std::vector<Vector>>& corrections,
                                                   Options const& options)
  {
    using TransmissionScalar = double;
    using BddcSubdomain = Kaskade::BDDC::Subdomain<1,double,double,Kaskade::BDDC::SpaceTransfer<1,double,TransmissionScalar>>;

    auto const& points = grid.points();
    int const intervals = points.size()-1;
    auto const& S = grid.integrationMatrix();
    typename Matrix::field_type norm = 0.0;

    for (auto& localCorrection : corrections[0])
      localCorrection = 0.0;

    Kaskade::BDDC::InterfaceAverages<1,int> interfaceAverages(data.sharedDofs,
                                                              data.subdomainSizes,
                                                              options.bddcInterfaceTypes);

    std::vector<int> activeIds(data.localDofs.size());
    std::iota(activeIds.begin(),activeIds.end(),0);
    bool useCgSolver = options.bddcUseCg != 0;
    bool verbose = options.bddcVerbose != 0;

    for (int i = 1; i <= intervals; ++i)
    {
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

      std::vector<BddcSubdomain> subdomains;
      subdomains.reserve(data.localDofs.size());
      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        subdomains.emplace_back(static_cast<int>(subdomain),localJ[subdomain],interfaceAverages);

      Kaskade::BDDC::BDDCSolver<BddcSubdomain> solver(subdomains,
                                                      interfaceAverages.coarseConstraints(),
                                                      activeIds,
                                                      useCgSolver,
                                                      verbose);
      solver.setRhs(rhs);

      double residual = std::numeric_limits<double>::infinity();
      int iteration = 0;
      for (; iteration < options.bddcIterations; ++iteration)
      {
        residual = solver.solve();
        if (residual < options.bddcTolerance)
          break;
      }

      for (size_t subdomain = 0; subdomain < subdomains.size(); ++subdomain)
        corrections[i][subdomain] = subdomains[subdomain].getSolution();

      if (verbose)
        std::cout << "    BDDC SDC interval " << i << "/" << intervals
                  << ": iterations=" << std::min(iteration+1,options.bddcIterations)
                  << ", residual=" << residual << "\n";

      for (size_t subdomain = 0; subdomain < data.localDofs.size(); ++subdomain)
        norm += (points[i]-points[i-1]) * (corrections[i][subdomain] * rhs[subdomain]);
    }

    return std::sqrt(std::max<typename Matrix::field_type>(0.0,norm/(points[intervals]-points[0])));
  }

  template <class Functional, class Spaces, class State, class Element, class Matrix, class Vector, class Options>
  State runBddcSdc(Functional& F,
                   Spaces const& spaces,
                   State state,
                   Element&,
                   BddcData<Matrix,Vector> const& data,
                   Matrix const& mass,
                   Matrix const& stiffness,
                   size_t nDofs,
                   int steps,
                   Options const& options)
  {
    using SemiLinearization = Kaskade::SemiLinearizationAtInner<Kaskade::SemiImplicitEulerStep<Functional>>;
    using Assembler = Kaskade::VariationalFunctionalAssembler<SemiLinearization>;
    using StateU = typename boost::fusion::result_of::value_at_c<typename State::Sequence,0>::type;

    if (options.algebraicAdaptivity)
      throw std::runtime_error("BDDC+SDC algebraic adaptivity is not enabled yet. Run without --algebraicAdaptivity.");

    auto const localMassMatrices = extractLocalMatrices(data,mass);
    auto const localStiffnessMatrices = extractLocalMatrices(data,stiffness);

    Assembler assembler(spaces);
    Kaskade::SemiImplicitEulerStep<Functional> equation(&F,options.dt);

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
      int sweep = 0;

      std::cout << "BDDC+SDC step " << step+1 << "/" << steps
                << ", time [" << t << ", " << stepEnd << "]\n";

      for (; sweep < options.maximumSdcSweeps; ++sweep)
      {
        std::vector<std::vector<Vector>> residuals(grid.points().N());
        computeBddcSdcResiduals(F,equation,state,collocationStates,grid,data,nDofs,
                                sweep,t,stepDt,assembler,options,residuals);

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

        std::vector<std::vector<Vector>> corrections(grid.points().N());
        for (auto& pointData : corrections)
          for (auto const& localDofs : data.localDofs)
          {
            pointData.emplace_back(localDofs.size());
            pointData.back() = 0.0;
          }

        sweepNorms.push_back(sdcIterationStepBddc(grid,Shat,data,localMassMatrices,localStiffnessMatrices,
                                                  residuals,massDifferences,corrections,options));

        for (int i = 1; i < grid.points().N(); ++i)
        {
          Vector globalCorrection = combineSubdomainVector(data,corrections[i],nDofs);
          for (size_t j = 0; j < nDofs; ++j)
            collocationStates[i].coefficients()[j] += globalCorrection[j];
        }

        if (sweepNorms.size() > 1)
        {
          double const c = sweepNorms.back()/sweepNorms[sweepNorms.size()-2];
          sdcContraction = std::sqrt(c*sdcContraction);
        }

        std::cout << "  sweep " << sweep+1
                  << ": ||du||=" << sweepNorms.back()
                  << ", contraction=" << sdcContraction << "\n";

        bool const reachedMinimumSweeps = sweep+1 >= options.minimumSdcSweeps;
        bool const smallCorrection = sweepNorms.back() < options.sdcTolerance;
        bool const reliableContraction = sdcContraction < 1.0;
        bool const estimatedSmall = reliableContraction
                                 && sweepNorms.back()*sdcContraction/(1.0-sdcContraction) <= options.sdcAbsoluteTolerance;
        if (reachedMinimumSweeps && (smallCorrection || estimatedSmall))
          break;

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

    return state;
  }

  template <class Functional, class Spaces, class State, class Element, class Matrix, class Vector, class Options>
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
    using BddcSubdomain = Kaskade::BDDC::Subdomain<1,double,double,Kaskade::BDDC::SpaceTransfer<1,double,TransmissionScalar>>;

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

      std::vector<BddcSubdomain> subdomains;
      subdomains.reserve(data.localMatrices.size());
      for (size_t subdomain = 0; subdomain < data.localMatrices.size(); ++subdomain)
        subdomains.emplace_back(static_cast<int>(subdomain),data.localMatrices[subdomain],interfaceAverages);

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

    return state;
  }
}

#endif
