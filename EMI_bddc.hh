#ifndef EMI_BDDC_HH
#define EMI_BDDC_HH

#include <algorithm>
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
