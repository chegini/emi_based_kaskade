#define FUSION_MAX_VECTOR_SIZE 20

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/timer/timer.hpp>
#include <cstdint>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/gridmanager.hh"
#include "fem/istlinterface.hh"
#include "fem/lagrangespace.hh"
#include "fem/partitionedspace.hh"
#include "fem/spaces.hh"
#include "fem/variables.hh"
#include "io/vtk.hh"
#include "io/vtkreader.hh"
#include "linalg/apcg.hh"
#include "linalg/direct.hh"
#include "linalg/dynamicMatrix.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "timestepping/semieuler.hh"
#include "timestepping/sdc.hh"
#include "utilities/kaskopt.hh"
#include "utilities/timing.hh"

#include "fem/diffops/membraneModels.hh"

#include "EMI_bddc.hh"
#include "EMI_model.hh"

using namespace Kaskade;

#ifndef SPACEDIM
#define SPACEDIM 3
#endif

struct EmiOptions
{
  std::string input = "./input/2Cells3d_2extra_mesh.vtu";
  std::string extraSet = "./input/2Cells3d_2extra_list_extracellular.txt";
  std::string intraSet = "./input/2Cells3d_2extra_list_intracellular.txt";
  std::string excitedSet = "./input/2Cells3d_2extra_early_excited.txt";
  std::string outputDir = "output";

  int refinements = 0;
  int order = 1;
  int maximumNumberOfTimeSteps = 0;
  int assemblyThreads = 1;
  int maxCGIter = 10000;
  int writeVTK = 1;
  int diagnostics = 0;
  int sdc = 0;
  int sdcCollocationPoints = 3;
  int sdcStartCollocationPoints = 3;
  int minimumSdcSweeps = 3;
  int maximumSdcSweeps = 5;
  int sdcSweepType = 1;
  int algebraicAdaptivity = 0;
  int bddc = 0;
  int bddcCompression = 0;
  int bddcCompressionBits = 16;
  int bddcIterations = 3000;
  int bddcInterfaceTypes = 7;
  int bddcVerbose = 0;
  int bddcUseCg = 1;

  double finalTime = 0.01;
  double dt = 0.01;
  double cgTol = 1e-8;
  double sdcTolerance = 1e-6;
  double sdcAbsoluteTolerance = 1e-12;
  double sdcInitialContraction = 0.2;
  double algebraicAdaptivityTolerance = 0.0;
  double bddcTolerance = 1e-8;
  double penalty = 1e6;
  double sigmaI = 3.0;
  double sigmaE = 20.0;
  double membraneCapacitance = 1.0;
  double gapConductance = 0.1;
  double extraConductance = 1e-7;

  bool direct = false;
  bool cgShift = true;
};

void readTagList(std::string const& filename, std::vector<int>& tags)
{
  std::ifstream file(filename);
  if (!file)
    throw std::runtime_error("Cannot open tag file: " + filename);

  int count = 0;
  file >> count;
  tags.resize(count);
  for (int& tag : tags)
    file >> tag;
}

template <typename Material>
struct InitialValue
{
  using Scalar = double;
  static constexpr int components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  InitialValue(Material const& material_, std::vector<int> excitedTags_)
  : material(material_), excitedTags(std::move(excitedTags_))
  {}

  template <class Cell>
  int order(Cell const&) const { return std::numeric_limits<int>::max(); }

  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> const&) const
  {
    Dune::FieldVector<double,SPACEDIM> zero(0);
    int const tag = material.value(cell,zero);
    if (std::find(excitedTags.begin(), excitedTags.end(), tag) != excitedTags.end())
      return 0.5;
    return 0.0;
  }

private:
  Material const& material;
  std::vector<int> excitedTags;
};

template <class VariableSet, class Element>
void writeState(VariableSet const& u, Element& uAll, int order, std::string const& filename)
{
  uAll = component<0>(u);
  writeVTK(uAll,filename,
           IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),
           "u");
}

template <class Matrix, class Vector>
void applyConstantShift(Matrix const& A, Vector& du, Vector const& rhs)
{
  // EMI systems with pure Neumann-type coupling can contain a constant-mode
  // ambiguity. After PCG, correct only that constant shift without changing
  // the resolved potential differences.
  Vector ones(du.size());
  Vector Aones(du.size());
  ones = 1.0;
  Aones = 0.0;

  A.umv(ones,Aones);

  double const numerator = ones.dot(rhs) - du.dot(Aones);
  double const denominator = ones.dot(Aones);
  if (std::abs(denominator) > 0.0)
  {
    ones *= numerator/denominator;
    du += ones;
  }
}

template <class Matrix, class Vector>
void solveLinearSystem(Matrix const& A, Vector& du, Vector& rhs, EmiOptions const& options)
{
  du = 0.0;

  if (options.direct)
  {
    DirectSolver<Vector,Vector> solver(A);
    solver.apply(du,rhs);
    return;
  }

  DefaultDualPairing<Vector,Vector> dp;
  Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);
  JacobiPreconditioner<Dune::MatrixAdapter<Matrix,Vector,Vector>> preconditioner(op,1.0);
  PCGEnergyErrorTerminationCriterion<double> term(options.cgTol,options.maxCGIter);
  Pcg<Vector,Vector> pcg(op,preconditioner,dp,term,0);
  Dune::InverseOperatorResult result;
  pcg.apply(du,rhs,result);

  if (options.cgShift)
    applyConstantShift(A,du,rhs);
}

template <class Matrix>
void printMatrixDiagnostics(Matrix const& A, std::string const& name)
{
  size_t nanEntries = 0;
  size_t infEntries = 0;
  size_t numericalNonzeros = 0;
  size_t zeroDiagonal = 0;
  size_t missingDiagonal = 0;
  double minValue = std::numeric_limits<double>::infinity();
  double maxValue = -std::numeric_limits<double>::infinity();
  double maxAbs = 0.0;

  for (auto row = A.begin(); row != A.end(); ++row)
  {
    bool sawDiagonal = false;
    double diagonal = 0.0;
    for (auto col = row.begin(); col != row.end(); ++col)
    {
      double const value = (*col)[0][0];
      if (std::isnan(value))
        ++nanEntries;
      if (std::isinf(value))
        ++infEntries;
      if (std::isfinite(value))
      {
        minValue = std::min(minValue,value);
        maxValue = std::max(maxValue,value);
        maxAbs = std::max(maxAbs,std::abs(value));
      }
      if (value != 0.0)
        ++numericalNonzeros;
      if (row.index() == col.index())
      {
        sawDiagonal = true;
        diagonal = value;
      }
    }

    if (!sawDiagonal)
      ++missingDiagonal;
    else if (diagonal == 0.0)
      ++zeroDiagonal;
  }

  std::cout << name << ": N=" << A.N()
            << " M=" << A.M()
            << " nnz=" << A.nonzeroes()
            << " numericalNonzeros=" << numericalNonzeros
            << " nan=" << nanEntries
            << " inf=" << infEntries
            << " missingDiagonal=" << missingDiagonal
            << " zeroDiagonal=" << zeroDiagonal
            << " min=" << minValue
            << " max=" << maxValue
            << " maxAbs=" << maxAbs << "\n";
}

template <class Vector>
void printVectorDiagnostics(Vector const& v, std::string const& name)
{
  size_t nanEntries = 0;
  size_t infEntries = 0;
  double minValue = std::numeric_limits<double>::infinity();
  double maxValue = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < v.N(); ++i)
  {
    double const value = v[i][0];
    if (std::isnan(value))
      ++nanEntries;
    if (std::isinf(value))
      ++infEntries;
    if (std::isfinite(value))
    {
      minValue = std::min(minValue,value);
      maxValue = std::max(maxValue,value);
    }
  }

  std::cout << name << ": N=" << v.N()
            << " nan=" << nanEntries
            << " inf=" << infEntries
            << " min=" << minValue
            << " max=" << maxValue << "\n";
}

template <class Matrix, class Vectors, class ReactionDerivatives, class Solver>
typename Matrix::field_type sdcIterationStepJacobi(SDCTimeGrid const& grid,
                                                   SDCTimeGrid::RealMatrix const& Shat,
                                                   Solver const& solve,
                                                   Matrix const& M,
                                                   Matrix const& stiffness,
                                                   Vectors const& residuals,
                                                   ReactionDerivatives const& reactionDerivatives,
                                                   Vectors const& massDifferences,
                                                   Vectors& corrections)
{
  auto const& points = grid.points();
  int const intervals = points.size()-1;
  using Vector = typename Vectors::value_type;

  corrections[0] = 0.0;
  Matrix J = M;
  Vector rhs(massDifferences[0].size());
  Vector tmp(massDifferences[0].size());
  auto const& S = grid.integrationMatrix();
  typename Matrix::field_type norm = 0.0;

  for (int i = 1; i <= intervals; ++i)
  {
    for (size_t row = 0; row < J.N(); ++row)
    {
      auto colJ = J[row].begin();
      auto const endJ = J[row].end();
      auto colM = M[row].begin();
      auto colA = stiffness[row].begin();
      auto colR = reactionDerivatives[i][row].begin();
      auto const endR = reactionDerivatives[i][row].end();

      while (colJ != endJ)
      {
        *colJ = *colM - Shat[i-1][i] * *colA;
        if (colR != endR && colJ.index() == colR.index())
        {
          *colJ -= std::min(0.5 * *colM, Shat[i-1][i] * *colR);
          ++colR;
        }
        ++colJ;
        ++colM;
        ++colA;
      }
    }

    rhs = massDifferences[i-1];
    M.umv(corrections[i-1],rhs);
    for (int j = 0; j <= intervals; ++j)
      rhs.axpy(S[i-1][j],residuals[j]);

    tmp = 0.0;
    for (int j = 0; j < i; ++j)
    {
      reactionDerivatives[j].usmv(Shat[i-1][j],corrections[j],rhs);
      tmp.axpy(Shat[i-1][j],corrections[j]);
    }
    stiffness.umv(tmp,rhs);

    corrections[i] = corrections[i-1];
    solve(J,corrections[i],rhs);
    norm += (points[i]-points[i-1]) * (corrections[i] * rhs);
  }

  return std::sqrt(std::max<typename Matrix::field_type>(0.0,norm/(points[intervals]-points[0])));
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

template <class Functional, class State, class StateU, class TimeGrid, class Assembler, class Vector>
void computeFullSdcResiduals(Functional& F,
                             SemiImplicitEulerStep<Functional>& equation,
                             State const& stateAtStart,
                             std::vector<StateU> const& collocationStates,
                             TimeGrid const& grid,
                             std::vector<size_t> const& expandedIndices,
                             size_t nDofs,
                             int sweep,
                             double t,
                             double dt,
                             Assembler& assembler,
                             EmiOptions const& options,
                             std::vector<Vector>& residuals)
{
  using SemiLinearization = SemiLinearizationAtInner<SemiImplicitEulerStep<Functional>>;

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
    for (size_t j = 0; j < expandedIndices.size(); ++j)
      residuals[i][j] = fullResidual[expandedIndices[j]];
    residuals[i] *= (1.0/dt);
  }
}

template <class Matrix, class StateU, class Vector>
void computeFullMassDifferences(Matrix const& M,
                                std::vector<StateU> const& collocationStates,
                                std::vector<size_t> const& expandedIndices,
                                std::vector<Vector>& massDifferences)
{
  for (int i = 0; i < static_cast<int>(collocationStates.size())-1; ++i)
    massDifferences[i] = 0.0;

  for (size_t localRow = 0; localRow < expandedIndices.size(); ++localRow)
  {
    size_t const row = expandedIndices[localRow];
    for (auto col = M[row].begin(); col != M[row].end(); ++col)
    {
      for (int i = 0; i < static_cast<int>(collocationStates.size())-1; ++i)
      {
        auto const diff = collocationStates[i].coefficients()[col.index()]
                        - collocationStates[i+1].coefficients()[col.index()];
        massDifferences[i][localRow] += *col * diff;
      }
    }
  }
}

template <class Functional, class Spaces, class State, class Element>
State runFullSdc(Functional& F,
                 Spaces const& spaces,
                 State state,
                 Element& uAll,
                 size_t nDofs,
                 int steps,
                 EmiOptions const& options)
{
  using SemiLinearization = SemiLinearizationAtInner<SemiImplicitEulerStep<Functional>>;
  using Assembler = VariationalFunctionalAssembler<SemiLinearization>;
  using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
  using StateU = typename boost::fusion::result_of::value_at_c<typename State::Sequence,0>::type;
  using DiagonalMatrix = Dune::BDMatrix<typename Matrix::block_type>;

  Assembler assembler(spaces);
  State zeroState(state);
  zeroState *= 0.0;

  SemiImplicitEulerStep<Functional> equation(&F,options.dt);
  equation.setTau(0.0);
  F.Mass_stiff(1);
  assembler.assemble(SemiLinearization(equation,state,state,zeroState),
                     Assembler::MATRIX|Assembler::RHS,
                     options.assemblyThreads);
  Matrix mass = assembler.template get<Matrix>(false);

  equation.setTau(1.0);
  F.Mass_stiff(0);
  assembler.assemble(SemiLinearization(equation,state,state,zeroState),
                     Assembler::MATRIX|Assembler::RHS,
                     options.assemblyThreads);
  Matrix stiffness = assembler.template get<Matrix>(false);

  printMatrixDiagnostics(mass,"sdc mass");
  printMatrixDiagnostics(stiffness,"sdc stiffness");
  auto const dofNeighborhood = buildMatrixNeighborhood(mass,stiffness);

  for (int step = 0; step < steps; ++step)
  {
    double const t = step * options.dt;
    double const stepEnd = std::min(options.finalTime,t+options.dt);
    double const stepDt = stepEnd - t;
    F.time(t);

    RadauTimeGrid grid(options.sdcStartCollocationPoints,t,stepEnd);
    std::vector<StateU> collocationStates(grid.points().N(),component<0>(state));
    std::vector<size_t> expandedIndices(nDofs);
    std::iota(expandedIndices.begin(),expandedIndices.end(),0);
    std::vector<size_t> compressedIndex(nDofs);
    std::iota(compressedIndex.begin(),compressedIndex.end(),0);

    double sdcContraction = options.sdcInitialContraction;
    std::vector<double> sweepNorms;
    int sweep = 0;

    std::cout << "SDC step " << step+1 << "/" << steps
              << ", time [" << t << ", " << stepEnd << "]\n";

    for (; sweep < options.maximumSdcSweeps; ++sweep)
    {
      size_t const activeDofs = expandedIndices.size();
      std::vector<Vector> residuals(grid.points().N(),Vector(activeDofs));
      computeFullSdcResiduals(F,equation,state,collocationStates,grid,expandedIndices,nDofs,
                              sweep,t,stepDt,assembler,options,residuals);

      SDCTimeGrid::RealMatrix Shat;
      if (options.sdcSweepType == 0)
        eulerIntegrationMatrix(grid,Shat);
      else if (options.sdcSweepType == 1)
        luIntegrationMatrix(grid,Shat);
      else
        throw std::runtime_error("sdcSweepType must be 0 (Euler) or 1 (LU)");

      Matrix activeMass(expandedIndices,compressedIndex,mass);
      Matrix activeStiffness(expandedIndices,compressedIndex,stiffness);
      std::vector<Vector> massDifferences(grid.points().N(),Vector(activeDofs));
      std::vector<Vector> corrections(grid.points().N(),Vector(activeDofs));
      computeFullMassDifferences(mass,collocationStates,expandedIndices,massDifferences);

      std::vector<DiagonalMatrix> reactionDerivatives(grid.points().N(),DiagonalMatrix(activeDofs));
      for (auto& derivative : reactionDerivatives)
        derivative = 0.0;

      auto solver = [&options](Matrix const& J, Vector& du, Vector& rhs) {
        solveLinearSystem(J,du,rhs,options);
      };

      sweepNorms.push_back(sdcIterationStepJacobi(grid,Shat,solver,activeMass,activeStiffness,
                                                  residuals,reactionDerivatives,
                                                  massDifferences,corrections));

      for (int i = 1; i < grid.points().N(); ++i)
        for (size_t j = 0; j < activeDofs; ++j)
          collocationStates[i].coefficients()[expandedIndices[j]] += corrections[i][j];

      State tmpState(state);
      component<0>(tmpState) = collocationStates.back();
      Vector tmp(nDofs);
      Vector sol(nDofs);
      tmp = 0.0;
      sol = 0.0;
      tmpState.write(sol.begin());
      mass.mv(sol,tmp);
      double const normU = std::sqrt(std::max(0.0,sol*tmp));

      if (sweepNorms.size() > 1)
      {
        double const c = sweepNorms.back()/sweepNorms[sweepNorms.size()-2];
        sdcContraction = std::sqrt(c*sdcContraction);
      }

      std::cout << "  sweep " << sweep+1
                << ": ||du||=" << sweepNorms.back()
                << ", ||u||=" << normU
                << ", contraction=" << sdcContraction
                << ", active dofs=" << activeDofs << "\n";

      bool const reachedMinimumSweeps = sweep+1 >= options.minimumSdcSweeps;
      bool const smallCorrection = sweepNorms.back() < options.sdcTolerance;
      bool const reliableContraction = sdcContraction < 1.0;
      bool const estimatedSmall = reliableContraction
                               && sweepNorms.back()*sdcContraction/(1.0-sdcContraction) <= options.sdcAbsoluteTolerance;
      if (reachedMinimumSweeps && (smallCorrection || estimatedSmall))
        break;

      if (options.algebraicAdaptivity && options.algebraicAdaptivityTolerance > 0.0)
      {
        std::set<size_t> nextIndices;
        for (size_t j = 0; j < activeDofs; ++j)
        {
          double duMax = 0.0;
          for (auto const& correction : corrections)
            duMax = std::max(duMax,std::abs(correction[j][0]));

          bool const selected = sdcContraction >= 1.0
                             || sdcContraction*duMax/(1.0-sdcContraction) > options.algebraicAdaptivityTolerance;
          if (selected)
          {
            size_t const globalDof = expandedIndices[j];
            nextIndices.insert(dofNeighborhood[globalDof].begin(),dofNeighborhood[globalDof].end());
          }
        }

        if (nextIndices.empty())
        {
          std::cout << "  AA selected no active dofs for the next sweep\n";
          break;
        }

        expandedIndices.assign(nextIndices.begin(),nextIndices.end());
        compressedIndex.assign(nDofs,nDofs);
        for (size_t i = 0; i < expandedIndices.size(); ++i)
          compressedIndex[expandedIndices[i]] = i;

        std::cout << "  AA selected active dofs for next sweep: " << expandedIndices.size() << "\n";
      }

      if (grid.points().N() < options.sdcCollocationPoints+1 && sweep+1 < options.maximumSdcSweeps)
      {
        SDCTimeGrid::RealMatrix prolongation;
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

    component<0>(state) = collocationStates.back();
    F.time(stepEnd);

    if (options.writeVTK)
      writeState(state,uAll,options.order,options.outputDir + "/emiSdcStep" + paddedString(step+1,3));
  }

  return state;
}

int main(int argc, char* argv[])
{
  using namespace boost::fusion;

  EmiOptions options;
  if (getKaskadeOptions(argc,argv,Options
    ("input", options.input, options.input, "VTU mesh with domain cell data")
    ("extra_set", options.extraSet, options.extraSet, "text file with extracellular material tags")
    ("intra_set", options.intraSet, options.intraSet, "text file with intracellular material tags")
    ("excited", options.excitedSet, options.excitedSet, "text file with initially excited material tags")
    ("dir", options.outputDir, options.outputDir, "output directory")
    ("refine", options.refinements, options.refinements, "uniform mesh refinements")
    ("order", options.order, options.order, "FE polynomial order")
    ("maximumNumberOfTimeSteps", options.maximumNumberOfTimeSteps, options.maximumNumberOfTimeSteps, "maximum number of time steps; 0 means automatic from finalTime/dt")
    ("finalTime", options.finalTime, options.finalTime, "final time")
    ("dt", options.dt, options.dt, "time step size")
    ("cgTol", options.cgTol, options.cgTol, "PCG tolerance")
    ("maxCGIter", options.maxCGIter, options.maxCGIter, "maximum PCG iterations")
    ("direct", options.direct, options.direct, "use direct solver instead of PCG")
    ("CG_shift", options.cgShift, options.cgShift, "apply constant shift correction after PCG")
    ("vtk", options.writeVTK, options.writeVTK, "write VTK output: 0=no, 1=yes")
    ("diagnostics", options.diagnostics, options.diagnostics, "print per-step vector diagnostics: 0=no, 1=yes")
    ("sdc", options.sdc, options.sdc, "use SDC time integrator: 0=no, 1=yes")
    ("sdcCollocationPoints", options.sdcCollocationPoints, options.sdcCollocationPoints, "target number of SDC collocation points")
    ("sdcStartCollocationPoints", options.sdcStartCollocationPoints, options.sdcStartCollocationPoints, "initial number of SDC collocation points")
    ("minimumSdcSweeps", options.minimumSdcSweeps, options.minimumSdcSweeps, "minimum number of SDC sweeps")
    ("maximumSdcSweeps", options.maximumSdcSweeps, options.maximumSdcSweeps, "maximum number of SDC sweeps")
    ("sdcSweepType", options.sdcSweepType, options.sdcSweepType, "SDC sweep matrix: 0=Euler, 1=LU")
    ("sdcTolerance", options.sdcTolerance, options.sdcTolerance, "stop SDC when correction norm is below this value")
    ("sdcAbsoluteTolerance", options.sdcAbsoluteTolerance, options.sdcAbsoluteTolerance, "absolute SDC error-estimate tolerance")
    ("sdcInitialContraction", options.sdcInitialContraction, options.sdcInitialContraction, "initial SDC contraction estimate")
    ("algebraicAdaptivity", options.algebraicAdaptivity, options.algebraicAdaptivity, "prepare/run algebraic adaptivity mode: 0=no, 1=yes")
    ("algebraicAdaptivityTolerance", options.algebraicAdaptivityTolerance, options.algebraicAdaptivityTolerance, "AA dof-selection tolerance; 0 disables selection")
    ("bddc", options.bddc, options.bddc, "use BDDC solver for EMI: 0=no, 1=yes")
    ("bddcCompression", options.bddcCompression, options.bddcCompression, "use quantized BDDC transfer: 0=no, 1=yes")
    ("bddcCompressionBits", options.bddcCompressionBits, options.bddcCompressionBits, "quantization bits for compressed BDDC transfer")
    ("bddcIterations", options.bddcIterations, options.bddcIterations, "maximum BDDC iterations per time step")
    ("bddcTolerance", options.bddcTolerance, options.bddcTolerance, "BDDC residual tolerance")
    ("bddcInterfaceTypes", options.bddcInterfaceTypes, options.bddcInterfaceTypes, "BDDC interface flags: 1=corner, 2=edge, 4=face, 7=all")
    ("bddcVerbose", options.bddcVerbose, options.bddcVerbose, "print BDDC iteration residuals: 0=no, 1=yes")
    ("bddcUseCg", options.bddcUseCg, options.bddcUseCg, "use CG in BDDC coarse solve path: 0=no, 1=yes")
    ("nThreads", options.assemblyThreads, options.assemblyThreads, "assembler threads")
    ("penalty", options.penalty, options.penalty, "boundary penalty")
    ("sigma_i", options.sigmaI, options.sigmaI, "intracellular conductivity")
    ("sigma_e", options.sigmaE, options.sigmaE, "extracellular conductivity")
    ("C_m", options.membraneCapacitance, options.membraneCapacitance, "membrane capacitance")
    ("R", options.gapConductance, options.gapConductance, "gap junction conductance")
    ("R_extra", options.extraConductance, options.extraConductance, "extracellular interface conductance")
  )) return 0;

  if (options.dt <= 0.0)
    throw std::runtime_error("dt must be positive");
  if (options.finalTime < 0.0)
    throw std::runtime_error("finalTime must be nonnegative");
  if (options.minimumSdcSweeps < 1)
    throw std::runtime_error("minimumSdcSweeps must be at least 1");
  if (options.maximumSdcSweeps < options.minimumSdcSweeps)
    throw std::runtime_error("maximumSdcSweeps must be >= minimumSdcSweeps");
  if (options.sdcStartCollocationPoints < 1 || options.sdcCollocationPoints < options.sdcStartCollocationPoints)
    throw std::runtime_error("Require 1 <= sdcStartCollocationPoints <= sdcCollocationPoints");
  if (options.bddcIterations < 1)
    throw std::runtime_error("bddcIterations must be at least 1");
  if (options.bddcTolerance <= 0.0)
    throw std::runtime_error("bddcTolerance must be positive");

  std::filesystem::create_directories(options.outputDir);

  std::cout << "Start EMI-only model\n";
  std::cout << "mesh: " << options.input << "\n";
  std::cout << "dt: " << options.dt
            << ", finalTime: " << options.finalTime
            << ", maximumNumberOfTimeSteps: ";
  if (options.maximumNumberOfTimeSteps > 0)
    std::cout << options.maximumNumberOfTimeSteps;
  else
    std::cout << "automatic";
  std::cout << "\n";

  std::vector<int> extraTags;
  std::vector<int> intraTags;
  std::vector<int> excitedTags;
  readTagList(options.extraSet,extraTags);
  readTagList(options.intraSet,intraTags);
  readTagList(options.excitedSet,excitedTags);

  std::cout << "extracellular tags: " << extraTags.size() << "\n";
  std::cout << "intracellular tags: " << intraTags.size() << "\n";
  std::cout << "excited tags: " << excitedTags.size() << "\n";

  constexpr int dim = SPACEDIM;
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;

  boost::timer::cpu_timer totalTimer;

  VTKReader vtk(options.input);
  GridManager<Grid> gridManager(vtk.createGrid<Grid>());
  gridManager.enforceConcurrentReads(true);
  gridManager.globalRefine(options.refinements);

  using MaterialSpace = FEFunctionSpace<DiscontinuousLagrangeMapper<double,LeafView>>;
  using Material = MaterialSpace::Element<1>::type;
  MaterialSpace materialSpace(gridManager,gridManager.grid().leafGridView(),0);
  Material material(materialSpace);
  vtk.getCoefficients("domain",material);

  Dune::FieldVector<double,dim> zero(0);
  std::map<int,Dune::FieldVector<double,1>> domainMap;
  for (int tag : extraTags)
    domainMap[tag] = 0.0;
  for (int tag : intraTags)
    domainMap[tag] = static_cast<double>(tag);

  // The piecewise-continuous mapper is what makes EMI possible with one scalar
  // variable: extracellular regions share one potential, while each
  // intracellular material tag gets its own disconnected potential component.
  FEFunctionSpace uSpace(gridManager,
                         PiecewiseContinuousLagrangeMapper(gridManager.grid().leafGridView(),
                                                           options.order,
                                                           [&](auto const& cell)
                                                           {
                                                             int const tag = material.value(cell,zero);
                                                             auto const found = domainMap.find(tag);
                                                             if (found != domainMap.end())
                                                               return found->second;
                                                             Dune::FieldVector<double,1> fallback(0.0);
                                                             fallback[0] = tag;
                                                             return fallback;
                                                           }));

  using OutputSpace = L2Space<Grid>;
  OutputSpace outputSpace(gridManager,gridManager.grid().leafGridView(),options.order);
  OutputSpace::Element_t<1> uAll(outputSpace);

  auto spaces = makeSpaceList(&uSpace);
  auto variableSetDesc = makeVariableSetDescription(
    spaces,
    boost::fusion::make_vector(Variable<SpaceIndex<0>,Components<1>>("u")));

  using VariableSetDesc = decltype(variableSetDesc);
  using Membrane = AlievPanfilov;
  using Functional = EMIModel<double,VariableSetDesc,Material,Grid,decltype(spaces),Membrane>;
  using SemiLinearization = SemiLinearizationAtInner<SemiImplicitEulerStep<Functional>>;
  using Assembler = VariationalFunctionalAssembler<SemiLinearization>;
  using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
  using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;

  Functional F(material,
               gridManager.grid(),
               spaces,
               options.penalty,
               options.sigmaI,
               options.sigmaE,
               options.membraneCapacitance,
               options.gapConductance,
               options.extraConductance);
  F.extracellularMaterials(extraTags);

  auto u = variableSetDesc.variableSet();
  F.template scaleInitialValue<0>(InitialValue(material,excitedTags),u);

  size_t const nDofs = variableSetDesc.degreesOfFreedom(0,Functional::AnsatzVars::noOfVariables);
  std::cout << "cells: " << gridManager.grid().size(0) << "\n";
  std::cout << "dofs: " << nDofs << "\n";

  if (options.writeVTK)
    writeState(u,uAll,options.order,options.outputDir + "/emiInitial");

  SemiImplicitEulerStep<Functional> equation(&F,options.dt);
  equation.setTau(options.dt);

  Assembler assembler(spaces);
  auto duState(u);
  auto stepState(u);
  duState *= 0.0;
  stepState *= 0.0;

  // assemble(f, flags, nThreads): the second argument is a bit mask, not the
  // thread count. Request MATRIX explicitly, otherwise nThreads=1 assembles
  // only Assembler::VALUE and leaves a structurally nonempty but zero matrix.
  assembler.assemble(SemiLinearization(equation,u,u,duState),Assembler::MATRIX|Assembler::RHS,options.assemblyThreads);
  Matrix lhs = assembler.template get<Matrix>(false);
  printMatrixDiagnostics(lhs,"lhs");

  int const requestedNumberOfTimeSteps = static_cast<int>(std::ceil(options.finalTime/options.dt));
  // maximumNumberOfTimeSteps is optional. If it is positive, use it as a
  // debugging safety cap; otherwise run all steps requested by finalTime/dt.
  int const steps = options.maximumNumberOfTimeSteps > 0
                  ? std::min(options.maximumNumberOfTimeSteps,requestedNumberOfTimeSteps)
                  : requestedNumberOfTimeSteps;
  std::cout << "time steps: " << steps << "\n";

  if (options.sdc && options.bddc)
  {
    std::cout << "time integrator: SDC with BDDC linear solves\n";

    equation.setTau(0.0);
    F.Mass_stiff(1);
    assembler.assemble(SemiLinearization(equation,u,u,duState),
                       Assembler::MATRIX|Assembler::RHS,
                       options.assemblyThreads);
    Matrix mass = assembler.template get<Matrix>(false);

    equation.setTau(1.0);
    F.Mass_stiff(0);
    assembler.assemble(SemiLinearization(equation,u,u,duState),
                       Assembler::MATRIX|Assembler::RHS,
                       options.assemblyThreads);
    Matrix stiffness = assembler.template get<Matrix>(false);

    printMatrixDiagnostics(mass,"sdc mass");
    printMatrixDiagnostics(stiffness,"sdc stiffness");

    auto bddcData = EmiBddc::buildBddcData<Grid,decltype(uSpace),Material,Matrix,Vector>(
      gridManager.grid(),uSpace,material,lhs,nDofs);
    u = EmiBddc::runBddcSdc(F,spaces,u,uAll,bddcData,mass,stiffness,nDofs,steps,options);
    writeState(u,uAll,options.order,options.outputDir + "/emiSDCBDDCLast");
    std::cout << "total cpu-time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
    std::cout << "End EMI-only model\n";
    return 0;
  }

  if (options.sdc)
  {
    std::cout << "time integrator: SDC\n";
    u = runFullSdc(F,spaces,u,uAll,nDofs,steps,options);
    writeState(u,uAll,options.order,options.outputDir + "/emiSDCLast");
    std::cout << "total cpu-time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
    std::cout << "End EMI-only model\n";
    return 0;
  }

  if (options.bddc)
  {
    std::cout << "time integrator: semi-implicit Euler with BDDC\n";
    auto bddcData = EmiBddc::buildBddcData<Grid,decltype(uSpace),Material,Matrix,Vector>(
      gridManager.grid(),uSpace,material,lhs,nDofs);
    if (options.bddcCompression)
    {
      using CompressedTransfer = Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,std::uint16_t,std::uint8_t>;
      u = EmiBddc::runBddc<CompressedTransfer>(F,spaces,u,uAll,bddcData,nDofs,steps,options);
    }
    else
      u = EmiBddc::runBddc<Kaskade::BDDC::SpaceTransfer<1,double,double>>(F,spaces,u,uAll,bddcData,nDofs,steps,options);
    writeState(u,uAll,options.order,options.outputDir + "/emiBDDCLast");
    std::cout << "total cpu-time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
    std::cout << "End EMI-only model\n";
    return 0;
  }

  std::cout << "time integrator: semi-implicit Euler\n";
  for (int step = 0; step < steps; ++step)
  {
    double const time = step*options.dt;
    F.time(time);

    assembler.assemble(SemiLinearization(equation,u,u,duState),Assembler::RHS,options.assemblyThreads);
    auto rhs = assembler.rhs();

    Vector rhsVector(nDofs);
    Vector stepVector(nDofs);
    rhs.write(rhsVector.begin());
    if (options.diagnostics)
      printVectorDiagnostics(rhsVector,"rhs");

    std::cout << "EMI step " << step+1 << "/" << steps << "\n";
    solveLinearSystem(lhs,stepVector,rhsVector,options);

    for (size_t i = 0; i < nDofs; ++i)
      at_c<0>(stepState.data).coefficients()[i] = stepVector[i];

    component<0>(u) += component<0>(stepState);

    if (options.writeVTK)
      writeState(u,uAll,options.order,options.outputDir + "/emiStep" + paddedString(step+1,3));
  }

  writeState(u,uAll,options.order,options.outputDir + "/emiLast");

  std::cout << "total cpu-time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End EMI-only model\n";
}
