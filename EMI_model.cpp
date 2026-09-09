#define FUSION_MAX_VECTOR_SIZE 20

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/timer/timer.hpp>

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
#include "utilities/kaskopt.hh"
#include "utilities/timing.hh"

#include "fem/diffops/membraneModels.hh"

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

  double finalTime = 0.01;
  double dt = 0.01;
  double cgTol = 1e-8;
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
