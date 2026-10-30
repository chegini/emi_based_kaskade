#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "../EMI_bddc_convergence.hh"

using Vector = std::array<double,2>;

double dot(Vector const& a, Vector const& b) { return a[0]*b[0]+a[1]*b[1]; }
Vector multiply(Vector const& x) { return {4*x[0]+x[1],x[0]+3*x[1]}; }

struct Run
{
  Vector x;
  EmiBddc::TrueResidualResult result;
  int checks;
};

// Solve a genuinely coupled SPD system with diagonal PCG. Reporting observes
// precisely the same accepted corrections and checks as the unreported run.
Run solve(double scale, std::string const& csv = {})
{
  Vector const b{scale,2*scale};
  Vector x{}, r=b, z{r[0]/4,r[1]/3}, p=z;
  double rho = dot(r,z);
  int updates=0, checks=0;
  std::optional<double> alpha;
  EmiBddc::ConvergenceReport report;
  std::ostringstream output;
  output.precision(12);
  auto const flags = output.flags();
  if (!csv.empty())
    assert(report.open(csv,1,1e-12,output));
  auto result = EmiBddc::solveToTrueResidual(10,1e-12,
    [&]()
    {
      Vector const ap=multiply(p);
      alpha=rho/dot(p,ap);
      for (int i=0; i<2; ++i)
      {
        x[i] += *alpha*p[i];
        r[i] -= *alpha*ap[i];
      }
      z={r[0]/4,r[1]/3};
      double const nextRho=dot(r,z);
      for (int i=0; i<2; ++i)
        p[i]=z[i]+(nextRho/rho)*p[i];
      rho=nextRho;
      ++updates;
      return 0.0; // The actual residual, not this internal norm, must stop CG.
    }, []() { return false; },
    [&]()
    {
      ++checks;
      Vector const ax=multiply(x);
      Vector const residual{b[0]-ax[0],b[1]-ax[1]};
      EmiBddc::TrueResidual measured{std::sqrt(dot(residual,residual)),std::sqrt(dot(b,b))};
      if (!csv.empty())
        report.record(updates,alpha,measured);
      return measured;
    });
  if (!csv.empty())
  {
    assert(report.close());
    assert(output.str().find("--")!=std::string::npos);
    assert(output.precision()==12 && output.flags()==flags);
    std::ifstream input(csv);
    std::string line;
    int lines=0;
    while (std::getline(input,line)) ++lines;
    assert(lines==result.calls+2); // Header, initial row, accepted updates.
  }
  return {x,result,checks};
}

int main(int argc, char** argv)
{
  assert(argc==2); // Caller supplies a temporary output directory.
  auto const directory=std::filesystem::path(argv[1]);
  std::filesystem::create_directories(directory);
  int index=0;
  for (double scale : {1e-20,1.0,1e20})
  {
    Run const baseline=solve(scale);
    Run const reported=solve(scale,(directory/("pcg-"+std::to_string(++index)+".csv")).string());
    assert(baseline.result.converged() && reported.result.converged());
    assert(baseline.x==reported.x);
    assert(baseline.checks==reported.checks && reported.checks==reported.result.calls+1);
    assert(baseline.result.calls==reported.result.calls && reported.result.calls==2);
    assert(baseline.result.measured.absolute==reported.result.measured.absolute);
    assert(std::abs(reported.x[0]/scale-1.0/11)<1e-12);
    assert(std::abs(reported.x[1]/scale-7.0/11)<1e-12);
  }

  std::ostringstream output;
  EmiBddc::ConvergenceReport zero;
  assert(zero.open((directory/"zero.csv").string(),2,1e-7,output));
  auto result=EmiBddc::solveToTrueResidual(2,1e-7,
    []() { assert(false); return 0.0; }, []() { return false; },
    [&]() { zero.record(0,std::nullopt,{0.0,0.0}); return EmiBddc::TrueResidual{0.0,0.0}; });
  assert(result.converged() && result.calls==0 && zero.close());

  EmiBddc::ConvergenceReport rejected;
  assert(rejected.open((directory/"rejected.csv").string(),3,1e-7,output));
  int checks=0;
  result=EmiBddc::solveToTrueResidual(2,1e-7,
    []() { return 0.0; }, []() { return true; },
    [&]() { ++checks; rejected.record(0,std::nullopt,{1.0,1.0}); return EmiBddc::TrueResidual{1.0,1.0}; });
  assert(result.status==EmiBddc::TrueResidualStatus::solverStopped && checks==1);
  assert(rejected.close()); // A rejected correction must not create a new row.

  EmiBddc::ConvergenceReport missing;
  assert(!missing.open((directory/"missing-parent"/"data.csv").string(),1,1e-7,output));
  std::cout << "BDDC convergence report checks passed; reporting preserves PCG results and check counts\n";
}
