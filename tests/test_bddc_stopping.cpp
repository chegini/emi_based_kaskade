#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "../EMI_bddc_stopping.hh"

int main()
{
  using namespace EmiBddc;
  // A small internal search norm must not stop a solve with a large true residual.
  std::vector<double> residuals = {1.0,1e-3,1e-6,1e-8};
  int updates = 0, checks = 0;
  auto result = solveToTrueResidual(10,1e-7,
    [&]() { ++updates; return 0.0; },
    []() { return false; },
    [&]() { ++checks; return TrueResidual{residuals.at(updates),1.0}; });
  assert(result.converged() && result.calls == 3 && checks == 4);
  assert(result.measured.relative() == residuals.back());

  // Scaling the entire system's RHS/residual leaves the relative decision unchanged.
  for (double scale : {1e-20,1.0,1e20})
  {
    updates = 0;
    result = solveToTrueResidual(10,1e-7,
      [&]() { ++updates; return 1.0; },
      []() { return false; },
      [&]() { return TrueResidual{scale*residuals.at(updates),scale}; });
    assert(result.converged() && result.calls == 3);
  }

  // An initially solved system, including zero RHS, needs no CG calls.
  updates = 0;
  result = solveToTrueResidual(10,1e-7,
    [&]() { ++updates; return 1.0; },
    []() { return false; },
    []() { return TrueResidual{0.0,0.0}; });
  assert(result.converged() && updates == 0 && result.calls == 0);
  assert((TrueResidual{0.5,1.0}.meets(0.5))); // equality is accepted
  assert((TrueResidual{1e-9,0.0}.meets(1e-7))); // documented zero-RHS absolute fallback
  assert((!TrueResidual{1e-6,0.0}.meets(1e-7)));

  // Never label an exhausted iteration budget as convergence.
  result = solveToTrueResidual(2,1e-7,
    []() { return 0.0; },
    []() { return false; },
    []() { return TrueResidual{1.0,1.0}; });
  assert(!result.converged() && result.calls == 2);
  assert(result.status == TrueResidualStatus::iterationLimit);

  // A rejected correction has not changed the accepted solution: retain its
  // measured residual, without inspecting a hypothetical rejected update.
  checks = 0;
  result = solveToTrueResidual(10,1e-7,
    []() { return std::numeric_limits<double>::infinity(); },
    []() { return true; },
    [&]() { ++checks; return TrueResidual{1.0,1.0}; });
  assert(result.status == TrueResidualStatus::solverStopped && result.calls == 1);
  assert(!result.converged() && checks == 1 && result.measured.absolute == 1.0);

  // NaN/Inf and invalid norms cannot silently satisfy the target.
  for (auto measured : {TrueResidual{std::numeric_limits<double>::quiet_NaN(),1.0},
                        TrueResidual{1.0,std::numeric_limits<double>::infinity()},
                        TrueResidual{-1.0,1.0},
                        TrueResidual{1.0,-1.0}})
  {
    result = solveToTrueResidual(10,1e-7,
      []() { return 0.0; },
      []() { return false; },
      [=]() { return measured; });
    assert(result.status == TrueResidualStatus::nonfiniteResidual && result.calls == 0);
  }
  updates = 0;
  result = solveToTrueResidual(10,1e-7,
    [&]() { ++updates; return 0.0; },
    []() { return false; },
    [&]() { return TrueResidual{updates ? std::numeric_limits<double>::quiet_NaN() : 1.0,1.0}; });
  assert(result.status == TrueResidualStatus::nonfiniteResidual && result.calls == 1);

  for (double tolerance : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()})
  {
    bool rejected = false;
    try
    {
      (void)solveToTrueResidual(10,tolerance,[]() { return 0.0; },
                               []() { return false; },[]() { return TrueResidual{0.0,1.0}; });
    }
    catch (std::invalid_argument const&)
    {
      rejected = true;
    }
    assert(rejected);
  }
  std::cout << "BDDC true-residual stopping checks passed\n";
}
