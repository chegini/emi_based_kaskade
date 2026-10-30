#ifndef EMI_BDDC_STOPPING_HH
#define EMI_BDDC_STOPPING_HH

#include <cmath>
#include <limits>
#include <stdexcept>

namespace EmiBddc
{
  // Euclidean norms of the assembled, unpreconditioned system residual and RHS.
  // With a zero initial global correction, ||b|| is the initial residual norm.
  struct TrueResidual
  {
    double absolute = std::numeric_limits<double>::infinity();
    double rhsNorm = 0.0;

    bool valid() const
    {
      return std::isfinite(absolute) && absolute >= 0.0
          && std::isfinite(rhsNorm) && rhsNorm >= 0.0;
    }

    double relative() const
    {
      // A zero RHS has no relative scale. Use the absolute residual in that case.
      return rhsNorm > 0.0 ? absolute/rhsNorm : absolute;
    }

    bool meets(double tolerance) const
    {
      return valid() && relative() <= tolerance;
    }
  };

  enum class TrueResidualStatus { converged, iterationLimit, solverStopped, nonfiniteResidual };

  struct TrueResidualResult
  {
    TrueResidual measured;
    TrueResidualStatus status = TrueResidualStatus::iterationLimit;
    int calls = 0;
    double internalNorm = std::numeric_limits<double>::infinity();

    bool converged() const { return status == TrueResidualStatus::converged; }
  };

  // Check the initial accepted iterate (including any interior elimination),
  // then every accepted solver update. The internal search-direction norm
  // never terminates this mode. A solver rejection leaves the previous iterate
  // unchanged, so its cached true residual remains the relevant measurement.
  // Check must produce the same measurement/decision on every MPI rank.
  template <class Solve, class Stopped, class Check>
  TrueResidualResult solveToTrueResidual(int limit, double tolerance,
                                          Solve solve, Stopped stopped, Check check)
  {
    if (limit < 1 || !std::isfinite(tolerance) || tolerance <= 0.0)
      throw std::invalid_argument("True-residual stopping needs a positive limit and tolerance");
    TrueResidualResult result;
    auto inspect = [&]()
    {
      result.measured = check();
      if (!result.measured.valid())
        result.status = TrueResidualStatus::nonfiniteResidual;
      else if (result.measured.meets(tolerance))
        result.status = TrueResidualStatus::converged;
      else
        return false;
      return true;
    };
    if (inspect())
      return result;
    for (int i = 0; i < limit; ++i)
    {
      result.internalNorm = solve();
      ++result.calls;
      if (stopped())
      {
        result.status = TrueResidualStatus::solverStopped;
        return result;
      }
      if (inspect())
        return result;
    }
    result.status = TrueResidualStatus::iterationLimit;
    return result;
  }
}

#endif
