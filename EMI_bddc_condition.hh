#ifndef EMI_BDDC_CONDITION_HH
#define EMI_BDDC_CONDITION_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace EmiBddc
{
  struct CgConditionEstimate
  {
    double lambdaMin = 0;
    double lambdaMax = 0;
    double condition = 0;
    std::size_t dimension = 0;
    char const* status = "insufficient-iterations";
  };

  /**
   * Reconstruct the Lanczos projection from accepted PCG steps only:
   * d[0]=1/alpha[0], d[i]=1/alpha[i]+beta[i]/alpha[i-1],
   * e[i-1]=sqrt(beta[i])/alpha[i-1]. beta[0] must be zero.
   *
   * For fixed SPD A and B, these are Ritz values of B^(1/2) A B^(1/2),
   * restricted to the Krylov space after the solver's interior elimination.
   * Their ratio is an RHS-dependent spectral condition estimate, not the
   * Euclidean singular-value condition of A or BA. Unseen modes, finite
   * precision and loss of CG orthogonality can make it underestimate the
   * full spectral range. No SPD proof is supplied by this diagnostic.
   * See Netlib Templates, "The Lanczos Connection":
   * https://www.netlib.org/linalg/html_templates/node118.html
   */
  class CgConditionEstimator
  {
  public:
    void append(double alpha, double beta)
    {
      if (!valid_)
        return;
      if (!std::isfinite(alpha) || alpha <= 0 || !std::isfinite(beta) || beta < 0
          || (diagonal_.empty() && beta != 0))
      {
        valid_ = false;
        return;
      }
      long double diagonal = 1.0L/alpha;
      if (!diagonal_.empty())
      {
        diagonal += static_cast<long double>(beta)/previousAlpha_;
        offDiagonal_.push_back(std::sqrt(static_cast<long double>(beta))/previousAlpha_);
      }
      diagonal_.push_back(diagonal);
      previousAlpha_ = alpha;
    }

    CgConditionEstimate estimate() const
    {
      CgConditionEstimate result;
      result.dimension = diagonal_.size();
      if (!valid_)
      {
        result.status = "invalid-cg-scalars";
        return result;
      }
      if (diagonal_.size() < 2)
        return result;

      // Scale the small symmetric tridiagonal, then use a Sturm LDL count
      // and bisection for its two extreme eigenvalues. This needs O(k)
      // storage, no dense global matrix, eigenvectors, MPI or extra solves.
      long double scale = 0;
      for (auto d : diagonal_)
        scale = std::max(scale,std::abs(d));
      for (auto e : offDiagonal_)
        scale = std::max(scale,std::abs(e));
      if (!std::isfinite(scale) || scale <= 0)
      {
        result.status = "invalid-tridiagonal";
        return result;
      }
      auto d = diagonal_;
      auto e = offDiagonal_;
      for (auto& value : d) value /= scale;
      for (auto& value : e) value /= scale;

      long double lower = std::numeric_limits<long double>::infinity();
      long double upper = -lower;
      for (std::size_t i = 0; i < d.size(); ++i)
      {
        long double radius = (i ? std::abs(e[i-1]) : 0)
                           + (i+1 < d.size() ? std::abs(e[i]) : 0);
        lower = std::min(lower,d[i]-radius);
        upper = std::max(upper,d[i]+radius);
      }
      long double const margin = 16*std::numeric_limits<long double>::epsilon()
                               * std::max(std::abs(lower),std::abs(upper));
      lower -= margin;
      upper += margin;
      auto countBelow = [&](long double x)
      {
        // A tiny negative pivot implements the limiting Sturm convention
        // at an exact zero while avoiding division by zero and overflow.
        long double const pivotMin = 16*std::numeric_limits<long double>::min();
        long double pivot = d[0]-x;
        std::size_t count = 0;
        for (std::size_t i = 0; i < d.size(); ++i)
        {
          if (i) pivot = d[i]-x-e[i-1]*(e[i-1]/pivot);
          if (std::abs(pivot) < pivotMin) pivot = -pivotMin;
          if (pivot < 0) ++count;
        }
        return count;
      };
      auto eigenvalue = [&](std::size_t index)
      {
        long double lo = lower, hi = upper;
        for (int iteration = 0; iteration < 256; ++iteration)
        {
          long double mid = lo+(hi-lo)/2;
          if (mid == lo || mid == hi
              || hi-lo <= 16*std::numeric_limits<long double>::epsilon()
                             * std::max(std::abs(lo),std::abs(hi)))
            return mid;
          if (countBelow(mid) <= index) lo = mid;
          else hi = mid;
        }
        return std::numeric_limits<long double>::quiet_NaN();
      };
      long double const minimum = eigenvalue(0);
      long double const maximum = eigenvalue(d.size()-1);
      if (!std::isfinite(minimum) || !std::isfinite(maximum))
      {
        result.status = "eigensolve-unresolved";
        return result;
      }
      if (minimum <= 0 || maximum < minimum)
      {
        result.status = "nonpositive-ritz-spectrum";
        return result;
      }
      result.lambdaMin = static_cast<double>(minimum*scale);
      result.lambdaMax = static_cast<double>(maximum*scale);
      result.condition = static_cast<double>(maximum/minimum);
      if (!std::isfinite(result.lambdaMin) || result.lambdaMin <= 0
          || !std::isfinite(result.lambdaMax) || !std::isfinite(result.condition))
      {
        result.status = "unrepresentable-spectrum";
        return result;
      }
      result.status = "estimated";
      return result;
    }

  private:
    bool valid_ = true;
    long double previousAlpha_ = 0;
    std::vector<long double> diagonal_;
    std::vector<long double> offDiagonal_;
  };
}

#endif
