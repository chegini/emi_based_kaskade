#ifndef EMI_LOCAL_SUPPORT_DIAGNOSTIC_HH
#define EMI_LOCAL_SUPPORT_DIAGNOSTIC_HH

#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace EmiBddc
{
  struct MaterialLocalSupport
  {
    size_t omittedEntries = 0;
    double retainedMaxAbs = 0.0;
    double omittedMaxAbs = 0.0;
    long double omittedAbsSumBound = 0.0L;

    long double roundoffBudget() const
    {
      return static_cast<long double>(std::numeric_limits<double>::epsilon())
             * retainedMaxAbs;
    }
    bool acceptable() const
    {
      return omittedEntries == 0 ||
        (retainedMaxAbs > 0.0 && omittedAbsSumBound <= roundoffBudget());
    }
  };

  // Matrix extraction retains every entry with both endpoints in localDofs,
  // including arbitrarily small entries. Bound only what extraction omits.
  // The absolute sum bounds the Frobenius and spectral norms of that change.
  // No absolute floor: a uniformly rescaled problem gets the same decision.
  // Round the positive sum upwards so individually small entries cannot hide
  // a large aggregate defect, even through accumulation rounding.
  template <class Matrix>
  MaterialLocalSupport inspectMaterialLocalSupport(Matrix const& matrix,
      std::set<size_t> const& localDofs)
  {
    MaterialLocalSupport result;
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const value = (*col)[0][0];
        if (!std::isfinite(value))
          throw std::runtime_error("Nonfinite BDDC material contribution entry.");
        if (value == 0.0) continue;
        double const magnitude = std::abs(value);
        if (localDofs.contains(row.index()) && localDofs.contains(col.index()))
        {
          if (magnitude > result.retainedMaxAbs)
            result.retainedMaxAbs = magnitude;
        }
        else
        {
          ++result.omittedEntries;
          if (magnitude > result.omittedMaxAbs)
            result.omittedMaxAbs = magnitude;
          result.omittedAbsSumBound = std::nextafter(
            result.omittedAbsSumBound + static_cast<long double>(magnitude),
            std::numeric_limits<long double>::infinity());
        }
      }
    return result;
  }

  // Failure details; this function never modifies entries or local spaces.
  template <class Matrix>
  std::string describeMissingMaterialSupport(Matrix const& matrix,
      std::set<size_t> const& localDofs, size_t subdomain, int tag)
  {
    size_t count = 0, largestRow = 0, largestCol = 0;
    double largestValue = 0.0, maxEntry = 0.0;
    std::set<size_t> missingDofs;
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const value = (*col)[0][0];
        if (std::abs(value) > maxEntry)
          maxEntry = std::abs(value);
        if (value == 0.0 ||
            (localDofs.contains(row.index()) && localDofs.contains(col.index())))
          continue;
        if (count == 0 || std::abs(value) > std::abs(largestValue))
        {
          largestRow = row.index();
          largestCol = col.index();
          largestValue = value;
        }
        ++count;
        if (!localDofs.contains(row.index())) missingDofs.insert(row.index());
        if (!localDofs.contains(col.index())) missingDofs.insert(col.index());
      }
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "legacy BDDC local space omits a nonzero material coupling endpoint."
            << " subdomain_index=" << subdomain << " material_tag=" << tag
            << " local_dofs=" << localDofs.size()
            << " missing_entries=" << count << " missing_dofs=" << missingDofs.size()
            << " largest_missing_row=" << largestRow
            << " largest_missing_col=" << largestCol
            << " largest_missing_value=" << largestValue
            << " row_present=" << localDofs.contains(largestRow)
            << " col_present=" << localDofs.contains(largestCol)
            << " contribution_max_abs=" << maxEntry
            << " missing_global_dofs_first16=";
    size_t shown = 0;
    for (size_t dof : missingDofs)
    {
      if (shown == 16) { message << "..."; break; }
      message << (shown++ ? "," : "") << dof;
    }
    return message.str();
  }

  template <class Matrix>
  MaterialLocalSupport checkMaterialLocalSupport(Matrix const& matrix,
      std::set<size_t> const& localDofs, size_t subdomain, int tag)
  {
    auto const support = inspectMaterialLocalSupport(matrix,localDofs);
    if (!support.acceptable())
    {
      std::ostringstream message;
      message << describeMissingMaterialSupport(matrix,localDofs,subdomain,tag)
              << std::setprecision(std::numeric_limits<long double>::max_digits10)
              << " omitted_abs_sum_bound=" << support.omittedAbsSumBound
              << " roundoff_budget=" << support.roundoffBudget();
      throw std::runtime_error(message.str());
    }
    return support;
  }
}
#endif
