#ifndef EMI_DIRECT_MATRIX_HH
#define EMI_DIRECT_MATRIX_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace EmiBddc
{
  struct SpdInputCheck
  {
    double maxAbs = 0;
    double maxAbsSkew = 0;
    double symmetryTolerance = 0;
    std::size_t lowerNonzeros = 0;
  };

  // Validate the full assembled scalar-block matrix before a solver uses only
  // its lower triangle. Rows must have ascending column indices, as in NumaBCRS.
  // This is not an SPD proof: MUMPS must still accept the positive-definite
  // factorization. No diagonal shift, dropping tolerance or symmetrization is
  // applied. The symmetry tolerance permits assembly-level roundoff only.
  template <class Matrix>
  SpdInputCheck checkSpdDirectInput(Matrix const& matrix)
  {
    if (matrix.N() == 0 || matrix.N() != matrix.M())
      throw std::runtime_error("SPD direct input must be nonempty and square");
    if (static_cast<std::size_t>(matrix.N()) >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::runtime_error("SPD direct input exceeds the MUMPS integer index range");
    SpdInputCheck check;
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
    {
      bool positiveDiagonal = false;
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const value = (*col)[0][0];
        if (!std::isfinite(value))
          throw std::runtime_error("Nonfinite SPD direct matrix entry");
        if (col.index() >= matrix.N())
          throw std::runtime_error("Invalid SPD direct matrix column index");
        check.maxAbs = std::max(check.maxAbs,std::abs(value));
        if (row.index() == col.index()) positiveDiagonal = value > 0;
        if (row.index() >= col.index() && value != 0) ++check.lowerNonzeros;

        auto const transposeRow = matrix[col.index()];
        auto transpose = transposeRow.begin();
        while (transpose != transposeRow.end() && transpose.index() < row.index())
          ++transpose;
        double const transposed = transpose != transposeRow.end()
                               && transpose.index() == row.index() ? (*transpose)[0][0] : 0;
        if (!std::isfinite(transposed))
          throw std::runtime_error("Nonfinite transposed SPD direct matrix entry");
        check.maxAbsSkew = std::max(check.maxAbsSkew,std::abs(value-transposed));
      }
      if (!positiveDiagonal)
        throw std::runtime_error("SPD direct input needs a positive diagonal in every row");
    }
    check.symmetryTolerance = 1e-12*check.maxAbs;
    if (check.maxAbsSkew > check.symmetryTolerance)
      throw std::runtime_error("SPD direct input is not symmetric within assembly roundoff");
    if (check.lowerNonzeros > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::runtime_error("SPD direct nonzeros exceed the MUMPS integer index range");
    return check;
  }

  // Only exact zeros and the redundant upper triangle are omitted. MUMPS sums
  // duplicate symmetric entries, so sending both nonzero halves would double
  // off-diagonals and change the operator. Scalar indices are checked above.
  template <class Matrix, class Visitor>
  void visitDirectLowerNonzeros(Matrix const& matrix, Visitor&& emit)
  {
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
      for (auto col = row.begin(); col != row.end(); ++col)
        if (row.index() >= col.index() && (*col)[0][0] != 0)
          emit(static_cast<int>(row.index()),static_cast<int>(col.index()),(*col)[0][0]);
  }
}

#endif
