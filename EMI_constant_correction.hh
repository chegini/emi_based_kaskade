#ifndef EMI_CONSTANT_CORRECTION_HH
#define EMI_CONSTANT_CORRECTION_HH
#include <cmath>
#include <limits>

namespace EmiBddc
{
  struct ConstantCorrection
  {
    double delta=0, residualSum=0, stiffness=0;
    bool valid=false;
  };

  // Exact one-dimensional energy correction for the CURRENT linear system:
  // delta = (1^T (b-A*x))/(1^T A*1). No mean of x or target solution is used.
  // Wider accumulators reduce cancellation in the weakly anchored mode.
  template<class Matrix,class Vector>
  ConstantCorrection proposeConstantCorrection(Matrix const& A, Vector const& b,
                                                Vector const& x)
  {
    long double numerator=0,denominator=0;
    for (auto row=A.begin(); row!=A.end(); ++row)
    {
      long double residual=b[row.index()][0];
      for (auto col=row.begin(); col!=row.end(); ++col)
      {
        long double const value=(*col)[0][0];
        residual-=value*x[col.index()][0];
        denominator+=value;
      }
      numerator+=residual;
    }
    ConstantCorrection result;
    result.residualSum=static_cast<double>(numerator);
    result.stiffness=static_cast<double>(denominator);
    if (!(denominator>0) || !std::isfinite(denominator) || !std::isfinite(numerator))
      return result;
    result.delta=static_cast<double>(numerator/denominator);
    result.valid=std::isfinite(result.delta) && std::isfinite(result.residualSum)
                 && std::isfinite(result.stiffness) && result.stiffness>0;
    return result;
  }
}
#endif
