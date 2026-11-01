#ifndef EMI_FGMRES_HH
#define EMI_FGMRES_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace EmiBddc::Fgmres
{
  using Vector = std::vector<double>;
  enum class Status { converged, iterationLimit, breakdown, nonfinite };
  struct Result
  {
    Vector x;
    Status status = Status::iterationLimit;
    int iterations = 0, cycles = 0, matvecs = 0;
    double residual = 0, rhsNorm = 0;
    std::size_t peakBasisBytes = 0;
  };

  inline double dot(Vector const& a, Vector const& b)
  {
    if (a.size()!=b.size()) throw std::invalid_argument("FGMRES vector size mismatch");
    long double sum=0;
    for (std::size_t i=0; i<a.size(); ++i) sum+=static_cast<long double>(a[i])*b[i];
    return static_cast<double>(sum);
  }
  inline double norm(Vector const& v)
  {
    double sum=0;
    for (double x:v) sum=std::hypot(sum,x);
    return sum;
  }
  inline bool finite(Vector const& v)
  {
    return std::all_of(v.begin(),v.end(),[](double x) {return std::isfinite(x);});
  }
  inline void axpy(Vector& y, double a, Vector const& x)
  {
    if (y.size()!=x.size()) throw std::invalid_argument("FGMRES vector size mismatch");
    for (std::size_t i=0; i<y.size(); ++i) y[i]+=a*x[i];
  }

  /** Restarted flexible RIGHT-preconditioned GMRES (Saad, 1993).
   * z_j=M_j(v_j), w=A*z_j; store each z_j and update x=x_base+Z*y.
   * Two-pass modified Gram-Schmidt, incremental Givens QR, explicit b-A*x
   * after EVERY candidate. The true norm alone accepts convergence; a small
   * Hessenberg residual or happy breakdown is not sufficient. No rho sign test.
   * applyA and applyM return new vectors. observe sees (iteration,cycle,
   * estimated absolute residual,true absolute residual). Runs on one process;
   * an application adapter supplies the distributed preconditioner.
   */
  template<class ApplyA, class ApplyM, class Observe>
  Result solve(Vector const& b, int restart, int limit, double tolerance,
               ApplyA applyA, ApplyM applyM, Observe observe)
  {
    if (b.empty() || restart<1 || limit<1 || !(tolerance>0) || !std::isfinite(tolerance))
      throw std::invalid_argument("Invalid FGMRES dimensions/limits/tolerance");
    Result out; out.x.assign(b.size(),0); out.rhsNorm=norm(b); out.residual=out.rhsNorm;
    auto meets=[&]() {return out.residual<=tolerance*(out.rhsNorm>0 ? out.rhsNorm : 1);};
    if (!finite(b) || !std::isfinite(out.rhsNorm)) {out.status=Status::nonfinite; return out;}
    observe(0,0,out.residual,out.residual);
    if (meets()) {out.status=Status::converged; return out;}
    Vector r=b;
    while (out.iterations<limit)
    {
      ++out.cycles;
      int const m=std::min(restart,limit-out.iterations);
      Vector base=out.x;
      std::vector<Vector> v, z;
      v.push_back(r);
      for (double& x:v.back()) x/=out.residual;
      std::vector<Vector> h(m,Vector(m+1,0));
      Vector cs(m),sn(m),g(m+1,0); g[0]=out.residual;
      for (int j=0; j<m; ++j)
      {
        z.push_back(applyM(v[j])); ++out.iterations;
        if (z.back().size()!=b.size()) throw std::invalid_argument("FGMRES preconditioner size");
        if (!finite(z.back())) {out.status=Status::nonfinite;return out;}
        Vector w=applyA(z.back()); ++out.matvecs;
        if (w.size()!=b.size()) throw std::invalid_argument("FGMRES operator size");
        if (!finite(w)) {out.status=Status::nonfinite;return out;}
        double const originalNorm=norm(w);
        for (int pass=0; pass<2; ++pass)
          for (int i=0; i<=j; ++i)
          {
            double const a=dot(v[i],w); h[j][i]+=a; axpy(w,-a,v[i]);
          }
        double const nextNorm=norm(w);
        bool const exhausted=nextNorm<=32*std::numeric_limits<double>::epsilon()*originalNorm;
        h[j][j+1]=nextNorm;
        out.peakBasisBytes=std::max(out.peakBasisBytes,(v.size()+z.size())*b.size()*sizeof(double));
        for (int i=0; i<j; ++i)
        {
          double const a=h[j][i], c=h[j][i+1];
          h[j][i]=cs[i]*a+sn[i]*c; h[j][i+1]=-sn[i]*a+cs[i]*c;
        }
        double const diagonal=std::hypot(h[j][j],h[j][j+1]);
        if (!std::isfinite(diagonal)) {out.status=Status::nonfinite;return out;}
        if (diagonal==0) {out.status=Status::breakdown;return out;}
        cs[j]=h[j][j]/diagonal; sn[j]=h[j][j+1]/diagonal;
        h[j][j]=diagonal; h[j][j+1]=0;
        g[j+1]=-sn[j]*g[j]; g[j]*=cs[j];
        Vector y(g.begin(),g.begin()+j+1);
        for (int i=j; i>=0; --i)
        {
          for (int k=i+1; k<=j; ++k) y[i]-=h[k][i]*y[k];
          if (h[i][i]==0) {out.status=Status::breakdown;return out;}
          y[i]/=h[i][i];
        }
        Vector trial=base;
        for (int i=0; i<=j; ++i) axpy(trial,y[i],z[i]);
        if (!finite(trial)) {out.status=Status::nonfinite;return out;}
        Vector residual=b;
        axpy(residual,-1,applyA(trial)); ++out.matvecs;
        double const trueNorm=norm(residual);
        if (!std::isfinite(trueNorm)) {out.status=Status::nonfinite;return out;}
        out.x=std::move(trial); r=std::move(residual); out.residual=trueNorm;
        observe(out.iterations,out.cycles,std::abs(g[j+1]),trueNorm);
        if (meets()) {out.status=Status::converged;return out;}
        if (exhausted) {out.status=Status::breakdown;return out;}
        if (j+1<m)
        {
          for (double& x:w) x/=nextNorm;
          v.push_back(std::move(w));
        }
      }
    }
    return out;
  }
}
#endif
