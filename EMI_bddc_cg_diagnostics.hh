#ifndef EMI_BDDC_CG_DIAGNOSTICS_HH
#define EMI_BDDC_CG_DIAGNOSTICS_HH

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace EmiBddc
{
  namespace CgDiagnostics
  {
    inline double norm(std::vector<double> const& x)
    {
      double result = 0;
      for (double value : x) result = std::hypot(result,value);
      return result;
    }

    inline double dot(std::vector<double> const& x, std::vector<double> const& y)
    {
      double result = 0;
      for (size_t i=0; i<x.size(); ++i) result += x[i]*y[i];
      return result;
    }

    inline double missing() { return std::numeric_limits<double>::quiet_NaN(); }

    // An undefined diagnostic must not look like perfect agreement.
    inline double ratio(double numerator, double denominator)
    {
      return std::isfinite(numerator) && std::isfinite(denominator) && denominator>0
           ? numerator/denominator : missing();
    }
  }

  /**
   * Read-only EMI CG observations on the original assembled system.
   *
   * Each owned subdomain contributes once. Primal vectors are assembled with
   * the application's partition-of-unity weights; restricted residuals are
   * summed WITHOUT another weight, since restriction already distributes the
   * dual residual. Norms/dots are evaluated once on root in global DOF space.
   * Shared-copy jumps are reported so projection cannot hide inconsistency.
   *
   * Two diagnostic sum reductions per event and root matrix-vector products
   * are additional work. They never call transfer methods or change solver
   * vectors/scalars, and must not be treated as compression timing data.
   */
  template<class Matrix, class Vector, class Storage, class Data>
  class CgDiagnosticReport
  {
  public:
    using Sum = std::function<void(std::vector<double>&)>;

    CgDiagnosticReport(Matrix const& matrix, Vector const& rhs, Storage const& storage,
                       Data const& data, bool root, Sum sum)
      : matrix_(matrix), rhs_(rhs), storage_(storage), data_(data),
        root_(root), sum_(std::move(sum)), n_(rhs.size())
    {}

    bool open(std::string const& path)
    {
      if (!root_) return true;
      output_.open(path);
      output_ << "stage,solver_call,restriction_pass,applied,true_residual_norm,"
                 "restriction_error_abs,restriction_error_rel,solution_change_since_restriction,"
                 "rho_internal,rho_true,rho_error_abs,rho_error_rel,rho_true_near_zero,"
                 "beta,alpha_candidate,curvature_internal,curvature_global,"
                 "conjugacy_previous,solution_jump,z_jump,q_jump\n";
      output_ << std::setprecision(17);
      return output_.good();
    }

    template<class Event>
    void observe(Event const& event)
    {
      // All ranks, including ranks owning no subdomains, follow these reductions.
      std::vector<double> values(4*n_,0), localU;
      for (size_t s=0; s<storage_.size(); ++s)
      {
        if (!storage_.has(s)) continue;
        auto const& sub = storage_[s];
        auto const x = sub.getSolution();
        auto const& rr = sub.getRestrictedResidual();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
        {
          size_t const g = data_.localDofs[s][i];
          double const w = data_.weights[s][i][0];
          localU.push_back(x[i][0]);
          values[g] += w*x[i][0];
          values[n_+g] += rr[i][0];
          if (!event.restriction)
          {
            values[2*n_+g] += w*sub.diagnosticPreconditionedResidual()[i][0];
            if (event.directionAvailable)
              values[3*n_+g] += w*sub.diagnosticSearchDirection()[i][0];
          }
        }
      }
      sum_(values);
      std::vector<double> jumps(4,0);
      for (size_t s=0; s<storage_.size(); ++s)
      {
        if (!storage_.has(s)) continue;
        auto const& sub = storage_[s];
        auto const x = sub.getSolution();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
        {
          size_t const g = data_.localDofs[s][i];
          double d = x[i][0]-values[g];
          jumps[0] += d*d;
          if (!event.restriction)
          {
            d = sub.diagnosticPreconditionedResidual()[i][0]-values[2*n_+g];
            jumps[1] += d*d;
            if (event.directionAvailable)
            {
              d = sub.diagnosticSearchDirection()[i][0]-values[3*n_+g];
              jumps[2] += d*d;
            }
          }
        }
      }
      if (call_ == event.solverCall && previousLocalU_.size()==localU.size())
        for (size_t i=0; i<localU.size(); ++i)
        {
          double const d = localU[i]-previousLocalU_[i];
          jumps[3] += d*d;
        }
      call_ = event.solverCall;
      if (event.restriction) previousLocalU_ = localU;
      sum_(jumps);
      if (!root_) return;

      using namespace CgDiagnostics;
      std::vector<double> x(values.begin(),values.begin()+n_);
      std::vector<double> rr(values.begin()+n_,values.begin()+2*n_);
      std::vector<double> z(values.begin()+2*n_,values.begin()+3*n_);
      std::vector<double> q(values.begin()+3*n_,values.end());
      // The solution is fixed during restriction refinement and direction
      // construction. Reuse A*x within a call, unless an unexpected change is
      // observed (which is also recorded in the solution-change diagnostic).
      if (referenceCall_!=event.solverCall || referenceX_!=x)
      {
        auto Ax = multiply(x);
        referenceResidual_.resize(n_);
        for (size_t i=0; i<n_; ++i)
          referenceResidual_[i] = rhs_[i][0]-Ax[i];
        referenceX_ = x;
        referenceCall_ = event.solverCall;
      }
      auto const& r = referenceResidual_;
      std::vector<double> error(n_);
      for (size_t i=0; i<n_; ++i) error[i] = rr[i]-r[i];
      double const rn = norm(r), en = norm(error);
      Row row;
      row.call = event.solverCall;
      row.pass = event.pass;
      row.restriction = event.restriction;
      row.values = {rn,en,rn>0 ? en/rn : en,std::sqrt(jumps[3]),
                    missing(),missing(),missing(),missing(),missing(),
                    missing(),missing(),missing(),missing(),missing(),
                    std::sqrt(jumps[0])/std::max(1.0,norm(x)),missing(),missing()};
      if (event.restriction)
      {
        write(row,false);
        return;
      }
      double const rhoTrue = dot(r,z);
      double const floor = 64*std::numeric_limits<double>::epsilon()*rn*norm(z);
      row.values[4] = event.rho;
      row.values[5] = rhoTrue;
      row.values[6] = std::abs(event.rho-rhoTrue);
      row.values[7] = std::abs(rhoTrue)>floor
                    ? ratio(row.values[6],std::abs(rhoTrue)) : missing();
      row.values[8] = !std::isfinite(rhoTrue) || std::abs(rhoTrue)<=floor;
      row.values[9] = event.beta;
      row.values[10] = event.curvature>0 ? event.rho/event.curvature : missing();
      row.values[11] = event.curvature;
      row.values[15] = std::sqrt(jumps[1])/std::max(1.0,norm(z));
      if (event.directionAvailable)
      {
        auto Aq = multiply(q);
        double const energy = dot(q,Aq);
        row.values[12] = energy;
        if (!previousQ_.empty() && energy>0 && previousEnergy_>0)
          row.values[13] = ratio(std::abs(dot(q,previousAq_)),
                                std::sqrt(energy)*std::sqrt(previousEnergy_));
        row.values[16] = std::sqrt(jumps[2])/std::max(1.0,norm(q));
        pendingQ_ = std::move(q);
        pendingAq_ = std::move(Aq);
        pendingEnergy_ = energy;
      }
      else
      {
        pendingQ_.clear();
        pendingAq_.clear();
      }
      pending_ = std::move(row);
      hasPending_ = true;
    }

    // Called after solve(), so rejected/stagnated directions do not become
    // the reference direction for the next conjugacy measurement.
    void finishCall(bool applied)
    {
      if (!root_ || !hasPending_) return;
      write(pending_,applied);
      if (applied && !pendingQ_.empty())
      {
        previousQ_ = pendingQ_;
        previousAq_ = pendingAq_;
        previousEnergy_ = pendingEnergy_;
      }
      hasPending_ = false;
    }

    bool good() const { return !root_ || output_.good(); }
    bool close()
    {
      if (!root_) return true;
      output_.flush();
      bool const ok = output_.good();
      output_.close();
      return ok && !output_.fail();
    }

  private:
    struct Row
    {
      bool restriction = false;
      int call = 0, pass = 0;
      std::vector<double> values;
    };

    std::vector<double> multiply(std::vector<double> const& x) const
    {
      Vector input(n_), output(n_);
      for (size_t i=0; i<n_; ++i) input[i][0] = x[i];
      matrix_.mv(input,output);
      std::vector<double> result(n_);
      for (size_t i=0; i<n_; ++i) result[i] = output[i][0];
      return result;
    }

    void write(Row const& row, bool applied)
    {
      output_ << (row.restriction ? "restriction" : "search") << ','
              << row.call << ',' << row.pass << ',';
      if (!row.restriction) output_ << applied;
      for (double value : row.values)
      {
        output_ << ',';
        if (std::isfinite(value)) output_ << value;
      }
      output_ << '\n';
      output_.flush(); // Keep the final diagnostic even when the solver aborts.
    }

    Matrix const& matrix_;
    Vector const& rhs_;
    Storage const& storage_;
    Data const& data_;
    bool root_;
    Sum sum_;
    size_t n_;
    int call_ = 0;
    int referenceCall_ = 0;
    std::ofstream output_;
    std::vector<double> previousLocalU_, previousQ_, previousAq_, pendingQ_, pendingAq_;
    std::vector<double> referenceX_, referenceResidual_;
    double previousEnergy_ = 0, pendingEnergy_ = 0;
    Row pending_;
    bool hasPending_ = false;
  };
}
#endif
