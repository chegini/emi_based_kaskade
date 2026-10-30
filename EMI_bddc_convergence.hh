#ifndef EMI_BDDC_CONVERGENCE_HH
#define EMI_BDDC_CONVERGENCE_HH

#include <fstream>
#include <iomanip>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>

#include "EMI_bddc_stopping.hh"

namespace EmiBddc
{
  // Rank-zero observer of existing true-residual checks. It neither applies
  // corrections nor decides convergence. Iteration zero is the accepted
  // initial correction after interior elimination; alpha is absent there.
  class ConvergenceReport
  {
  public:
    bool open(std::string const& path, int timeStep, double target, std::ostream& output)
    {
      timeStep_ = timeStep;
      target_ = target;
      output_ = &output;
      file_.open(path);
      if (!file_)
        return false;
      file_ << std::setprecision(17)
            << "time_step,iteration,step_length,residual_norm,rhs_norm,relative_residual,target\n";
      output << "BDDC true-residual convergence step " << timeStep << ":\n"
             << "  iteration 0 = initial accepted correction after interior elimination; alpha is absent.\n"
             << "  Other rows follow the update using alpha; rel is absolute when RHS is zero.\n"
             << "  iteration    step length   ||b-A*du||_2       ||b||_2            rel\n";
      return static_cast<bool>(file_);
    }

    void record(int iteration, std::optional<double> alpha, TrueResidual const& residual)
    {
      file_ << timeStep_ << ',' << iteration << ',';
      if (alpha)
        file_ << *alpha;
      file_ << ',' << residual.absolute << ',' << residual.rhsNorm << ','
            << residual.relative() << ',' << target_ << '\n';
      // Retain partial histories when a solve fails; I/O errors are checked
      // collectively by the caller before advancing the timestep state.
      file_.flush();
      std::ostringstream row;
      row << "  " << std::setw(9) << iteration << std::scientific << std::setprecision(6);
      if (alpha)
        row << std::setw(15) << *alpha;
      else
        row << std::setw(15) << "--";
      row << std::setw(15) << residual.absolute << std::setw(15) << residual.rhsNorm
          << std::setw(15) << residual.relative() << '\n';
      *output_ << row.str() << std::flush;
    }

    bool close()
    {
      file_.flush();
      bool const written = static_cast<bool>(file_);
      file_.close();
      return written && !file_.fail();
    }

  private:
    std::ofstream file_;
    std::ostream* output_ = nullptr;
    int timeStep_ = 0;
    double target_ = 0.0;
  };
}

#endif
