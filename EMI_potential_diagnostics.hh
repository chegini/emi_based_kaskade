#ifndef EMI_POTENTIAL_DIAGNOSTICS_HH
#define EMI_POTENTIAL_DIAGNOSTICS_HH

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

namespace EmiBddc
{
  // Read-only, full-precision global algebraic DOFs, not interpolated VTK points.
  // The initial row has zero correction/RHS; these are not an initial solve.
  template<class Matrix, class State, class Vector>
  bool writePotentialState(std::string const& directory, int step, Matrix const& A,
                           State const& state, Vector const& du, Vector const& rhs)
  {
    auto const path = std::filesystem::path(directory) /
                      ("potential-state-step" + std::to_string(step) + ".csv");
    if (std::filesystem::exists(path)) return false;
    std::ofstream out(path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "dof,u,correction,rhs,residual\n";
    for (auto row=A.begin(); row!=A.end(); ++row)
    {
      auto const i=row.index();
      long double residual=rhs[i][0];
      for (auto col=row.begin(); col!=row.end(); ++col)
        residual-=static_cast<long double>((*col)[0][0])*du[col.index()][0];
      out << i << ',' << state[i][0] << ',' << du[i][0] << ','
          << rhs[i][0] << ',' << residual << '\n';
    }
    out.close();
    if (!out) return false;
    if (step==0)
    {
      auto const opPath=std::filesystem::path(directory)/"potential-operator.csv";
      if (std::filesystem::exists(opPath)) return false;
      std::vector<long double> rows(A.N(),0), columns(A.N(),0), diagonal(A.N(),0);
      for (auto row=A.begin(); row!=A.end(); ++row)
        for (auto col=row.begin(); col!=row.end(); ++col)
        {
          auto const value=(*col)[0][0];
          rows[row.index()]+=value;
          columns[col.index()]+=value;
          if (row.index()==col.index()) diagonal[row.index()]=value;
        }
      std::ofstream op(opPath);
      op << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "dof,a_one,a_transpose_one,diagonal\n";
      for (size_t i=0;i<rows.size();++i)
        op << i << ',' << rows[i] << ',' << columns[i] << ',' << diagonal[i] << '\n';
      op.close();
      if (!op) return false;
    }
    return true;
  }

  template<class Matrix, class State, class Vector>
  void capturePotentialState(bool enabled, bool rankZero, bool mpiEnabled,
                             std::string const& directory, int step, Matrix const& A,
                             State const& state, Vector const& du, Vector const& rhs)
  {
    if (!enabled) return;
    int ready=1;
    if (rankZero)
    {
      try { ready=writePotentialState(directory,step,A,state,du,rhs) ? 1 : 0; }
      catch (std::exception const&) { ready=0; }
    }
#ifdef KASKADE_HAVE_MPI
    if (mpiEnabled) MPI_Bcast(&ready,1,MPI_INT,0,MPI_COMM_WORLD);
#else
    (void)mpiEnabled;
#endif
    if (!ready) throw std::runtime_error("Cannot write fresh potential diagnostics in "+directory);
  }
}
#endif
