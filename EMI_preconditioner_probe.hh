#ifndef EMI_PRECONDITIONER_PROBE_HH
#define EMI_PRECONDITIONER_PROBE_HH

#include "EMI_fgmres_mpi.hh"
#include <fstream>
#include <iomanip>
#include <string>

namespace EmiBddc
{
  // Diagnostic only: one independent application to the first Arnoldi input.
  // No Krylov update, PDE step or production convergence claim is made.
  template<class ApplyA, class LocalPreconditioner>
  void probeFirstPreconditioner(Fgmres::Vector rhs, bool mpi,
                               ApplyA applyA, LocalPreconditioner applyLocal,
                               std::string const& filename)
  {
    int rank=0;
#ifdef KASKADE_HAVE_MPI
    if (mpi) MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
    if (mpi) throw std::invalid_argument("Probe MPI requested in serial build");
#endif
    auto synchronizeError=[&](int error)
    {
#ifdef KASKADE_HAVE_MPI
      if (mpi) MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
      if (error) throw std::runtime_error("Frozen BDDC probe failed; inspect rank error output");
    };
    synchronizeError(rhs.empty() || rhs.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()));
    Fgmres::Vector v(rhs.size(),0),local,z(rhs.size(),0);
    double const rhsNorm=Fgmres::norm(rhs);
    synchronizeError(!Fgmres::finite(rhs) || !std::isfinite(rhsNorm) || rhsNorm==0);
    if (rank==0)
      for (std::size_t i=0;i<v.size();++i) v[i]=rhs[i]/rhsNorm;
#ifdef KASKADE_HAVE_MPI
    if (mpi) MPI_Bcast(v.data(),static_cast<int>(v.size()),MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    int error=0;
    try
    {
      local=applyLocal(v);
      if (local.size()!=v.size() || !Fgmres::finite(local)) error=1;
    }
    catch (std::exception const& e)
    {std::cerr << "Frozen probe local error: " << e.what() << '\n';error=1;}
    synchronizeError(error);
#ifdef KASKADE_HAVE_MPI
    if (mpi) MPI_Reduce(local.data(),z.data(),static_cast<int>(v.size()),MPI_DOUBLE,
                        MPI_SUM,0,MPI_COMM_WORLD);
    else
#endif
      z=std::move(local);
    if (rank==0)
    {
      try
      {
        auto az=applyA(z);
        if (az.size()!=v.size() || !Fgmres::finite(z) || !Fgmres::finite(az))
          throw std::runtime_error("Nonfinite/mismatched frozen correction or operator image");
        std::ofstream out(filename);
        out << "dof,rhs,v,z,Az\n" << std::setprecision(17);
        for (std::size_t i=0;i<v.size();++i)
          out << i << ',' << rhs[i] << ',' << v[i] << ',' << z[i] << ',' << az[i] << '\n';
        out.close();
        if (!out) throw std::runtime_error("Cannot write frozen probe CSV");
      }
      catch (std::exception const& e)
      {std::cerr << "Frozen probe root error: " << e.what() << '\n';error=1;}
    }
    synchronizeError(error);
  }
}
#endif
