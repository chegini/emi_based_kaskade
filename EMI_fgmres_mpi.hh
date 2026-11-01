#ifndef EMI_FGMRES_MPI_HH
#define EMI_FGMRES_MPI_HH

#include <cstdint>
#include <iostream>
#include <limits>

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif
#include "EMI_fgmres.hh"

namespace EmiBddc
{
  /** Logical tree-equivalent payload, not hardware network traffic.
   * Bcast/Reduce: (P-1)*payload; Allreduce: 2*(P-1)*payload.
   * Each collective is counted once globally. This ledger covers outer solver
   * communication only, separate from BDDC's interface/coarse/codebook traffic.
   */
  struct OuterSolverTraffic
  {
    std::uint64_t bytes=0, collectives=0;
    void add(std::size_t payload, int ranks, bool allreduce=false)
    {
      ++collectives;
      bytes+=static_cast<std::uint64_t>(payload)*(ranks-1)*(allreduce ? 2 : 1);
    }
  };

  /** Root-owned Arnoldi and exact operator; BDDC is collective/distributed.
   * Workers receive each normalized Arnoldi vector and reduce their weighted
   * local preconditioned contributions. No CG positivity test or update occurs.
   * Matvec/orthogonalization on root add no MPI messages; their CPU/memory cost
   * remains real. This is a validation adapter, not a scalable Krylov layout.
   */
  template<class ApplyA, class LocalPreconditioner, class Observe>
  Fgmres::Result rootFgmres(Fgmres::Vector const& rhs, bool mpi, int restart,
                            int limit, double tolerance, ApplyA applyA,
                            LocalPreconditioner applyLocal, Observe observe,
                            OuterSolverTraffic& traffic)
  {
    if (rhs.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::invalid_argument("FGMRES vector exceeds MPI count limit");
    int rank=0, ranks=1;
#ifdef KASKADE_HAVE_MPI
    if (mpi) {MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);}
#else
    if (mpi) throw std::invalid_argument("FGMRES MPI requested in serial build");
#endif
    auto command=[&](int tag)
    {
#ifdef KASKADE_HAVE_MPI
      if (mpi)
      {
        MPI_Bcast(&tag,1,MPI_INT,0,MPI_COMM_WORLD);
        traffic.add(sizeof(int),ranks);
      }
#endif
      return tag;
    };
    auto apply=[&](Fgmres::Vector v)
    {
#ifdef KASKADE_HAVE_MPI
      if (mpi)
      {
        MPI_Bcast(v.data(),static_cast<int>(v.size()),MPI_DOUBLE,0,MPI_COMM_WORLD);
        traffic.add(v.size()*sizeof(double),ranks);
      }
#endif
      Fgmres::Vector local;
      int error=0;
      try
      {
        local=applyLocal(v);
        if (local.size()!=v.size() || !Fgmres::finite(local)) error=1;
      }
      catch (std::exception const& e)
      {
        std::cerr << "FGMRES local preconditioner error: " << e.what() << '\n'; error=1;
      }
#ifdef KASKADE_HAVE_MPI
      if (mpi)
      {
        MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        traffic.add(sizeof(int),ranks,true);
      }
#endif
      if (error) return Fgmres::Vector(v.size(),std::numeric_limits<double>::quiet_NaN());
      Fgmres::Vector global(v.size(),0);
#ifdef KASKADE_HAVE_MPI
      if (mpi)
      {
        MPI_Reduce(local.data(),global.data(),static_cast<int>(v.size()),MPI_DOUBLE,
                   MPI_SUM,0,MPI_COMM_WORLD);
        traffic.add(v.size()*sizeof(double),ranks);
      }
      else
#endif
        global=std::move(local);
      return global;
    };
    Fgmres::Result result;
    if (rank==0)
    {
      try
      {
        result=Fgmres::solve(rhs,restart,limit,tolerance,applyA,
          [&](Fgmres::Vector const& v) {command(1);return apply(v);},observe);
      }
      catch (std::exception const& e)
      {
        std::cerr << "FGMRES root error: " << e.what() << '\n';
        result.status=Fgmres::Status::nonfinite;
        result.x.assign(rhs.size(),0); result.rhsNorm=Fgmres::norm(rhs);
        result.residual=std::numeric_limits<double>::infinity();
      }
      command(0);
    }
    else
    {
      while (command(0)==1) apply(Fgmres::Vector(rhs.size(),0));
    }
#ifdef KASKADE_HAVE_MPI
    if (mpi)
    {
      double info[7]={static_cast<double>(result.status),double(result.iterations),
        double(result.cycles),double(result.matvecs),result.residual,result.rhsNorm,
        double(result.peakBasisBytes)};
      MPI_Bcast(info,7,MPI_DOUBLE,0,MPI_COMM_WORLD);
      traffic.add(sizeof(info),ranks);
      result.status=static_cast<Fgmres::Status>(static_cast<int>(info[0]));
      result.iterations=static_cast<int>(info[1]); result.cycles=static_cast<int>(info[2]);
      result.matvecs=static_cast<int>(info[3]); result.residual=info[4];result.rhsNorm=info[5];
      result.peakBasisBytes=static_cast<std::size_t>(info[6]);
      result.x.resize(rhs.size());
      MPI_Bcast(result.x.data(),static_cast<int>(rhs.size()),MPI_DOUBLE,0,MPI_COMM_WORLD);
      traffic.add(rhs.size()*sizeof(double),ranks);
    }
#else
    (void)traffic; (void)ranks;
#endif
    return result;
  }
}
#endif
