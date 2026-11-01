// Tiny MPI adapter test; no mesh, PDE or Kaskade library required.
#include "EMI_fgmres_mpi.hh"
#include <iostream>

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);
  int rank=0,ranks=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  using namespace EmiBddc;
  auto identity=[](Fgmres::Vector const& x) {return x;};
  int calls=0,observations=0;
  auto local=[&](Fgmres::Vector x) {
    ++calls; for (auto& a:x) a=-a/ranks;return x;
  };
  OuterSolverTraffic traffic;
  auto result=rootFgmres(Fgmres::Vector{1,2,3},true,30,10,1e-12,identity,local,
                        [&](int,int,double,double) {++observations;},traffic);
  // Per application: command int, vector bcast/reduce, error int Allreduce.
  // End: command int, seven-double info, solution vector. No matvec/dot collectives.
  std::uint64_t expected=(ranks-1)*(4ULL+24+24+8+4+56+24);
  bool good=result.status==Fgmres::Status::converged && result.iterations==1 && calls==1
         && traffic.bytes==expected && traffic.collectives==7
         && observations==(rank==0 ? 2 : 0);
  for (int i=0;i<3;++i) good=good && std::abs(result.x[i]-(i+1))<1e-12;
  OuterSolverTraffic badTraffic;
  auto bad=rootFgmres(Fgmres::Vector{1,2,3},true,30,10,1e-7,identity,
    [&](Fgmres::Vector x) {if(rank==ranks-1) x[0]=NAN;return x;},
    [](int,int,double,double) {},badTraffic);
  good=good && bad.status==Fgmres::Status::nonfinite;
  int status=good ? 0 : 1,global=0;
  MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0) std::cout<<(global ? "FAIL" : "PASS")
    <<": root FGMRES MPI adapter, negative preconditioner, replicated result, failure synchronization, payload counts\n";
  MPI_Finalize();return global;
}
