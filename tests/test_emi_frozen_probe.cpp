// Run serial, or compile with KASKADE_HAVE_MPI and run with two MPI ranks.
#include "EMI_preconditioner_probe.hh"
#include <sstream>
#include <cstdio>

int main(int argc,char** argv)
{
  if (argc!=2) return 2; // CSV scratch path supplied by caller.
  bool mpi=false;int rank=0,ranks=1;
#ifdef KASKADE_HAVE_MPI
  MPI_Init(&argc,&argv);mpi=true;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#endif
  using namespace EmiBddc;
  int calls=0,matvecs=0;
  auto a=[&](Fgmres::Vector x) {++matvecs;x[0]*=10;return x;};
  auto local=[&](Fgmres::Vector x) {++calls;for(auto& v:x)v*=-2.0/ranks;return x;};
  probeFirstPreconditioner(Fgmres::Vector{3,4},mpi,a,local,argv[1]);
  bool good=calls==1 && matvecs==(rank==0 ? 1 : 0);
  if (rank==0)
  {
    std::ifstream in(argv[1]);std::string line;
    std::getline(in,line);good=good && line=="dof,rhs,v,z,Az";
    for (int i=0;i<2;++i)
    {
      std::getline(in,line);std::replace(line.begin(),line.end(),',',' ');
      std::istringstream row(line);int id;double rhs,v,z,az;
      row>>id>>rhs>>v>>z>>az;
      double expected=i==0 ? 0.6 : 0.8;
      good=good && bool(row) && id==i && rhs==i+3 && std::abs(v-expected)<1e-15
           && std::abs(z+2*expected)<1e-14
           && std::abs(az+2*expected*(i==0 ? 10 : 1))<1e-13;
    }
    std::remove(argv[1]);
  }
  // Nonfinite output on one worker must fail collectively, not hang in Reduce.
  bool failed=false;
  try
  {
    probeFirstPreconditioner(Fgmres::Vector{3,4},mpi,a,
      [&](Fgmres::Vector x) {if(rank==ranks-1)x[0]=NAN;return x;},argv[1]);
  }
  catch (std::runtime_error const&) {failed=true;}
  good=good && failed;
  // A root-only filesystem failure must also propagate to every rank.
  failed=false;
  try {probeFirstPreconditioner(Fgmres::Vector{3,4},mpi,a,local,std::string(argv[1])+"/missing.csv");}
  catch (std::runtime_error const&) {failed=true;}
  good=good && failed;
  int status=good ? 0 : 1;
#ifdef KASKADE_HAVE_MPI
  MPI_Allreduce(MPI_IN_PLACE,&status,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  if(rank==0)std::cout<<(status ? "FAIL" : "PASS")
    <<": frozen first-input BDDC adapter, normalization, assembly, operator image, collective errors\n";
#ifdef KASKADE_HAVE_MPI
  MPI_Finalize();
#endif
  return status;
}
