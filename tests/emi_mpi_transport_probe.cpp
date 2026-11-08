// Transport integrity check only: no mesh, solver, compression, or benchmark claim.
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

unsigned char symbol(int sender, int receiver, int index)
{
  return static_cast<unsigned char>((sender*31+receiver*17+index*13+index/257)%251);
}

void requireAll(bool ok, MPI_Comm comm, char const* operation)
{
  int local=ok, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,comm);
  if (!global) {
    std::cerr << "FAIL: " << operation << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

void broadcasts(MPI_Comm comm, char const* label)
{
  int rank=0,size=0;
  MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
  for (int count: {1,65536,3046272,16777216}) {
    if (!rank) std::cout << label << " broadcasts bytes=" << count << std::endl;
    std::vector<unsigned char> data(count);
    for (int root=0;root<size;++root) {
      std::fill(data.begin(),data.end(),255);
      if (rank==root) for (int i=0;i<count;++i) data[i]=symbol(root,0,i);
      MPI_Bcast(data.data(),count,MPI_UNSIGNED_CHAR,root,comm);
      bool ok=true;
      for (int i=0;i<count;++i) if (data[i]!=symbol(root,0,i)) {ok=false;break;}
      requireAll(ok,comm,"broadcast data mismatch");
    }
  }
}

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);
  int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  requireAll(argc==3 && size==std::atoi(argv[1]),MPI_COMM_WORLD,"unexpected rank count/arguments");
  MPI_Comm local;
  MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
  int leader=rank;
  MPI_Allreduce(MPI_IN_PLACE,&leader,1,MPI_INT,MPI_MIN,local);
  std::vector<int> leaders(size);
  MPI_Allgather(&leader,1,MPI_INT,leaders.data(),1,MPI_INT,MPI_COMM_WORLD);
  int nodes=std::set<int>(leaders.begin(),leaders.end()).size();
  requireAll(nodes==std::atoi(argv[2]),MPI_COMM_WORLD,"unexpected shared-memory node count");
  if (!rank) {
    char library[MPI_MAX_LIBRARY_VERSION_STRING]; int length=0;
    MPI_Get_library_version(library,&length);
    std::cout.write(library,length); std::cout << std::endl;
  }
  // No TCP BTL is enabled in the first, single-node stage of the launcher.
  broadcasts(local,"shared-node");
  MPI_Barrier(MPI_COMM_WORLD);
  broadcasts(MPI_COMM_WORLD,"world");
  for (int count: {65536,3046272,16777216}) {
    if (!rank) std::cout << "world ring bytes=" << count << std::endl;
    int next=(rank+1)%size, previous=(rank+size-1)%size;
    std::vector<unsigned char> send(count),recv(count);
    for (int i=0;i<count;++i) send[i]=symbol(rank,next,i);
    MPI_Sendrecv(send.data(),count,MPI_UNSIGNED_CHAR,next,7,
                 recv.data(),count,MPI_UNSIGNED_CHAR,previous,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    bool ok=true;
    for (int i=0;i<count;++i) if (recv[i]!=symbol(previous,rank,i)) {ok=false;break;}
    requireAll(ok,MPI_COMM_WORLD,"ring data mismatch");
  }
  if (!rank) std::cout << "world variable-size destination exchange" << std::endl;
  std::vector<int> sc(size),rc(size),sd(size),rd(size);
  for (int peer=0;peer<size;++peer) sc[peer]=(rank==peer ? 0 : 1048576+rank*13+peer*7);
  MPI_Alltoall(sc.data(),1,MPI_INT,rc.data(),1,MPI_INT,MPI_COMM_WORLD);
  int sn=0,rn=0;
  for (int peer=0;peer<size;++peer) {sd[peer]=sn;rd[peer]=rn;sn+=sc[peer];rn+=rc[peer];}
  std::vector<unsigned char> send(sn),recv(rn);
  for (int peer=0;peer<size;++peer)
    for (int i=0;i<sc[peer];++i) send[sd[peer]+i]=symbol(rank,peer,i);
  MPI_Alltoallv(send.data(),sc.data(),sd.data(),MPI_UNSIGNED_CHAR,
                recv.data(),rc.data(),rd.data(),MPI_UNSIGNED_CHAR,MPI_COMM_WORLD);
  bool ok=true;
  for (int peer=0;peer<size;++peer)
    for (int i=0;i<rc[peer];++i) if (recv[rd[peer]+i]!=symbol(peer,rank,i)) {ok=false;break;}
  requireAll(ok,MPI_COMM_WORLD,"destination data mismatch");
  if (!rank) std::cout << "PASS: transport payload integrity; ranks=" << size << " nodes=" << nodes << std::endl;
  MPI_Comm_free(&local);
  MPI_Finalize();
}
