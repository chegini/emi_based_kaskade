#!/usr/bin/env python3
"""Run the actual Kaskade R/P exchange methods with small deterministic MPI payloads.

Only the subdomain/solver algebra is replaced; production routing, headers, dispatch,
traffic accumulation and receive order are compiled verbatim. No Dune dependency.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def method(source, name):
    anchor = '  void SharedMemoryDomain<Subdomain,Storage>::'+name+'()'
    start = source.index(anchor)
    opening = source.index('{', start)
    end, depth = opening+1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return 'template <class Subdomain, class Storage>\n'+source[start:end]


def routing_setter(source):
    start = source.index('    void setMpiDestinationRouting(bool enabled)')
    opening = source.index('{', start)
    end, depth = opening+1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[start:end]


PREAMBLE = r'''
#include "mg/bddcMpi.hh"
#include "mg/bddcMpiRouting.hh"
#include <cstring>
#include <utility>
namespace Kaskade::BDDC {
struct SmallSubdomain {
  using TransmissionType=std::vector<double>;
  int id;
  std::vector<int> peers;
  std::vector<std::pair<int,TransmissionType>> restricted,prolongated;
  std::vector<int> const& neighbors() const { return peers; }
  void getInterfaceResidualUpdate(int to,TransmissionType& out) const {
    out={double(id),double(to),1e16,1.,-1e16};
    out.resize(5+(id+to)%4,double(id-to));
  }
  void getInterfaceValues(int to,TransmissionType& out) const {
    getInterfaceResidualUpdate(to,out); out.push_back(0.125);
  }
  void setInterfaceResidualUpdate(int from,TransmissionType const& in) {
    restricted.emplace_back(from,in);
  }
  void setInterfaceValues(int from,TransmissionType const& in) {
    prolongated.emplace_back(from,in);
  }
};
template<class Subdomain,class Storage>
struct SharedMemoryDomain {
  using Subdom=Subdomain;
  Storage* subdomains;
  std::vector<bool> activeSubdomains;
  MpiDomainLayout mpiLayout;
  MpiNodeTopology nodeTopology;
  MpiRoutingTraffic routedTraffic;
  bool destinationRouting=false;
  int restrictionTraffic=0,prolongationTraffic=0,transfers=0;
  SharedMemoryDomain(Storage& s,bool routed,bool master):subdomains(&s),
      activeSubdomains(s.size(),true) {
    mpiLayout.configure(true,s.size(),master);
    setMpiDestinationRouting(routed);
  }
  bool mpiEnabled() const { return mpiLayout.enabled(); }
  bool mpiRankZero() const { return mpiLayout.rank()==0; }
  int size() const { return subdomains->size(); }
  int mpiSize() const { return mpiLayout.size(); }
  int mpiOwner(int id) const { return mpiLayout.owner(id); }
  bool mpiOwns(int id) const { return mpiLayout.owns(id); }
  bool activeOnThisRank(int id) const { return activeSubdomains[id] && mpiOwns(id); }
  void mpiRestrictExchange();
  void mpiProlongateExchange();
  // PRODUCTION_ROUTING_SETTER
};
'''

CHECKS = r'''
} // namespace
using namespace Kaskade::BDDC;
void require(bool ok,char const* message) { if(!ok) throw std::runtime_error(message); }
void check(bool master,int ranks) {
  int n=3*ranks+1;
  std::vector<SmallSubdomain> a(n),b(n);
  for(int from=0;from<n;++from) {
    a[from].id=from;
    for(int to=0;to<n;++to)
      if(from!=to && (from+to)%3!=0) a[from].peers.push_back(to);
  }
  b=a;
  SharedMemoryDomain<SmallSubdomain,decltype(a)> legacy(a,false,master),routed(b,true,master);
  for(int pass=0;pass<4;++pass) {
    for(int i=0;i<n;++i) {
      bool active=pass!=3 && (pass!=2 || i%2==0);
      legacy.activeSubdomains[i]=routed.activeSubdomains[i]=active;
    }
    legacy.mpiRestrictExchange(); routed.mpiRestrictExchange();
    legacy.mpiProlongateExchange(); routed.mpiProlongateExchange();
    for(int i=0;i<n;++i)
      require(a[i].restricted==b[i].restricted && a[i].prolongated==b[i].prolongated,
              "Production R/P payload or accumulation order changed");
    require(legacy.transfers==routed.transfers &&
            legacy.restrictionTraffic==routed.restrictionTraffic &&
            legacy.prolongationTraffic==routed.prolongationTraffic,"Logical traffic changed");
  }
  std::array<std::uint64_t,6> total{};
  MPI_Allreduce(routed.routedTraffic.bytes.data(),total.data(),6,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
  require(total[0]+total[1]+total[2]==std::uint64_t(legacy.restrictionTraffic),"R accounting");
  require(total[3]+total[4]+total[5]==std::uint64_t(legacy.prolongationTraffic),"P accounting");
  require(total[0]>0 && total[3]>0,"Multiple local subdomains not exercised");
  if(routed.nodeTopology.nodes()==1) require(total[2]==0 && total[5]==0,"False inter-node bytes");
  if(ranks-(master?1:0)>1) require(total[1]+total[2]>0,"No cross-rank messages exercised");
  bool rejected=false;
  try { routed.setMpiDestinationRouting(false); }
  catch(std::logic_error const&) { rejected=true; }
  require(rejected,"Changing routing mid-solve must be rejected");
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank,ranks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    check(false,ranks);
    if(ranks>1) check(true,ranks);
    if(rank==0) std::cout<<"PASS: production R/P exchanges; exact payload/order/traffic parity; "
                            "multiple subdomains, inactive sets and idle master; ranks="<<ranks<<"\n";
    MPI_Finalize(); return 0;
  } catch(std::exception const& e) {
    std::cerr<<e.what()<<"\n"; MPI_Abort(MPI_COMM_WORLD,1); return 1;
  }
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kaskade', type=Path, default=Path(os.environ.get('KASKADE7', '/kaskade7')))
    parser.add_argument('--ranks', type=int, nargs='+', default=[1, 2, 4])
    a = parser.parse_args()
    if any(n < 1 for n in a.ranks):
        parser.error('MPI ranks must be positive')
    source = (a.kaskade/'mg/bddc.hpp').read_text()
    with tempfile.TemporaryDirectory(prefix='bddc-routing-') as folder:
        cpp, binary = Path(folder)/'check.cpp', Path(folder)/'check'
        preamble = PREAMBLE.replace('// PRODUCTION_ROUTING_SETTER',
                                   routing_setter((a.kaskade/'mg/bddc.hh').read_text()))
        cpp.write_text(preamble + method(source, 'mpiRestrictExchange') + '\n' +
                       method(source, 'mpiProlongateExchange') + CHECKS)
        subprocess.run(shlex.split(os.environ.get('MPICXX', 'mpicxx')) +
                       ['-std=c++20', '-O2', '-Wall', '-Wextra', '-Werror', '-DKASKADE_HAVE_MPI',
                        '-I', str(a.kaskade), str(cpp), '-o', str(binary)], check=True)
        for ranks in a.ranks:
            subprocess.run(shlex.split(os.environ.get('MPIEXEC', 'mpirun')) +
                           ['-n', str(ranks), str(binary)], check=True, timeout=120)


if __name__ == '__main__':
    main()
