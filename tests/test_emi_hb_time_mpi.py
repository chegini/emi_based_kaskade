#!/usr/bin/env python3
"""Compile real R/P routing with complete packet fields; no Dune/PDE mock timing."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import test_bddc_mpi_exchange as base


def build(kaskade,exe):
    header=(kaskade/'mg/bddc.hh').read_text()
    fields=header[header.index('    std::uint64_t packetOriginalBytes='):header.index('    MpiNodeTopology nodeTopology;',header.index('    std::uint64_t packetOriginalBytes='))]
    preamble=base.PREAMBLE.replace('#include <cstring>','#include <cstring>\n#include "mg/bddcPacketPacking.hh"')
    preamble=preamble.replace('  int restrictionTraffic=0,prolongationTraffic=0,transfers=0;',
        '  int restrictionTraffic=0,prolongationTraffic=0,transfers=0;\n  bool losslessPacking=false;\n'+fields)
    preamble=preamble.replace('  bool mpiEnabled() const',
        '  bool storageHas(int) const { return true; }\n  int mpiRank() const { return mpiLayout.rank(); }\n  bool mpiEnabled() const')
    preamble=preamble.replace('// PRODUCTION_ROUTING_SETTER',base.routing_setter(header))
    checks=base.CHECKS.replace('  for(int pass=0;pass<4;++pass) {','''  auto before=TransferTiming::snapshot();
  for(int pass=0;pass<4;++pass) {''')
    checks=checks.replace('  std::array<std::uint64_t,6> total{};', '''  auto delta=TransferTiming::difference(TransferTiming::snapshot(),before);
  if(TransferTiming::mode()) {
    std::uint64_t calls=0;
    for(auto const& phase:delta)
      calls+=phase[static_cast<unsigned>(TransferTiming::Stage::PacketExchange)][0].calls;
    require(calls==8,"Missing or extra destination-exchange timer calls");
  }
  if(routed.nodeTopology.nodes()==1)
    require(routed.packetInterNodeSentBytes==0,"Same-node packets counted as inter-node");
  std::array<std::uint64_t,6> total{};''')
    source=(kaskade/'mg/bddc.hpp').read_text()
    with tempfile.TemporaryDirectory() as tmp:
        cpp=Path(tmp)/'exchange.cpp'
        cpp.write_text(preamble+base.method(source,'mpiRestrictExchange')+'\n'+base.method(source,'mpiProlongateExchange')+checks)
        subprocess.run(shlex.split(os.environ.get('MPICXX','mpicxx'))+[
            '-std=c++20','-O2','-Wall','-Wextra','-Werror','-DKASKADE_HAVE_MPI',
            '-I',str(kaskade),str(cpp),'-o',str(exe)],check=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kaskade',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--build-only',action='store_true');p.add_argument('--ranks',type=int,default=4)
    a=p.parse_args();build(a.kaskade,a.executable)
    if not a.build_only:
        for mode in ['0','2']:
            subprocess.run(shlex.split(os.environ.get('MPIEXEC','mpiexec'))+['-n',str(a.ranks),str(a.executable)],
                env=dict(os.environ,KASKADE_TRANSFER_TIMING=mode,KASKADE_TIMESTEP_BYTES='1'),check=True,timeout=120)
