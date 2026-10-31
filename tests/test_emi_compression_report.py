#!/usr/bin/env python3
"""Check production report aggregation with simulated MPI ranks, without Dune."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def extract(source, anchor):
    start = source.index(anchor)
    opening = source.index('{', start)
    end, depth = opening+1, 1
    while depth:
        depth += (source[end]=='{')-(source[end]=='}')
        end += 1
    return source[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
                        default=Path(__file__).resolve().parents[1]/'EMI_bddc.hh')
    args = parser.parse_args()
    source = args.source.read_text()
    helper = extract(source,'  inline bool mpiRankZero(')
    report = extract(source,'  struct CompressionReport')+';'
    preamble = r'''
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#define KASKADE_HAVE_MPI
int simulatedRank = 0;
constexpr int MPI_COMM_WORLD=0, MPI_UNSIGNED_LONG_LONG=1, MPI_DOUBLE=2, MPI_SUM=3, MPI_MAX=4;
int MPI_Initialized(int* value) { *value=1; return 0; }
int MPI_Comm_rank(int,int* rank) { *rank=simulatedRank; return 0; }
// printGlobal is compiled but is not invoked by this local aggregation harness.
int MPI_Allreduce(void const* src,void* dst,int n,int type,int,int) {
  std::memcpy(dst,src,n*(type==MPI_DOUBLE ? sizeof(double) : sizeof(unsigned long long)));
  return 0;
}
namespace EmiBddc {
'''
    checks = r'''
} // namespace EmiBddc
int main() {
  using EmiBddc::CompressionReport;
  std::array<int,6> traffic{0,0,0,150,250,0};
  std::uint64_t globalRaw=0,globalWire=0,globalRawReceived=0,globalReceived=0;
  for (simulatedRank=0;simulatedRank<8;++simulatedRank) {
    CompressionReport report;
    auto add = [&]<class Report>(Report& r) {
      if constexpr (requires { r.addUncompressed(traffic,true); })
        r.addUncompressed(traffic,true);
      else
        r.addUncompressed(traffic);
    };
    add(report);
    globalRaw+=report.originalSentBytes; globalWire+=report.sentBytes;
    globalRawReceived+=report.originalReceivedBytes; globalReceived+=report.receivedBytes;
  }
  if(globalRaw!=400 || globalWire!=400 || globalRawReceived!=400 || globalReceived!=400) {
    std::cerr << "MPI-global uncompressed counters must be contributed once, not once per rank\n";
    return 1;
  }
  CompressionReport serial;
  serial.addUncompressed(traffic);
  serial.addUncompressed(traffic);
  if(serial.originalSentBytes!=800 || serial.sentBytes!=800) return 2;
  std::cout << "Compression report checks passed (8 MPI ranks, serial accumulation).\n";
}
'''
    with tempfile.TemporaryDirectory(prefix='emi-report-check-') as folder:
        cpp = Path(folder)/'check.cpp'
        binary = Path(folder)/'check'
        cpp.write_text(preamble+helper+'\n'+report+'\n'+checks)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2',str(cpp),'-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True)


if __name__ == '__main__':
    main()
