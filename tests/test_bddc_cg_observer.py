#!/usr/bin/env python3
"""Compile the production solve() against a small orchestration harness.

No PDE runs: compare arithmetic/update traces with diagnostics disabled/enabled,
including rejection paths. The container run remains the MPI integration check.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def block(source, anchor):
    start = source.index(anchor)
    opening = source.index('{', start)
    depth, end = 1, opening+1
    while depth:
        depth += (source[end] == '{')-(source[end] == '}')
        end += 1
    return source[start:end]


PREAMBLE = r'''
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
namespace Dune { template<class T,int n> using FieldVector=std::array<T,n>; }
namespace Kaskade::BDDC {
struct Timings {
 static Timings& instance() { static Timings t; return t; }
 void start(char const*) {} void stop(char const*) {}
};
struct MockDomain {
 std::vector<double> trace;
 double rho=2,curvature=5;
 int size() const { return 2; }
 bool mpiRankZero() const { return true; }
 void updateIteration(int n) { trace.push_back(100+n); }
 void restrict() { trace.push_back(200); }
 Dune::FieldVector<double,2> prolongate() { trace.push_back(300); return {rho,4}; }
 void updateSearchDirection(double beta) { trace.push_back(beta); }
 Dune::FieldVector<double,1> updateSearchDirectionNorm() { return {curvature}; }
 void updateSolution(double alpha) { trace.push_back(alpha); }
};
struct Coarse { int count=0; void solve() { ++count; } };
template<class Subdomain,class Storage> struct BDDCSolver {
 enum class SolveStatus { updated,stagnated,breakdown };
 struct CgStep { double alpha=0,beta=0; bool applied=false; };
 MockDomain domain;
 Coarse coarseSolver;
 int systemIteration=0,restrictionPasses_=2;
 bool cg_solver=true,BDDC_verbose=false,stagnate=false;
 bool flexibleCg_=false;
 double sigma_pre=0;
 SolveStatus lastStatus_=SolveStatus::updated;
 CgStep lastCgStep_;
 bool hasStagnated(double,double) { return stagnate; }
 double solve();
'''
CHECKS = r'''
};
}
'''
MAIN = r'''
int main() {
 using Solver=Kaskade::BDDC::BDDCSolver<int,int>;
 for (int scenario=0;scenario<5;++scenario) {
   Solver a,b;
   std::vector<Solver::DiagnosticEvent> events;
   b.diagnosticObserver_=[&](auto const& event){events.push_back(event);};
   double firstA=a.solve(),firstB=b.solve();
   assert(firstA==firstB);
   assert(events.size()==3 && events[0].restriction && events[0].pass==1);
   assert(events[1].restriction && events[1].pass==2);
   assert(!events[2].restriction && events[2].rho==2 && events[2].curvature==4);
   assert(a.lastCgStep_.applied && b.lastCgStep_.applied);
   if(scenario==1) a.domain.rho=b.domain.rho=-1;
   if(scenario==2) a.domain.curvature=b.domain.curvature=-1;
   if(scenario==3) a.stagnate=b.stagnate=true;
   if(scenario==4) a.cg_solver=b.cg_solver=false;
   double secondA=a.solve(),secondB=b.solve();
   assert(secondA==secondB);
   assert(a.domain.trace==b.domain.trace && a.coarseSolver.count==b.coarseSolver.count);
   assert(a.sigma_pre==b.sigma_pre && a.systemIteration==b.systemIteration);
   assert(a.lastStatus_==b.lastStatus_);
   assert(a.lastCgStep_.applied==b.lastCgStep_.applied);
   assert(a.lastCgStep_.alpha==b.lastCgStep_.alpha && a.lastCgStep_.beta==b.lastCgStep_.beta);
   assert(events.size()==6);
   if(scenario==1) assert(!events.back().directionAvailable);
 }
 std::cout<<"Production CG observer checks passed (arithmetic unchanged, all passes and rejection paths).\n";
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kaskade-root', type=Path, default=Path('/kaskade7'))
    args = parser.parse_args()
    hh = (args.kaskade_root/'mg/bddc.hh').read_text()
    hpp = (args.kaskade_root/'mg/bddc.hpp').read_text()
    event = block(hh, '    struct DiagnosticEvent')+';'
    method = block(hpp, '  template <class Subdomain, class Storage>\n'
                        '  double BDDCSolver<Subdomain,Storage>::solve()')
    with tempfile.TemporaryDirectory(prefix='emi-cg-observer-') as d:
        cpp, executable = Path(d)/'check.cpp', Path(d)/'check'
        cpp.write_text(PREAMBLE+event+'\nstd::function<void(DiagnosticEvent const&)> diagnosticObserver_;\n'
                       +CHECKS+'namespace Kaskade::BDDC {\n'+method+'\n}\n'+MAIN)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2',str(cpp),'-o',str(executable)],check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == '__main__':
    main()
