#!/usr/bin/env python3
"""Compile production applyPreconditioner against an ownership/state harness."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
from test_bddc_cg_observer import block

PREAMBLE = r'''
#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
struct V:std::vector<double> {
 using std::vector<double>::vector;
 V& operator+=(V const& b) {for(size_t i=0;i<size();++i) (*this)[i]+=b[i];return *this;}
};
template<class F> void parallelFor(int a,int b,F f) {for(int i=a;i<b;++i) f(i);}
struct Sub {
 V u,du; int resets=0;
 void getSolution(V& x) {x=u;}
 void getCorrection(V& x) {x=du;}
};
struct Domain {
 std::vector<Sub>* s; int restrictions=0,prolongations=0,iteration=-1;
 bool distributed=false;
 int size() {return s->size();}
 bool storageHas(int) {return true;}
 bool mpiEnabled() {return distributed;}
 bool mpiOwns(int i) {return i==0;}
 void updateIteration(int n) {iteration=n;}
 void restrict() {++restrictions;}
 void prolongate() {++prolongations;} // A negative rho is deliberately not inspected.
};
struct Coarse {int calls=0;void solve() {++calls;}};
struct Solver {
 using XVector=V;
 std::vector<Sub>* subdomains;Domain domain;Coarse coarseSolver;int restrictionPasses_=1;
 int restartCount=0;
 Solver(std::vector<Sub>& s):subdomains(&s),domain{&s} {}
 void restartRhs(std::vector<V> const& rhs) {
   ++restartCount;
   for(size_t i=0;i<rhs.size();++i) if(!domain.mpiEnabled() || domain.mpiOwns(i)) {
     ++(*subdomains)[i].resets;
     (*subdomains)[i].u=rhs[i]; (*subdomains)[i].du=V(rhs[i].size(),-2.0);
   }
 }
'''
MAIN = r'''
};
int main() {
 std::vector<Sub> s(2);Solver solver(s);std::vector<V> rhs{V{3,4},V{5,6}},out;
 solver.applyPreconditioner(rhs,out);
 assert(out[0][0]==1 && out[0][1]==2 && out[1][0]==3 && out[1][1]==4);
 assert(s[0].u[0]==3); // No CG solution update. Initial interior part is included in output.
 assert(solver.domain.restrictions==1 && solver.domain.prolongations==1 && solver.coarseSolver.calls==1);
 rhs[0][0]=9; solver.applyPreconditioner(rhs,out);
 assert(out[0][0]==7 && s[0].resets==2 && solver.restartCount==2);
 solver.domain.distributed=true;out.clear();solver.applyPreconditioner(rhs,out);
 assert(out[0].size()==2 && out[1].empty() && s[1].resets==2);
 std::vector<Sub> one(1);Solver single(one);single.applyPreconditioner({V{7}},out);
 assert(out[0][0]==7 && single.domain.restrictions==0 && single.coarseSolver.calls==0);
 std::cout<<"BDDC preconditioner application checks passed (interior contribution, no CG update, fresh RHS, ownership, one subdomain).\n";
}
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kaskade-root',type=Path,default=Path('/kaskade7'))
    args=p.parse_args()
    source=(args.kaskade_root/'mg/bddc.hh').read_text()
    method=block(source,'    void applyPreconditioner(')
    with tempfile.TemporaryDirectory() as temp:
        cpp=Path(temp)/'check.cpp';exe=Path(temp)/'check'
        cpp.write_text(PREAMBLE+method+MAIN)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2',str(cpp),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)

if __name__=='__main__':main()
