#!/usr/bin/env python3
"""Compile production BDDC methods against a tiny variable-SPD-preconditioner problem.

No Dune/MPI/PDE dependencies. This tests the recurrence, not EMI convergence.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
from test_bddc_cg_observer import block

PREAMBLE = r'''
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
namespace Dune {
template<class T,int N> struct FieldVector:std::array<T,N> {
 FieldVector() { this->fill(0); }
};
}
namespace Kaskade::BDDC {
struct Vec:std::vector<double> {
 using std::vector<double>::vector;
 Vec& operator=(double a) { for(auto& v:*this) v=a; return *this; }
 Vec& operator+=(Vec const& x) { axpy(1,x);return *this; }
 void axpy(double a,Vec const& x) { for(size_t i=0;i<size();++i) (*this)[i]+=a*x[i]; }
};
struct Matrix {
 void mv(Vec const& x,Vec& y) const { y[0]=4*x[0]+x[1];y[1]=x[0]+3*x[1]; }
 void smv(double a,Vec const& x,Vec& y) const { mv(x,y);for(auto& v:y)v*=a; }
 void usmv(double a,Vec const& x,Vec& y) const { Vec z(2);mv(x,z);y.axpy(a,z); }
};
struct TransferMock {
 using QuantizedScalar=std::uint32_t;
 int exchanges=0;
 void postProlongate(Vec const&,Vec&) { ++exchanges; }
 void resetPublishedResidual() {}
};
template<int m,class Scalar_,class CoarseScalar,class Transfer> struct Subdomain {
 using Scalar=Scalar_; using Vector=Vec;
 int nx=2,cg_itr=0;
 bool cg_solver=true,cgUseOriginalResidual=true,flexibleCg=false;
 Matrix A; Transfer trans;
 Vec f{2.,-1.},r=f,u=Vec(2),rRes=f,du=Vec(2),r_bar=Vec(2);
 Vec q=Vec(2),q_pre=Vec(2),Aq=Vec(2),Adu=Vec(2),previousCgResidual=Vec(2);
 Vec rawCorrection=Vec(2),residual=Vec(2),restrictedResidual=Vec(2);
 double alpha=0;
 void imposeDCq(Vec&) {}
 void solveDirichlet(Vec&,Vec const&) {}
 Dune::FieldVector<double,2> prolongate();
 Dune::FieldVector<double,1> updateSearchDirectionNorm();
 void updateSearchDirection(double);
 void updateSolution(Scalar);
 void resetForNewRhs();
};
using S=Subdomain<1,double,double,TransferMock>;
struct Timings {
 static Timings& instance() { static Timings t;return t; }
 void start(char const*) {} void stop(char const*) {}
};
struct Domain {
 S s; bool variable=false,invalid=false;
 int restrictions=0,prolongations=0;
 int size() const { return 2; } // avoid the production one-subdomain direct-solve shortcut
 bool mpiRankZero() const { return true; }
 void updateIteration(int k) { s.cg_itr=k; }
 void restrict() { ++restrictions;s.rRes=s.r; }
 Dune::FieldVector<double,2> prolongate() {
   ++prolongations;
   // SPD M_k^{-1}: identity first, then a noncommuting SPD matrix.
   s.du=s.r;
   if(variable && s.cg_itr>0) {
     s.du[0]=3*s.r[0]+.7*s.r[1];s.du[1]=.7*s.r[0]+s.r[1];
   }
   auto p=s.prolongate();
   if(invalid && s.cg_itr>0) p[1]=std::numeric_limits<double>::quiet_NaN();
   return p;
 }
 void updateSearchDirection(double b) { s.updateSearchDirection(b); }
 auto updateSearchDirectionNorm() { return s.updateSearchDirectionNorm(); }
 void updateSolution(double a) { s.updateSolution(a); }
};
struct Coarse { void solve() {} };
template<class Subdomain,class Storage> struct BDDCSolver {
 enum class SolveStatus { updated,stagnated,breakdown };
 struct CgStep { double alpha=0,beta=0;bool applied=false; };
 Domain domain; Coarse coarseSolver;
 int systemIteration=0,restrictionPasses_=1;
 bool cg_solver=true,BDDC_verbose=false,flexibleCg_=false;
 double sigma_pre=0;
 SolveStatus lastStatus_=SolveStatus::updated; CgStep lastCgStep_;
 bool hasStagnated(double,double) { return false; }
 double solve();
'''

MAIN = r'''
double norm(Vec const& v) { return std::sqrt(v[0]*v[0]+v[1]*v[1]); }
double dot(Vec const& a,Vec const& b) { return a[0]*b[0]+a[1]*b[1]; }
void close(double a,double b) { if(std::abs(a-b)>1e-12)throw std::runtime_error("scalar mismatch"); }
void test() {
 using Solver=BDDCSolver<int,int>;
 // A fixed preconditioner reproduces PCG; a changing one requires a different beta.
 for(bool variable:{false,true}) {
   Solver cg,fcg;
   cg.domain.variable=fcg.domain.variable=variable;
   fcg.flexibleCg_=fcg.domain.s.flexibleCg=true;
   close(cg.solve(),fcg.solve());
   Vec first=fcg.domain.s.q;
   cg.solve();fcg.solve();
   Vec aq(2);fcg.domain.s.A.mv(fcg.domain.s.q,aq);
   close(dot(first,aq),0); // actual A-conjugacy, not just recurrence arithmetic
   if(norm(fcg.domain.s.r)>1e-12)throw std::runtime_error("Flexible CG failed changing SPD test");
   if(variable && norm(cg.domain.s.r)<1e-3)throw std::runtime_error("Test failed to distinguish standard CG");
   if(!variable) { close(cg.lastCgStep_.beta,fcg.lastCgStep_.beta);close(norm(cg.domain.s.r),0); }
   assert(fcg.domain.restrictions==2 && fcg.domain.prolongations==2);
   assert(cg.domain.restrictions==2 && cg.domain.prolongations==2);
 }
 // Production local dot product and residual history: unquantized r is used.
 S s;s.flexibleCg=true;s.cg_itr=1;s.r={1.,2.};s.rRes={100.,-200.};
 s.previousCgResidual={3.,4.};s.du={.5,.25};
 auto p=s.prolongate();close(p[0],1.);close(p[1],-1.5);
 assert(s.previousCgResidual==s.r);
 s.resetForNewRhs();close(norm(s.previousCgResidual),0);
 assert(s.flexibleCg && s.cgUseOriginalResidual && s.cg_itr==0);
 // Negative PR beta is valid and must not be silently replaced by FR or clamped.
 Solver negative;negative.flexibleCg_=negative.domain.s.flexibleCg=true;
 negative.systemIteration=1;negative.sigma_pre=2;
 negative.domain.s.r={1.,0.};negative.domain.s.previousCgResidual={2.,0.};
 negative.solve();close(negative.lastCgStep_.beta,-.5);
 // A nonfinite PR numerator rejects before the solution update.
 Solver bad;bad.flexibleCg_=bad.domain.s.flexibleCg=true;
 bad.solve();auto u=bad.domain.s.u;bad.domain.invalid=true;bad.solve();
 assert(bad.lastStatus_==Solver::SolveStatus::breakdown && bad.domain.s.u==u && !bad.lastCgStep_.applied);
}
}
int main() {
 try { Kaskade::BDDC::test();std::cout<<"Production flexible-CG checks passed: variable SPD, fixed-PCG equivalence, conjugacy, reset, negative beta, rejection, one pass.\n"; }
 catch(std::exception const& e) { std::cerr<<e.what()<<'\n';return 1; }
}
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kaskade-root',type=Path,default=Path('/kaskade7'))
    args=p.parse_args()
    hh=(args.kaskade_root/'mg/bddc.hh').read_text()
    hpp=(args.kaskade_root/'mg/bddc.hpp').read_text()
    event=block(hh,'    struct DiagnosticEvent')+';'
    methods=[]
    prefix='  template <int m, class Scalar_, class CoarseScalar, class Transfer>\n'
    for signature in ['Dune::FieldVector<double,2> Subdomain<m,Scalar_,CoarseScalar,Transfer>::prolongate()',
                      'Dune::FieldVector<double,1> Subdomain<m,Scalar_,CoarseScalar,Transfer>::updateSearchDirectionNorm()',
                      'void Subdomain<m,Scalar_,CoarseScalar,Transfer>::updateSearchDirection(double beta)',
                      'void Subdomain<m,Scalar_,CoarseScalar,Transfer>::updateSolution(Scalar alpha_)',
                      'void Subdomain<m,Scalar_,CoarseScalar,Transfer>::resetForNewRhs()']:
        methods.append(block(hpp,prefix+'  '+signature))
    methods.append(block(hpp,'  template <class Subdomain, class Storage>\n  double BDDCSolver<Subdomain,Storage>::solve()'))
    cpp=PREAMBLE+event+'\nstd::function<void(DiagnosticEvent const&)> diagnosticObserver_;\n};\n'+ '\n'.join(methods)+MAIN
    with tempfile.TemporaryDirectory(prefix='bddc-flexible-') as folder:
        src=Path(folder)/'test.cpp';exe=Path(folder)/'test';src.write_text(cpp)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++20','-O2',str(src),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)

if __name__=='__main__':main()
