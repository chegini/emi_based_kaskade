#!/usr/bin/env python3
"""Compile the production prolongate method without Dune; MPI/PDE checks are separate."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

PREAMBLE = r'''
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace Dune {
template<class T,int N> struct FieldVector : std::array<T,N> {
  FieldVector() { this->fill(0); }
};
}
namespace Kaskade::BDDC {
struct Vec : std::vector<double> {
  using std::vector<double>::vector;
  Vec& operator=(double value) { for(auto& x:*this) x=value; return *this; }
  Vec& operator+=(Vec const& x) {
    for(size_t i=0;i<size();++i) (*this)[i]+=x[i]; return *this;
  }
};
struct Matrix {
  void mv(Vec const& x,Vec& y) const { y[0]=4*x[0]+x[1]; y[1]=x[0]+3*x[1]; }
  void smv(double a,Vec const& x,Vec& y) const { mv(x,y); for(auto& v:y) v*=a; }
};
struct TransferMock {
  int exchanges=0;
  void postProlongate(Vec const&,Vec&) { ++exchanges; }
};
template<int m,class Scalar_,class CoarseScalar,class Transfer>
struct Subdomain {
  using Vector=Vec;
  int nx=2,cg_itr=0;
  bool cg_solver=true,cgUseOriginalResidual=false,constrained=false;
  bool flexibleCg=false;
  Vec previousCgResidual=Vec(2);
  Transfer trans;
  Matrix A;
  Vec r{2.,-1.},rRes{2.125,-1.25},du{.2,-.1},r_bar=Vec(2);
  Vec q{.4,.3},Aq=Vec(2),Adu=Vec(2);
  void imposeDCq(Vec& x) { if(constrained) x[0]=0; }
  void solveDirichlet(Vec&,Vec const&) {}
  Dune::FieldVector<double,2> prolongate();
};
'''
CHECKS = r'''
void close(double a,double b) {
  if(!std::isfinite(a) || std::abs(a-b)>2e-14) throw std::runtime_error("wrong scalar");
}
using S=Subdomain<1,double,double,TransferMock>;
void test() {
  for(int iteration : {0,1}) for(bool constrained : {false,true}) {
    S legacy; legacy.cg_itr=iteration; legacy.constrained=constrained;
    auto original=legacy; original.cgUseOriginalResidual=true;
    auto r=original.r;
    auto p=legacy.prolongate(),t=original.prolongate();
    close(p[0],constrained ? .125 : .55);
    close(t[0],constrained ? .1 : .5);
    close(p[1],t[1]);
    if(original.du!=legacy.du || original.r_bar!=legacy.r_bar ||
       original.q!=legacy.q || original.Aq!=legacy.Aq ||
       original.r!=r || original.rRes!=legacy.rRes ||
       original.trans.exchanges!=1 || legacy.trans.exchanges!=1)
      throw std::runtime_error("dot product changed vectors or transfer count");
    // When restriction is exact, both choices coincide.
    S a; a.cg_itr=iteration; a.rRes=a.r; S b=a; b.cgUseOriginalResidual=true;
    close(a.prolongate()[0],b.prolongate()[0]);
    // Steepest descent deliberately retains its existing scalar.
    S sd; sd.cg_solver=false; S sdOn=sd; sdOn.cgUseOriginalResidual=true;
    auto s=sd.prolongate(),sOn=sdOn.prolongate();
    close(s[0],sOn[0]); close(s[1],sOn[1]);
  }
  // Two local copies of a shared DOF: original residual contributions sum,
  // while the continuous correction is the same on both copies.
  S a,b; a.r={2.,0.}; b.r={3.,0.}; a.du={.25,0.}; b.du=a.du;
  a.rRes={100.,0.}; b.rRes={-80.,0.};
  a.cgUseOriginalResidual=b.cgUseOriginalResidual=true;
  close(a.prolongate()[0]+b.prolongate()[0],(2.+3.)*.25);
}
}
int main() {
  try { Kaskade::BDDC::test(); std::cout<<"Original-residual CG rho checks passed.\n"; }
  catch(std::exception const& e) { std::cerr<<e.what()<<"\n"; return 1; }
}
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kaskade-root',type=Path,default=Path('/kaskade7'))
    args=p.parse_args()
    header=(args.kaskade_root/'mg/bddc.hh').read_text()
    assert 'bool cgUseOriginalResidual = false;' in header
    assert 'void setCgUseOriginalResidual(bool enabled) { cgUseOriginalResidual = enabled; }' in header
    source=(args.kaskade_root/'mg/bddc.hpp').read_text()
    anchor=('  template <int m, class Scalar_, class CoarseScalar, class Transfer>\n'
            '  Dune::FieldVector<double,2> Subdomain<m,Scalar_,CoarseScalar,Transfer>::prolongate()')
    assert source.count(anchor)==1
    start=source.index(anchor); end=source.index('{',start)+1; depth=1
    while depth:
        depth+=(source[end]=='{')-(source[end]=='}'); end+=1
    with tempfile.TemporaryDirectory(prefix='bddc-rho-') as folder:
        cpp=Path(folder)/'check.cpp'; exe=Path(folder)/'check'
        cpp.write_text(PREAMBLE+source[start:end]+CHECKS)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2',str(cpp),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=='__main__': main()
