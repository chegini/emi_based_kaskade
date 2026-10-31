#!/usr/bin/env python3
"""Compile the production updateSolution method against a tiny matrix/vector harness.

This isolates residual bookkeeping from Dune and factorization dependencies. The
container EMI runs remain the integration check for actual subdomains and MPI.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

PREAMBLE = r'''
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace Kaskade::BDDC {
struct Vec : std::vector<double> {
  using std::vector<double>::vector;
  void axpy(double a, Vec const& x) {
    for (std::size_t i=0;i<size();++i) (*this)[i] += a*x[i];
  }
};
struct Matrix {
  std::array<std::array<double,2>,2> a{{{{4,1}},{{1,3}}}};
  void mv(Vec const& x, Vec& y) const {
    for(int i=0;i<2;++i) { y[i]=0; for(int j=0;j<2;++j) y[i]+=a[i][j]*x[j]; }
  }
  void usmv(double scale, Vec const& x, Vec& y) const {
    Vec product(2); mv(x,product); y.axpy(scale,product);
  }
};
struct ExactTransfer {};
struct QuantizedTransfer { using QuantizedScalar = std::uint32_t; };
template<int m, class Scalar_, class CoarseScalar, class Transfer>
struct Subdomain {
  using Scalar = Scalar_;
  using Vector = Vec;
  Matrix A;
  Vec f{2.0,-1.0}, u{0.3,-0.2}, r=Vec(2), rRes=Vec(2);
  Vec q{0.2,-0.1}, q_pre=Vec(2), Aq=Vec(2), du{0.2,-0.1}, Adu=Vec(2);
  bool cg_solver = true;
  double alpha = 0;
  void imposeDCq(Vector&) {} // All DOFs are active in this regression.
  void updateSolution(Scalar alpha_);
};
'''
CHECKS = r'''
Vec trueResidual(auto const& s) {
  Vec result=s.f; s.A.usmv(-1,s.u,result); return result;
}
void check(Vec const& actual, Vec const& expected, char const* label) {
  for(int i=0;i<2;++i)
    if (!std::isfinite(actual[i]) || std::abs(actual[i]-expected[i])>2e-14)
      throw std::runtime_error(label);
}
void testExactCg() {
  Subdomain<1,double,double,ExactTransfer> s;
  for(int step=0;step<2;++step) {
    s.rRes=trueResidual(s); s.updateSolution(0.4); check(s.r,trueResidual(s),"Exact CG residual changed");
  }
}
void testCompressedCg() {
  Subdomain<1,double,double,QuantizedTransfer> s;
  for(int step=0;step<2;++step) {
    s.rRes=trueResidual(s);
    // Deliberately inject transfer error into the averaged residual. It must
    // not become an error in the stored residual of the original local system.
    s.rRes[0]+=0.125; s.rRes[1]-=0.25;
    s.updateSolution(0.4);
    check(s.r,trueResidual(s),"Compressed CG fed quantized residual error into the original RHS");
    check(s.q_pre,s.q,"CG search history changed");
  }
}
void testSteepestDescent() {
  Subdomain<1,double,double,QuantizedTransfer> s; s.cg_solver=false;
  s.r=trueResidual(s); s.A.mv(s.du,s.Adu); s.updateSolution(0.4);
  check(s.r,trueResidual(s),"Steepest descent residual changed");
}
} // namespace Kaskade::BDDC
int main() {
  try {
    using namespace Kaskade::BDDC;
    testExactCg(); testSteepestDescent(); testCompressedCg();
    std::cout << "BDDC residual update checks passed (exact CG, compressed CG, steepest descent).\n";
  } catch(std::exception const& error) {
    std::cerr << error.what() << "\n"; return 1;
  }
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kaskade-root', type=Path, default=Path('/kaskade7'))
    args = parser.parse_args()
    source = (args.kaskade_root/'mg/bddc.hpp').read_text()
    anchor = ('  template <int m, class Scalar_, class CoarseScalar, class Transfer>\n'
              '  void Subdomain<m,Scalar_,CoarseScalar,Transfer>::updateSolution(Scalar alpha_)')
    if source.count(anchor) != 1:
        raise SystemExit('Production method anchor changed; review this regression harness.')
    start = source.index(anchor)
    opening = source.index('{',start)
    depth = 1
    end = opening+1
    while depth:
        depth += (source[end]=='{')-(source[end]=='}')
        end += 1
    method = source[start:end]
    with tempfile.TemporaryDirectory(prefix='bddc-residual-check-') as folder:
        cpp = Path(folder)/'check.cpp'
        executable = Path(folder)/'check'
        cpp.write_text(PREAMBLE+method+'\n'+CHECKS)
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2',str(cpp),'-o',str(executable)],check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == '__main__':
    main()
