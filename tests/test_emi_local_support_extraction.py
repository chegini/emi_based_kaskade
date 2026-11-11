#!/usr/bin/env python3
"""Exercise the production local-matrix extraction function without Dune/MPI."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

from test_emi_compression_report import extract


def main():
    repo = Path(__file__).resolve().parents[1]
    source = (repo/'EMI_bddc.hh').read_text()
    data = 'template<class Matrix,class Vector>\n'+extract(source, 'struct BddcData')+';'
    function = 'template<class Matrix,class Vector>\n'+extract(
        source, 'void setLocalMatricesFromContributions(')
    preamble = r"""
#include <algorithm>
#include <iostream>
#include "EMI_local_support_diagnostic.hh"
#include "tests/emi_material_test_matrix.hh"
namespace Kaskade::BDDC { struct LocalDof { int s, i; }; }
namespace EmiBddc {
bool mpiRankZero(bool) { return true; }
"""
    checks = r"""
}
struct Matrix : MaterialTestMatrix
{
  using MaterialTestMatrix::MaterialTestMatrix;
  Matrix() = default;
  Matrix(std::vector<size_t> const& indices, Matrix const& source)
    : MaterialTestMatrix(indices.size())
  {
    for (size_t i=0; i<indices.size(); ++i)
    {
      auto row = source[indices[i]];
      for (auto col=row.begin(); col!=row.end(); ++col)
      {
        auto found = std::find(indices.begin(),indices.end(),col.index());
        if (found!=indices.end()) (*this)[i][found-indices.begin()] = *col;
      }
    }
  }
};
struct Vector : std::vector<std::array<double,1>>
{
  using std::vector<std::array<double,1>>::vector;
  Vector& operator=(double value) { for(auto& v:*this) v[0]=value; return *this; }
};
void require(bool ok,char const* msg) { if(!ok) throw std::runtime_error(msg); }
int main()
{
  EmiBddc::BddcData<Matrix,Vector> data;
  data.tags={255,256};
  data.localDofs={{0,1,2},{2,3}};
  data.dofSubdomains.resize(4);
  std::vector<Matrix> contributions{Matrix(4),Matrix(4)};
  auto& m=contributions[0];
  for(size_t i=0;i<3;++i) m[i][i][0][0]=0.015854964663192563;
  m[0][1][0][0]=m[1][0][0][0]=1e-200;
  for(size_t i=0;i<3;++i)
    m[3][i][0][0]=m[i][3][0][0]=8.1211812265281262e-21;
  contributions[1][2][2][0][0]=0.2;
  contributions[1][3][3][0][0]=0.3;
  std::ostringstream log;
  auto* previous=std::cout.rdbuf(log.rdbuf());
  EmiBddc::setLocalMatricesFromContributions(data,contributions,false,false);
  std::cout.rdbuf(previous);
  require(data.localDofs==std::vector<std::vector<size_t>>{{0,1,2},{2,3}},
          "Support changed");
  require(data.localMatrices[0].N()==3 && data.localMatrices[1].N()==2,
          "Local matrix dimensions changed");
  require(data.localMatrices[0][0][1][0][0]==1e-200,"Tiny supported entry was dropped");
  require(data.localMatrices[0][1][0][0][0]==1e-200,"Symmetry changed");
  require(contributions[0][3][0][0][0]==8.1211812265281262e-21,
          "Global-coordinate contribution mutated");
  require(data.sharedDofs.size()==1 && data.sharedDofs[0].size()==2,
          "Shared interface metadata changed");
  require(data.weights[0][2][0]==0.5 && data.weights[1][0][0]==0.5,
          "Partition weights changed");
  require(log.str().find("omitted_entries=6")!=std::string::npos &&
          log.str().find("global_operator=unchanged")!=std::string::npos,
          "Roundoff omission was not reported in quiet mode");
  m[3][0][0][0]=1e-6;
  bool rejected=false;
  try { EmiBddc::setLocalMatricesFromContributions(data,contributions,false,false); }
  catch(std::runtime_error const& e)
  {
    rejected=std::string(e.what()).find("roundoff_budget=")!=std::string::npos;
  }
  require(rejected,"Production extraction accepted significant missing coupling");
  // With no out-of-space entries, extraction must preserve all retained bits
  // and must not report any omitted roundoff.
  for(size_t i=0;i<3;++i) m[3][i][0][0]=m[i][3][0][0]=0;
  log.str("");
  previous=std::cout.rdbuf(log.rdbuf());
  EmiBddc::setLocalMatricesFromContributions(data,contributions,false,false);
  std::cout.rdbuf(previous);
  require(log.str().empty(),"Exact support unexpectedly reported roundoff");
  require(data.localMatrices[0][0][1][0][0]==1e-200,"Exact support changed entries");
  std::cout << "Production local extraction passed: bounded omission, unchanged retained entries/metadata, significant-defect rejection\n";
}
"""
    with tempfile.TemporaryDirectory(prefix='emi-local-extraction-') as tmp:
        cpp, binary = Path(tmp)/'check.cpp', Path(tmp)/'check'
        cpp.write_text(preamble+data+'\n'+function+'\n'+checks)
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++'))+
                       ['-std=c++20', '-O2', '-Wall', '-Wextra', '-Werror',
                        '-I', str(repo), str(cpp), '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    main()
