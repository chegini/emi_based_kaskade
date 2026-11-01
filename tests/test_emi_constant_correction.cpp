#include "../EMI_constant_correction.hh"
#include "../EMI_bddc_stopping.hh"
#include "emi_material_test_matrix.hh"
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
using Vec=std::vector<std::array<double,1>>;
double residualNorm(MaterialTestMatrix const& A,Vec const& b,Vec const& x)
{
  double norm=0;
  for(auto row=A.begin();row!=A.end();++row) {
    double r=b[row.index()][0];
    for(auto col=row.begin();col!=row.end();++col) r-=(*col)[0][0]*x[col.index()][0];
    norm+=r*r;
  }
  return std::sqrt(norm);
}
int main()
{
  MaterialTestMatrix A(2);
  A[0][0][0][0]=3; A[0][1][0][0]=-1;
  A[1][0][0][0]=-1; A[1][1][0][0]=2;
  Vec b{{4},{2}},x{{.2},{-.4}},before=x;
  auto p=EmiBddc::proposeConstantCorrection(A,b,x);
  assert(p.valid && std::abs(p.delta-2)<1e-14);
  assert(x==before); // proposal is read-only, not a mean removal
  for(auto& v:x) v[0]+=p.delta;
  assert(std::abs((x[0][0]-x[1][0])-(before[0][0]-before[1][0]))<1e-14);
  auto remaining=EmiBddc::proposeConstantCorrection(A,b,x);
  assert(std::abs(remaining.residualSum)<1e-14);
  // An exact constant error in a weakly grounded Laplacian is removed.
  A[0][0][0][0]=1.00001; A[1][1][0][0]=1;
  b={{{.00001}},{{0}}};x={{{1.125}},{{1.125}}};
  p=EmiBddc::proposeConstantCorrection(A,b,x);
  assert(p.valid && std::abs(p.delta+.125)<1e-10);
  // A true null mode must not be regularized silently.
  A[0][0][0][0]=1;
  assert(!EmiBddc::proposeConstantCorrection(A,b,x).valid);
  A[0][0][0][0]=std::numeric_limits<double>::quiet_NaN();
  assert(!EmiBddc::proposeConstantCorrection(A,b,x).valid);
  // Energy projection can increase the Euclidean residual: the production
  // path must recheck its actual stopping criterion before advancing state.
  A=0; A[0][0][0][0]=100; A[1][1][0][0]=1;
  b={{{1e-8}},{{-2e-8}}};x={{{0}},{{0}}};
  auto initial=residualNorm(A,b,x);
  p=EmiBddc::proposeConstantCorrection(A,b,x);
  for(auto& v:x) v[0]+=p.delta;
  auto final=residualNorm(A,b,x);
  assert(final>initial);
  double tolerance=(initial+final)/2;
  assert((EmiBddc::TrueResidual{initial,1}).meets(tolerance));
  assert(!(EmiBddc::TrueResidual{final,1}).meets(tolerance));
  std::cout<<"Constant-mode correction tests passed.\n";
}
