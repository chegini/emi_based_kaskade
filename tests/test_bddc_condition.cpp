#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "../EMI_bddc_condition.hh"

namespace
{
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

  double dot(Vector const& a, Vector const& b)
  {
    double sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i) sum += a[i]*b[i];
    return sum;
  }

  // Independent textbook PCG on a small explicitly stored SPD problem.
  // All tested RHS entries excite the full spectrum unless specified otherwise.
  std::vector<std::pair<double,double>> pcg(Matrix const& a, Vector const& b,
                                           Vector const& inverseDiagonal)
  {
    Vector r = b, z(b.size()), direction(b.size()), ad(b.size());
    double previousRz = 0;
    std::vector<std::pair<double,double>> scalars;
    for (std::size_t k = 0; k < b.size(); ++k)
    {
      for (std::size_t i = 0; i < b.size(); ++i) z[i] = inverseDiagonal[i]*r[i];
      double const rz = dot(r,z);
      double const beta = k ? rz/previousRz : 0;
      for (std::size_t i = 0; i < b.size(); ++i)
        direction[i] = z[i]+beta*direction[i];
      for (std::size_t i = 0; i < b.size(); ++i)
      {
        ad[i] = 0;
        for (std::size_t j = 0; j < b.size(); ++j) ad[i] += a[i][j]*direction[j];
      }
      double const alpha = rz/dot(direction,ad);
      scalars.emplace_back(alpha,beta);
      for (std::size_t i = 0; i < b.size(); ++i) r[i] -= alpha*ad[i];
      if (dot(r,r) < 1e-26*dot(b,b)) break;
      previousRz = rz;
    }
    return scalars;
  }

  void close(double actual, double expected)
  {
    assert(std::isfinite(actual));
    assert(std::abs(actual/expected-1) < 1e-10);
  }

  void check(std::vector<std::pair<double,double>> const& scalars,
             double minimum, double maximum, double scale = 1)
  {
    EmiBddc::CgConditionEstimator estimator;
    for (auto [alpha,beta] : scalars) estimator.append(alpha/scale,beta);
    auto result = estimator.estimate();
    assert(std::string(result.status) == "estimated");
    assert(result.dimension == scalars.size());
    close(result.lambdaMin,minimum*scale);
    close(result.lambdaMax,maximum*scale);
    close(result.condition,maximum/minimum);
  }
}

int main()
{
  using namespace EmiBddc;
  Matrix const a = {{4,1},{1,3}};
  // Noncommuting diagonal preconditioner: tests the preconditioned spectrum,
  // not A's eigenvalues or the singular values of the nonsymmetric product BA.
  auto const two = pcg(a,{1,2},{0.25,0.5});
  assert(two.size() == 2);
  check(two,1.25-std::sqrt(0.1875),1.25+std::sqrt(0.1875));
  for (double scale : {1e-200,1e-20,1e20,1e200})
    check(two,1.25-std::sqrt(0.1875),1.25+std::sqrt(0.1875),scale);

  Matrix diagonal(4,Vector(4,0));
  diagonal[0][0] = 1; diagonal[1][1] = 3;
  diagonal[2][2] = 10; diagonal[3][3] = 30;
  auto const four = pcg(diagonal,{1,2,3,4},{2,1,0.5,0.25});
  assert(four.size() == 4);
  check(four,2,7.5);
  auto const otherRhs = pcg(diagonal,{4,3,2,1},{2,1,0.5,0.25});
  check(otherRhs,2,7.5); // a fresh estimator carries no previous RHS history

  CgConditionEstimator partial;
  for (std::size_t i = 0; i < 2; ++i)
    partial.append(four[i].first,four[i].second);
  auto const partialResult = partial.estimate();
  assert(std::string(partialResult.status) == "estimated");
  assert(partialResult.lambdaMin >= 2-1e-12 && partialResult.lambdaMax <= 7.5+1e-12);
  assert(partialResult.condition < 7.5/2); // unseen modes are not an exact condition number

  CgConditionEstimator none;
  assert(none.estimate().dimension == 0);
  assert(std::string(none.estimate().status) == "insufficient-iterations");
  auto const oneMode = pcg(diagonal,{1,0,0,0},{2,1,0.5,0.25});
  assert(oneMode.size() == 1);
  none.append(oneMode[0].first,oneMode[0].second);
  assert(std::string(none.estimate().status) == "insufficient-iterations");

  // A repeated eigenvalue must not divide by zero in the Sturm count.
  CgConditionEstimator repeated;
  repeated.append(0.5,0); repeated.append(0.5,0);
  auto const equal = repeated.estimate();
  assert(std::string(equal.status) == "estimated");
  close(equal.lambdaMin,2); close(equal.lambdaMax,2); close(equal.condition,1);

  double const nan = std::numeric_limits<double>::quiet_NaN();
  double const inf = std::numeric_limits<double>::infinity();
  for (auto [alpha,beta] : std::vector<std::pair<double,double>>{
      {0,0},{-1,0},{nan,0},{inf,0},{1,-1},{1,nan},{1,inf},{1,0.1}})
  {
    CgConditionEstimator invalid;
    invalid.append(alpha,beta);
    invalid.append(1,0); // an invalid history cannot silently become valid
    assert(std::string(invalid.estimate().status) == "invalid-cg-scalars");
  }
  std::cout << "BDDC CG/Lanczos condition estimate checks passed\n";
}
