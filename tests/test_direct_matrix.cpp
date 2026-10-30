#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "../EMI_direct_matrix.hh"
#include "emi_material_test_matrix.hh"

int main()
{
  using namespace EmiBddc;
  MaterialTestMatrix matrix(3);
  matrix[0][0][0][0] = 4;
  matrix[1][1][0][0] = 3;
  matrix[2][2][0][0] = 2;
  matrix[0][1][0][0] = matrix[1][0][0][0] = -1;
  matrix[1][2][0][0] = matrix[2][1][0][0] = 1e-300;
  auto const checked = checkSpdDirectInput(matrix);
  assert(checked.lowerNonzeros == 5 && checked.maxAbsSkew == 0);
  std::vector<std::tuple<int,int,double>> entries;
  visitDirectLowerNonzeros(matrix,[&](int i, int j, double value)
  {
    assert(i >= j && value != 0);
    entries.emplace_back(i,j,value);
  });
  assert(entries.size() == 5);
  assert((entries[4] == std::tuple<int,int,double>(2,2,2)));
  MaterialTestMatrix reconstructed(3);
  for (auto [i,j,value] : entries)
  {
    reconstructed[i][j][0][0] = value;
    reconstructed[j][i][0][0] = value;
  }
  // Mirroring the emitted half reconstructs exactly the full operator;
  // nonzero entries, however tiny, are retained without thresholding.
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      assert(reconstructed[i][j][0][0] == matrix[i][j][0][0]);

  auto rejects = [](auto const& input)
  {
    try { (void)checkSpdDirectInput(input); }
    catch (std::runtime_error const&) { return true; }
    return false;
  };
  assert(rejects(MaterialTestMatrix()));
  auto bad = matrix;
  bad[0][1][0][0] = -2;
  assert(rejects(bad));
  bad = matrix; bad.erase(1,0);
  assert(rejects(bad)); // missing transposed nonzero is not silently discarded
  bad = matrix; bad.erase(0,0);
  assert(rejects(bad));
  for (double value : {0.0,-1.0,std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()})
  {
    bad = matrix; bad[1][1][0][0] = value;
    assert(rejects(bad));
  }
  bad = matrix; bad[0][2][0][0] = std::numeric_limits<double>::quiet_NaN();
  assert(rejects(bad));
  auto sparse = matrix;
  sparse.erase(0,2); sparse.erase(2,0);
  assert(checkSpdDirectInput(sparse).lowerNonzeros == 5);

  // Symmetry and positive diagonals are necessary but not an SPD proof.
  // The real backend must reject a failed positive-definite factorization.
  MaterialTestMatrix indefinite(2);
  indefinite[0][0][0][0] = indefinite[1][1][0][0] = 1;
  indefinite[0][1][0][0] = indefinite[1][0][0][0] = 2;
  assert(checkSpdDirectInput(indefinite).maxAbsSkew == 0);
  std::cout << "SPD direct symmetry and exact lower-triangle checks passed\n";
}
