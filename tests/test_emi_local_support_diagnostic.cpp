#include "../EMI_local_support_diagnostic.hh"
#include "emi_material_test_matrix.hh"
#include <iostream>

void require(bool ok, char const* message)
{
  if (!ok) throw std::runtime_error(message);
}

template <class Action>
void rejects(Action action)
{
  bool failed = false;
  try { action(); }
  catch (std::runtime_error const&) { failed = true; }
  require(failed,"Significant or invalid support defect was accepted");
}

int main()
{
  MaterialTestMatrix matrix(4);
  matrix[0][0][0][0] = 8.0;
  matrix[0][2][0][0] = 1e-200; // Tiny missing entries must still be reported.
  matrix[3][1][0][0] = -0.125;
  auto text = EmiBddc::describeMissingMaterialSupport(matrix,{0,1},2,108);
  for (auto token : {"subdomain_index=2", "material_tag=108",
                    "missing_entries=2", "missing_dofs=2",
                    "largest_missing_row=3", "largest_missing_col=1",
                    "largest_missing_value=-0.125", "row_present=0",
                    "col_present=1", "contribution_max_abs=8",
                    "missing_global_dofs_first16=2,3"})
    if (text.find(token) == std::string::npos)
      throw std::runtime_error(std::string("Missing diagnostic: ")+token);
  if (matrix[0][2][0][0] != 1e-200 || matrix[3][1][0][0] != -0.125)
    throw std::runtime_error("Diagnostic changed the matrix");
  auto covered = EmiBddc::describeMissingMaterialSupport(matrix,{0,1,2,3},2,108);
  if (covered.find("missing_entries=0 missing_dofs=0") == std::string::npos)
    throw std::runtime_error("Covered matrix reported missing entries");
  rejects([&] { EmiBddc::checkMaterialLocalSupport(matrix,{0,1},2,108); });

  // Kermit-shaped defect: six off-support entries, with the reported maximum.
  // Check scale independence and preserve tiny entries INSIDE the local space.
  for (double scale : {1e-150,1.0,1e150})
  {
    MaterialTestMatrix roundoff(4);
    roundoff[0][0][0][0] = 0.015854964663192563*scale;
    roundoff[0][1][0][0] = 1e-150*scale;
    for (size_t i = 0; i != 3; ++i)
    {
      roundoff[3][i][0][0] = 8.1211812265281262e-21*scale;
      roundoff[i][3][0][0] = 8.1211812265281262e-21*scale;
    }
    auto check = EmiBddc::checkMaterialLocalSupport(roundoff,{0,1,2},115,255);
    require(check.omittedEntries == 6,"Wrong omission count");
    require(check.omittedAbsSumBound >= 6.0L*roundoff[3][0][0][0],
            "Reported perturbation bound is too small");
    require(check.omittedAbsSumBound <= check.roundoffBudget(),"Accepted excessive perturbation");
    require(roundoff[0][1][0][0] == 1e-150*scale,"Tiny supported entry changed");
    require(roundoff[3][0][0][0] == 8.1211812265281262e-21*scale,
            "Source contribution changed");
    roundoff[3][0][0][0] = 1e-8*scale;
    rejects([&] { EmiBddc::checkMaterialLocalSupport(roundoff,{0,1,2},115,255); });
  }

  // Each entry is less than epsilon, but their absolute sum exceeds it.
  // Opposite signs must not cancel in the bound.
  MaterialTestMatrix aggregate(4);
  aggregate[0][0][0][0] = 1;
  for (size_t i = 0; i != 3; ++i)
  {
    aggregate[i][3][0][0] = std::numeric_limits<double>::epsilon()/4;
    aggregate[3][i][0][0] = -std::numeric_limits<double>::epsilon()/4;
  }
  rejects([&] { EmiBddc::checkMaterialLocalSupport(aggregate,{0,1,2},0,100); });

  MaterialTestMatrix tiny(2);
  tiny[1][1][0][0] = 1e-200;
  rejects([&] { EmiBddc::checkMaterialLocalSupport(tiny,{0},0,100); });
  auto exact = EmiBddc::checkMaterialLocalSupport(tiny,{0,1},0,100);
  require(exact.omittedEntries == 0 && exact.retainedMaxAbs == 1e-200,
          "Small but supported matrix was not preserved");
  for (double bad : {std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()})
  {
    tiny[1][1][0][0] = bad;
    rejects([&] { EmiBddc::checkMaterialLocalSupport(tiny,{0},0,100); });
    rejects([&] { EmiBddc::checkMaterialLocalSupport(tiny,{0,1},0,100); });
  }
  std::cout << "Local support checks passed: bounded roundoff, scaling, retained values, aggregate and significant-defect rejection\n";
}
