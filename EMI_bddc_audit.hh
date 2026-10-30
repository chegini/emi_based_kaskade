#ifndef EMI_BDDC_AUDIT_HH
#define EMI_BDDC_AUDIT_HH

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "EMI_bddc.hh"
#include "EMI_bddc_components.hh"

namespace EmiBddc
{
  // All coefficients are one for a constant field in the nodal Lagrange space
  // used by this application. The ratios are diagnostics, not stopping rules.
  template <class Matrix, class Vector>
  void auditConstantMode(Matrix const& matrix, Vector const& rhs, bool mpiEnabled)
  {
    if (!mpiRankZero(mpiEnabled))
      return;
    if (matrix.N() == 0)
      return;
    Vector ones(matrix.N()), image(matrix.N());
    ones = 1.0;
    image = 0.0;
    matrix.umv(ones,image);
    double matrixInf = 0.0, imageInf = 0.0;
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
    {
      double rowSum = 0.0;
      for (auto col = row->begin(); col != row->end(); ++col)
        rowSum += std::abs((*col)[0][0]);
      matrixInf = std::max(matrixInf,rowSum);
      imageInf = std::max(imageInf,std::abs(image[row.index()][0]));
    }
    double const n = static_cast<double>(matrix.N());
    double const rhsNorm = std::sqrt(std::max(0.0,rhs.dot(rhs)));
    double const rhsSum = ones.dot(rhs);
    std::ostringstream report;
    report << std::setprecision(12)
           << "BDDC constant-mode audit (current LHS and initial RHS):\n"
           << "  ||A*1||_2: " << std::sqrt(std::max(0.0,image.dot(image))) << '\n'
           << "  ||A*1||_inf: " << imageInf << '\n'
           << "  ||A||_inf: " << matrixInf << '\n'
           << "  ||A*1||_inf / ||A||_inf: "
           << (matrixInf > 0.0 ? imageInf/matrixInf : 0.0) << '\n'
           << "  constant Rayleigh quotient (1^T*A*1)/(1^T*1): "
           << ones.dot(image)/n << '\n'
           << "  initial RHS constant component 1^T*b: " << rhsSum << '\n'
           << "  normalized initial RHS constant component: "
           << (rhsNorm > 0.0 ? std::abs(rhsSum)/(std::sqrt(n)*rhsNorm) : 0.0)
           << '\n';
    std::cout << report.str();
  }

  namespace AuditDetail
  {
    template <class Range>
    std::string joined(Range const& values)
    {
      std::ostringstream result;
      bool first = true;
      for (auto value : values)
      {
        if (!first)
          result << ';';
        result << value;
        first = false;
      }
      return result.str();
    }
  }

  template <class Grid, class Space, class Material>
  std::vector<int> electricalCellLabels(Grid const& grid, Space const& space,
                                        Material const& material,
                                        std::vector<int> const& extraTags)
  {
    // Match the application's electrical continuity labels: all extra tags
    // map to zero, while each intracellular material keeps its own tag.
    std::set<int> const extracellular(extraTags.begin(),extraTags.end());
    std::vector<int> cellLabels(grid.size(0));
    Dune::FieldVector<double,Grid::dimension> zero(0.0);
    for (auto const& cell : Dune::elements(grid.leafGridView()))
    {
      size_t const index = space.indexSet().index(cell);
      int const tag = material.value(cell,zero);
      cellLabels.at(index) = extracellular.count(tag) ? 0 : tag;
    }
    return cellLabels;
  }

  // Change only the coarse partition; the assembled operator, physical BCs,
  // overlap, and diagonal weights remain unchanged. Each component has one
  // electrical label because the original-cell graph does not join membrane
  // traces and the mapper supplies exactly one label per potential DOF.
  template <class Matrix, class Vector, class Grid, class Space, class Material>
  void splitInterfaceGroups(BddcData<Matrix,Vector>& data,
                            Grid const& grid, Space const& space,
                            Material const& material, std::vector<int> const& extraTags,
                            bool mpiEnabled)
  {
    auto const cellLabels = electricalCellLabels(grid,space,material,extraTags);
    auto const groups = findInterfaceComponents(data,cellLabels);
    data.interfaceComponents = interfaceComponentLabels(data,groups);
    size_t count = 0;
    for (auto const& [owners,group] : groups)
    {
      (void)owners;
      count += group.components.size();
    }
    if (mpiRankZero(mpiEnabled))
      std::cout << "BDDC interface splitting: " << groups.size()
                << " owner groups -> " << count << " cell-connected electrical components\n";
  }

  // Multiple original-cell components flag a candidate for enrichment,
  // not a proven convergence defect. Use the same graph for audit and setup.
  template <class Matrix, class Vector, class Grid, class Space, class Material>
  void auditInterfaceGroups(BddcData<Matrix,Vector> const& data,
                            Grid const& grid, Space const& space,
                            Material const& material, std::vector<int> const& extraTags,
                            int interfaceTypes, std::string const& directory,
                            bool mpiEnabled)
  {
    if (!mpiRankZero(mpiEnabled))
      return;
    using Group = std::vector<int>;
    auto const cellLabels = electricalCellLabels(grid,space,material,extraTags);
    auto const groups = findInterfaceComponents(data,cellLabels);

    Kaskade::BDDC::InterfaceAverages<1,int> interfaces(
      data.sharedDofs,data.subdomainSizes,interfaceTypes,data.interfaceComponents);
    auto const constraints = interfaces.coarseConstraints();
    std::set<Group> enabledGroups;
    std::vector<std::vector<size_t>> constraintsBySubdomain(data.tags.size());
    for (size_t id = 0; id < constraints.size(); ++id)
    {
      enabledGroups.insert(constraints[id].first);
      for (int subdomain : constraints[id].first)
        constraintsBySubdomain.at(subdomain).push_back(id);
    }

    // C*1 must agree on every copy of a coarse constraint. This checks the
    // actual constraint rows, not whether their energy-minimizing basis
    // reproduces a near-null vector (that is a separate spectral question).
    std::vector<double> constraintMin(constraints.size(),
                                     std::numeric_limits<double>::infinity());
    std::vector<double> constraintMax(constraints.size(),
                                     -std::numeric_limits<double>::infinity());
    for (size_t subdomain = 0; subdomain < data.tags.size(); ++subdomain)
    {
      auto const C = interfaces.template coarseConstraint<double>(subdomain);
      if (C.N() != constraintsBySubdomain[subdomain].size())
        throw std::runtime_error("BDDC audit: constraint row mapping mismatch");
      for (auto row = C.begin(); row != C.end(); ++row)
      {
        double sum = 0.0;
        for (auto col = row->begin(); col != row->end(); ++col)
          sum += (*col)[0][0];
        size_t const id = constraintsBySubdomain[subdomain][row.index()];
        constraintMin[id] = std::min(constraintMin[id],sum);
        constraintMax[id] = std::max(constraintMax[id],sum);
      }
    }
    double maxConstantMismatch = 0.0;
    size_t zeroConstantConstraints = 0;
    for (size_t id = 0; id < constraints.size(); ++id)
    {
      maxConstantMismatch = std::max(maxConstantMismatch,
                                    std::abs(constraintMax[id]-constraintMin[id]));
      if (constraintMin[id] == 0.0 && constraintMax[id] == 0.0)
        ++zeroConstantConstraints;
    }

    std::filesystem::create_directories(directory);
    std::ofstream groupFile(directory + "/interface_groups.csv");
    std::ofstream componentFile(directory + "/interface_components.csv");
    if (!groupFile || !componentFile)
      throw std::runtime_error("cannot write BDDC interface audit in " + directory);
    groupFile << "group_id,type,enabled,subdomains,material_tags,global_dofs,"
                 "cell_components,electrical_labels\n";
    componentFile << "group_id,component_id,global_dofs,electrical_labels\n";
    size_t groupId = 0, splitGroups = 0, mixedGroups = 0;
    size_t splitEnabled = 0, mixedEnabled = 0, componentConstraints = 0;
    size_t ambiguousDofs = 0;
    for (auto const& [owners,group] : groups)
    {
      auto const& dofs = group.dofs;
      auto const& labels = group.electricalLabels;
      auto const& components = group.components;
      ambiguousDofs += group.ambiguousDofs;
      bool const enabled = enabledGroups.count(owners) != 0;
      if (components.size() > 1)
      {
        ++splitGroups;
        splitEnabled += enabled;
      }
      if (labels.size() > 1)
      {
        ++mixedGroups;
        mixedEnabled += enabled;
      }
      if (enabled)
        componentConstraints += components.size();
      std::vector<int> materialTags;
      for (int subdomain : owners)
        materialTags.push_back(data.tags.at(subdomain));
      char const* type = owners.size() == 2 ? "face"
                       : dofs.size() == 1 ? "corner" : "edge";
      groupFile << groupId << ',' << type << ',' << enabled << ','
                << AuditDetail::joined(owners) << ',' << AuditDetail::joined(materialTags)
                << ',' << dofs.size() << ',' << components.size() << ','
                << AuditDetail::joined(labels) << '\n';
      size_t componentId = 0;
      for (auto const& componentDofs : components)
      {
        std::set<int> componentLabels;
        for (size_t dof : componentDofs)
          for (size_t cell : data.dofCells.at(dof))
            componentLabels.insert(cellLabels.at(cell));
        componentFile << groupId << ',' << componentId++ << ',' << componentDofs.size()
                      << ',' << AuditDetail::joined(componentLabels) << '\n';
      }
      ++groupId;
    }
    std::cout << "BDDC interface audit (electrical label 0 = shared extracellular space):\n"
              << "  groups: " << groups.size() << ", enabled constraints: "
              << constraints.size() << '\n'
              << "  groups with multiple cell-connected components: " << splitGroups
              << " (enabled: " << splitEnabled << ")\n"
              << "  groups mixing electrical labels: " << mixedGroups
              << " (enabled: " << mixedEnabled << ")\n"
              << "  shared DOFs with zero/multiple electrical labels: " << ambiguousDofs << '\n'
              << "  hypothetical constraints after splitting enabled groups: "
              << componentConstraints << '\n'
              << "  maximum disagreement in C*1 across subdomains: "
              << maxConstantMismatch << '\n'
              << "  constraints invisible to a constant field: " << zeroConstantConstraints << '\n'
              << "  group/component details: " << directory << '\n';
  }
}

#endif
