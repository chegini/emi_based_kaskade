#ifndef EMI_BDDC_COMPONENTS_HH
#define EMI_BDDC_COMPONENTS_HH

#include <algorithm>
#include <cstddef>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace EmiBddc
{
  struct InterfaceGroup
  {
    std::vector<size_t> dofs;
    std::vector<std::vector<size_t>> components;
    std::set<int> electricalLabels;
    size_t ambiguousDofs = 0;
  };

  // Within each owner group, union DOFs sharing an original volume cell.
  // Membrane matrix couplings are deliberately absent from this graph: they
  // connect distinct electrical traces, rather than a continuous potential.
  // All extracellular tags have label zero, so their continuous traces join.
  template <class Data>
  std::map<std::vector<int>,InterfaceGroup> findInterfaceComponents(
    Data const& data, std::vector<int> const& cellLabels)
  {
    std::map<std::vector<int>,InterfaceGroup> groups;
    for (auto const& copies : data.sharedDofs)
    {
      if (copies.size() < 2)
        continue;
      std::vector<int> owners;
      for (auto const& copy : copies)
        owners.push_back(copy.s);
      std::sort(owners.begin(),owners.end());
      auto const& first = copies.front();
      groups[owners].dofs.push_back(data.localDofs.at(first.s).at(first.i));
    }
    for (auto& [owners,group] : groups)
    {
      (void)owners;
      std::vector<size_t> parent(group.dofs.size());
      std::iota(parent.begin(),parent.end(),0);
      auto root = [&](size_t node)
      {
        while (parent[node] != node)
        {
          parent[node] = parent[parent[node]];
          node = parent[node];
        }
        return node;
      };
      std::unordered_map<size_t,size_t> firstDofInCell;
      for (size_t local = 0; local < group.dofs.size(); ++local)
      {
        std::set<int> dofLabels;
        for (size_t cell : data.dofCells.at(group.dofs[local]))
        {
          int const label = cellLabels.at(cell);
          dofLabels.insert(label);
          group.electricalLabels.insert(label);
          auto const [it,inserted] = firstDofInCell.emplace(cell,local);
          if (!inserted)
            parent[root(local)] = root(it->second);
        }
        if (dofLabels.size() != 1)
          ++group.ambiguousDofs;
      }
      std::map<size_t,std::vector<size_t>> components;
      for (size_t local = 0; local < group.dofs.size(); ++local)
        components[root(local)].push_back(group.dofs[local]);
      for (auto& [representative,dofs] : components)
      {
        (void)representative;
        group.components.push_back(std::move(dofs));
      }
    }
    return groups;
  }

  // Labels are indexed by sharedDofs entry, not global DOF. They must be
  // identical on every MPI rank, so assignment uses ordered owner/component
  // traversal of the replicated partition metadata.
  template <class Data>
  std::vector<int> interfaceComponentLabels(
    Data const& data, std::map<std::vector<int>,InterfaceGroup> const& groups)
  {
    std::vector<int> byGlobal(data.dofCells.size(),-1);
    int next = 0;
    for (auto const& [owners,group] : groups)
    {
      (void)owners;
      if (group.ambiguousDofs != 0)
        throw std::runtime_error("BDDC interface splitting requires one electrical label per DOF");
      for (auto const& dofs : group.components)
      {
        for (size_t dof : dofs)
          byGlobal.at(dof) = next;
        ++next;
      }
    }
    std::vector<int> labels(data.sharedDofs.size(),-1);
    for (size_t i = 0; i < data.sharedDofs.size(); ++i)
    {
      auto const& copies = data.sharedDofs[i];
      if (copies.size() < 2)
        continue;
      auto const& first = copies.front();
      labels[i] = byGlobal.at(data.localDofs.at(first.s).at(first.i));
      if (labels[i] < 0)
        throw std::runtime_error("BDDC shared DOF has no interface component");
    }
    return labels;
  }
}

#endif
