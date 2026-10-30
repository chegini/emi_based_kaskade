#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "../EMI_bddc_components.hh"

struct LocalDof
{
  int s, i;
};

struct Data
{
  std::vector<std::vector<LocalDof>> sharedDofs;
  std::vector<std::vector<size_t>> localDofs;
  std::vector<std::vector<size_t>> dofCells;
};

int main()
{
  // DOFs 0--2 connect transitively through two extracellular cells. DOF 3
  // is a disconnected extracellular patch with identical owners. DOFs 4--5
  // are duplicated intracellular traces at the same geometric interface;
  // only original-cell incidence is used, so membrane coupling cannot join them.
  Data data;
  data.localDofs = {{0,1,2,3,4,5},{0,1,2,3,4,5},{2,3}};
  data.dofCells = {{0},{0,1},{1},{2},{3},{4}};
  for (int i = 0; i < 6; ++i)
    data.sharedDofs.push_back({{0,i},{1,i}});
  std::vector<int> const cellLabels = {0,0,0,3,5};
  auto const groups = EmiBddc::findInterfaceComponents(data,cellLabels);
  assert(groups.size() == 1);
  auto const& group = groups.at({0,1});
  assert(group.components.size() == 4);
  assert(group.ambiguousDofs == 0);
  assert(group.electricalLabels == std::set<int>({0,3,5}));
  auto const labels = EmiBddc::interfaceComponentLabels(data,groups);
  assert(labels[0] == labels[1] && labels[1] == labels[2]);
  assert(std::set<int>(labels.begin(),labels.end()).size() == 4);

  // Permuting local-copy order and shared-entry traversal preserves the partition.
  Data reordered = data;
  std::reverse(reordered.sharedDofs.begin(),reordered.sharedDofs.end());
  for (auto& copies : reordered.sharedDofs)
    std::reverse(copies.begin(),copies.end());
  auto const otherGroups = EmiBddc::findInterfaceComponents(reordered,cellLabels);
  auto const other = EmiBddc::interfaceComponentLabels(reordered,otherGroups);
  for (size_t i = 0; i < labels.size(); ++i)
    for (size_t j = 0; j < labels.size(); ++j)
      assert((labels[i] == labels[j]) ==
             (other[labels.size()-1-i] == other[labels.size()-1-j]));

  // Owner groups remain separate even when DOFs share an original cell.
  data.sharedDofs[2].push_back({2,0});
  auto const owners = EmiBddc::findInterfaceComponents(data,cellLabels);
  auto const ownerLabels = EmiBddc::interfaceComponentLabels(data,owners);
  assert(owners.size() == 2);
  assert(ownerLabels[1] != ownerLabels[2]);

  // Refuse splitting if a mapper accidentally joins two electrical labels.
  data.dofCells[0].push_back(3);
  bool rejected = false;
  try
  {
    auto const invalid = EmiBddc::findInterfaceComponents(data,cellLabels);
    (void)EmiBddc::interfaceComponentLabels(data,invalid);
  }
  catch (std::runtime_error const&)
  {
    rejected = true;
  }
  assert(rejected);
  Data empty;
  assert(EmiBddc::findInterfaceComponents(empty,{}).empty());
  assert(EmiBddc::interfaceComponentLabels(empty,{}).empty());
  std::cout << "BDDC interface component checks passed\n";
}
