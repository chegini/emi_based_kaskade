#include "mg/bddcRestrictionCapture.hh"
#include <array>
#include <cassert>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
  assert(argc == 2);
  namespace fs = std::filesystem;
  fs::path root(argv[1]);
  assert(!fs::exists(root));
  std::vector<std::array<double,1>> values{{{0.0}},{{-1e-12}},{{2.0}}};
  auto before = values;
  std::vector<bool> excluded{false,false,true};
  std::vector<std::pair<int,std::vector<int>>> interfaces{{2,{0,1}},{3,{1,2}}};
  std::uint64_t invocation=0;
  unsetenv("KASKADE_RESTRICTION_CAPTURE");
  Kaskade::BDDC::captureRestrictionBuffer<1>(7,invocation,values,excluded,interfaces);
  assert(invocation == 0 && !fs::exists(root));
  setenv("KASKADE_RESTRICTION_CAPTURE",argv[1],1);
  for (int k=0; k<50; ++k)
    Kaskade::BDDC::captureRestrictionBuffer<1>(7,invocation,values,excluded,interfaces);
  assert(values == before && invocation == 50);
  assert(fs::exists(root/"sub7-call32.csv"));
  assert(!fs::exists(root/"sub7-call33.csv"));
  assert(fs::exists(root/"sub7-call50.csv"));
  assert(fs::exists(root/"sub7-topology.csv"));
  std::uint64_t duplicate=0;
  bool refused=false;
  try {
    Kaskade::BDDC::captureRestrictionBuffer<1>(7,duplicate,values,excluded,interfaces);
  } catch (std::runtime_error const&) { refused=true; }
  assert(refused && values == before);
  unsetenv("KASKADE_RESTRICTION_CAPTURE");
  std::cout << "Restriction capture checks passed; input buffers unchanged.\n";
}
