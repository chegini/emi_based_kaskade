#include "mg/bddcTransferAudit.hh"
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>

int main(int argc, char** argv)
{
  if (argc!=2) return 2;
  setenv("KASKADE_TRANSFER_AUDIT_DIR",argv[1],1);
  namespace A=Kaskade::BDDC::TransferAudit;
  std::vector<std::uint32_t> values{0,65535,4294967295u};
  A::vector("restrict","send",0,1,"nodal",values);
  A::vector("restrict","send",0,1,"nodal",std::vector<int>{99});
  std::ifstream in(std::filesystem::path(argv[1])/"restrict_send_0_1_nodal.txt");
  std::size_t count=0;in>>count;assert(count==values.size());
  for (auto expected:values) {std::uint64_t actual=0;in>>actual;assert(actual==expected);}
  std::vector<double> v{0.12345678901234567,-1e-120,1e120};
  A::vector("prolongate","recv",1,0,"values",v);
  std::ifstream floats(std::filesystem::path(argv[1])/"prolongate_recv_1_0_values.txt");
  floats>>count;assert(count==v.size());
  for (double expected:v) {double actual=0;floats>>actual;assert(actual==expected);}
  A::wire("restrict","send",0,1,values);
  std::ifstream wire(std::filesystem::path(argv[1])/"restrict_send_0_1_wire.txt");
  wire>>count;assert(count==values.size()*sizeof(std::uint32_t));
  auto raw=reinterpret_cast<unsigned char const*>(values.data());
  for(std::size_t i=0;i<count;++i){int actual=-1;wire>>actual;assert(actual==raw[i]);}
  A::plan<double>("restrict","send",0,1,{0,2},{1},{{0,2}},{{1.0,1.0}});
  for (auto const* stage : {"coarse","offsets","neighbors"})
  {
    std::ifstream plan(std::filesystem::path(argv[1])/
                       (std::string("restrict_send_0_1_")+stage+".txt"));
    std::size_t n=0;int first=-1,last=-1;plan>>n>>first>>last;
    assert(n==2 && first==0 && last==2);
  }
  std::cout<<"Transfer audit capture passed: first observation retained, integer/float round trip, packet bytes.\n";
}
