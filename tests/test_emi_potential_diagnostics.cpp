#include "../EMI_potential_diagnostics.hh"
#include "emi_material_test_matrix.hh"
#include <array>
#include <cassert>
#include <sstream>

int main(int argc,char** argv)
{
  assert(argc==2);
  std::string directory=argv[1];
  MaterialTestMatrix A(2);
  A[0][0][0][0]=2; A[0][1][0][0]=-1;
  A[1][0][0][0]=-2; A[1][1][0][0]=3;
  using V=std::vector<std::array<double,1>>;
  V u{{3.141592653589793},{-0.125}},du{{1},{2}},rhs{{3},{4}},zero(2);
  auto before=u;
  EmiBddc::capturePotentialState(false,true,false,directory,0,A,u,du,rhs);
  assert(!std::filesystem::exists(directory+"/potential-state-step0.csv"));
  assert(EmiBddc::writePotentialState(directory,0,A,u,zero,zero));
  assert(EmiBddc::writePotentialState(directory,1,A,u,du,rhs));
  assert(!EmiBddc::writePotentialState(directory,1,A,u,du,rhs));
  assert(u==before && du[0][0]==1 && rhs[0][0]==3);
  std::ifstream op(directory+"/potential-operator.csv");
  std::string line;
  std::getline(op,line); std::getline(op,line); assert(line=="0,1,0,2");
  std::getline(op,line); assert(line=="1,1,2,3");
  std::ifstream state(directory+"/potential-state-step1.csv");
  std::getline(state,line); std::getline(state,line);
  std::stringstream row(line); std::string field;
  std::getline(row,field,','); assert(field=="0");
  std::getline(row,field,','); assert(std::stod(field)==u[0][0]);
  std::getline(row,field,','); assert(field=="1");
  std::getline(row,field,','); assert(field=="3");
  std::getline(row,field,','); assert(field=="3");
}
