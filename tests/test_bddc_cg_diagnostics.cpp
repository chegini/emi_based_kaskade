#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../EMI_bddc_cg_diagnostics.hh"

using Vec = std::vector<std::array<double,1>>;
struct Matrix
{
  void mv(Vec const& x, Vec& y) const
  {
    y[0][0]=2*x[0][0]-x[1][0];
    y[1][0]=-x[0][0]+4*x[1][0]-x[2][0];
    y[2][0]=-x[1][0]+2*x[2][0];
  }
};
struct Sub
{
  Vec x, rr, z, q;
  Vec getSolution() const { return x; }
  Vec const& getRestrictedResidual() const { return rr; }
  Vec const& diagnosticPreconditionedResidual() const { return z; }
  Vec const& diagnosticSearchDirection() const { return q; }
};
struct Storage : std::vector<Sub>
{
  using std::vector<Sub>::vector;
  bool remote = false;
  bool has(size_t s) const { return !(remote && s==1); }
};
struct Data
{
  std::vector<std::vector<size_t>> localDofs{{0,1},{1,2}};
  std::vector<Vec> weights{Vec{{1},{0.25}},Vec{{0.75},{1}}};
};
struct Event
{
  bool restriction;
  int solverCall, pass;
  double rho=0,beta=0,curvature=0;
  bool directionAvailable=false;
};
using Rows = std::vector<std::vector<std::string>>;
Rows read(std::string const& path)
{
  std::ifstream f(path);
  Rows result;
  std::string line,cell;
  while (std::getline(f,line))
  {
    std::stringstream stream(line);
    std::vector<std::string> row;
    while (std::getline(stream,cell,',')) row.push_back(cell);
    if (!line.empty() && line.back()==',') row.push_back("");
    if (row.size()!=21) throw std::runtime_error("Incorrect CSV width");
    result.push_back(row);
  }
  return result;
}
void close(double actual,double expected)
{
  if (!std::isfinite(actual) || std::abs(actual-expected)>1e-13)
    throw std::runtime_error("Incorrect CG diagnostic");
}
int main(int argc,char** argv)
{
  if (argc!=2) return 2;
  std::filesystem::path dir(argv[1]);
  Matrix A;
  Vec rhs{{1},{-2},{3}};
  Data data;
  Storage s{Sub{Vec(2),Vec{{1},{-0.5}},Vec{{1},{1}},Vec{{1},{0}}},
            Sub{Vec(2),Vec{{-1.5},{3}},Vec{{1},{1}},Vec{{0},{0}}}};
  auto original=s;
  auto noSum=[](std::vector<double>&) {};
  auto path=(dir/"diagnostic.csv").string();
  EmiBddc::CgDiagnosticReport report(A,rhs,s,data,true,noSum);
  if (!report.open(path)) return 3;
  report.observe(Event{true,1,1});
  report.observe(Event{false,1,1,2,0,2,true});
  report.finishCall(true);
  assert(s[0].x==original[0].x && s[1].rr==original[1].rr);
  assert(s[0].q==original[0].q && s[1].z==original[1].z);
  // Shared dual contributions must be summed once, without primal weights.
  s[0].rr[1][0]=-0.4; s[1].rr[0][0]=-1.2;
  report.observe(Event{true,2,1});
  s[0].rr=original[0].rr; s[1].rr=original[1].rr;
  report.observe(Event{true,2,2});
  s[0].q=Vec{{0.5},{1}}; s[1].q=Vec{{1},{0}};
  report.observe(Event{false,2,2,3,0.2,3.5,true});
  report.finishCall(false);
  s[0].q=Vec{{0},{1}}; s[1].q=Vec{{1},{0}};
  report.observe(Event{true,3,1});
  report.observe(Event{false,3,1,2,0.1,4,true});
  report.finishCall(true);
  // Zero reference scalar must be reported as undefined, not zero error.
  s[0].z=Vec{{2},{1}}; s[1].z=Vec{{1},{0}};
  report.observe(Event{false,4,1,1,0,0,false});
  report.finishCall(false);
  assert(report.close());
  Rows rows=read(path);
  assert(rows.size()==9);
  close(std::stod(rows[1][5]),0);
  close(std::stod(rows[2][9]),2);
  assert(rows[2][17].empty()); // No previous direction on first call.
  close(std::stod(rows[3][5]),0.4);
  close(std::stod(rows[3][6]),0.4/std::sqrt(14.));
  close(std::stod(rows[4][5]),0);
  close(std::stod(rows[4][7]),0);
  close(std::stod(rows[5][11]),0.5);
  close(std::stod(rows[5][17]),0);
  close(std::stod(rows[7][17]),1/std::sqrt(8.)); // Rejected q was not retained.
  assert(rows[8][11].empty() && rows[8][12]=="1" && rows[8][17].empty());

  // Emulate a second owner contributing the remote subdomain, not a second
  // copy of the full vectors. The root diagnostics must match serial values.
  s=original; s.remote=true;
  bool search=false;
  auto sumRemote=[&](std::vector<double>& v)
  {
    if(v.size()==12)
    {
      v[4]+=-1.5; v[5]+=3;
      if(search) { v[7]+=0.75; v[8]+=1; }
    }
  };
  auto remotePath=(dir/"remote.csv").string();
  EmiBddc::CgDiagnosticReport remote(A,rhs,s,data,true,sumRemote);
  assert(remote.open(remotePath));
  remote.observe(Event{true,1,1});
  search=true;
  remote.observe(Event{false,1,1,2,0,2,true});
  remote.finishCall(true);
  assert(remote.close());
  auto remoteRows=read(remotePath);
  assert(remoteRows[1]==rows[1] && remoteRows[2]==rows[2]);

  // A shared-copy discrepancy is visible even if averaging produces a vector.
  s=original; s[1].x[0][0]=0.1;
  auto jumpPath=(dir/"jump.csv").string();
  EmiBddc::CgDiagnosticReport jump(A,rhs,s,data,true,noSum);
  assert(jump.open(jumpPath));
  jump.observe(Event{true,1,1});
  s[1].x[0][0]=0.2; // Deliberate change between passes must be visible.
  jump.observe(Event{true,1,2});
  assert(jump.close());
  auto jumpRows=read(jumpPath);
  assert(std::stod(jumpRows[1][18])>0);
  close(std::stod(jumpRows[2][7]),0.1);
}
