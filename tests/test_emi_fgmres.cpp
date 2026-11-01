#include "EMI_fgmres_mpi.hh"
#include <cassert>
#include <iostream>

using namespace EmiBddc::Fgmres;
void require(bool value) {if (!value) throw std::runtime_error("FGMRES numerical regression");}

int main()
{
  auto identity=[](Vector const& x) {return x;};
  auto negative=[](Vector x) {for (auto& a:x) a=-a;return x;};
  auto observe=[](int,int,double,double) {};
  Vector b{1,2,3};
  auto exact=solve(b,30,50,1e-12,identity,negative,observe);
  require(exact.status==Status::converged && exact.iterations==1);
  for (size_t i=0;i<b.size();++i) require(std::abs(exact.x[i]-b[i])<1e-12);

  auto matrix=[](Vector const& x) {
    return Vector{4*x[0]+x[1],-x[0]+3*x[1]+x[2],2*x[1]+5*x[2]};
  };
  Vector known{1,-2,3}; b=matrix(known);
  int applications=0;
  auto variable=[&](Vector x) {
    ++applications;
    for (size_t i=0;i<x.size();++i)
      x[i]*=(applications%2 ? -1.0 : 1.0)*(i+applications+1);
    return x;
  };
  auto flex=solve(b,3,10,1e-12,matrix,variable,observe);
  require(flex.status==Status::converged && applications==flex.iterations);
  for (size_t i=0;i<known.size();++i) require(std::abs(flex.x[i]-known[i])<1e-10);

  auto diagonal=[](Vector x) {for (size_t i=0;i<x.size();++i) x[i]*=i+1;return x;};
  auto restarted=solve(Vector{1,1,1},1,150,1e-10,diagonal,identity,observe);
  require(restarted.status==Status::converged && restarted.cycles>1);
  auto limited=solve(Vector{1,1,1},1,1,1e-12,diagonal,identity,observe);
  require(limited.status==Status::iterationLimit && limited.iterations==1);
  auto zero=solve(Vector(3,0),30,50,1e-7,identity,identity,observe);
  require(zero.status==Status::converged && zero.iterations==0);
  auto breakdown=solve(b,30,50,1e-7,matrix,[](Vector x) {return Vector(x.size(),0);},observe);
  require(breakdown.status==Status::breakdown && breakdown.iterations==1);
  auto bad=solve(b,30,50,1e-7,matrix,[](Vector x) {x[0]=NAN;return x;},observe);
  require(bad.status==Status::nonfinite);

  // An inaccurate/changing matvec must not pass merely because the Arnoldi
  // least-squares residual is zero. The explicit true residual rejects it.
  int mv=0;
  auto inconsistent=[&](Vector x) {++mv;for (auto& a:x) a*=mv%2 ? 1 : 2;return x;};
  auto falseHappy=solve(Vector{1,1},3,3,1e-7,inconsistent,identity,observe);
  require(falseHappy.status==Status::breakdown && falseHappy.residual>0.5);

  EmiBddc::OuterSolverTraffic traffic;
  auto serial=EmiBddc::rootFgmres(b,false,30,50,1e-12,matrix,negative,observe,traffic);
  require(serial.status==Status::converged && traffic.bytes==0 && traffic.collectives==0);
  traffic.add(100,8);traffic.add(100,8,true);
  require(traffic.bytes==2100 && traffic.collectives==2);
  std::cout << "FGMRES checks passed: negative/variable preconditioning, restart, true residual, "
               "zero RHS, breakdown/nonfinite guards, root adapter and payload accounting.\n";
}
