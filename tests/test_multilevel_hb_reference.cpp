// Independent seven-node reference; no coefficients computed by the inverse under test.
#include "mg/bddcLiftingPlan.hh"
#include <iostream>
#include <stdexcept>
using P=Kaskade::BDDC::IntegerLiftingPlan<double>;
void check(bool v){if(!v)throw std::runtime_error("HB hand-reference mismatch");}
int main(){
 P p;p.nodeCount=7;p.coarse={0,2,4,6};p.fine={1,3,5};
 p.neighbors={{0,2},{2,4},{4,6}};p.weights={{1,1},{1,1},{1,1}};
 std::vector<std::vector<int>> g={{1},{0,2},{1,3},{2,4},{3,5},{4,6},{5}};
 std::vector<double> x={0,1,4,9,16,25,36};
 check(p.floating<1>(x,false)==std::vector<double>({0,-1,4,-1,16,-1,36}));
 p.extendHierarchy(g,0);
 check(p.levels==2 && p.coarse==std::vector<int>({0,4}));
 check(p.fine==std::vector<int>({1,3,5,2,6}));
 check(p.neighbors==std::vector<std::vector<int>>({{0,2},{2,4},{4,6},{0,4},{4}}));
 check(p.weights==std::vector<std::vector<double>>({{1,1},{1,1},{1,1},{1,1},{1}}));
 std::vector<double> expected={0,-1,-4,-1,16,-1,20};
 check(p.floating<1>(x,false)==expected);check(p.floating<1>(expected,true)==x);
 check(p.floating<1>(std::vector<double>(7,3),false)==std::vector<double>({3,0,0,0,3,0,0}));
 // The right endpoint has only one coarse neighbor: a ramp gives nonzero detail.
 check(p.floating<1>({0,1,2,3,4,5,6},false)==std::vector<double>({0,0,0,0,4,0,2}));
 std::cout<<"PASS: hand-computed square/ramp/constant coefficients, both stages, exact inverse, neighborhoods/weights. Boundary ramp detail=2 (one-sided constant prediction).\n";
}
