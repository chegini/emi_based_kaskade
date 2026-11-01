#!/usr/bin/env python3
"""Exercise the production shadow method/report on a known two-subdomain system."""
import argparse
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
from test_bddc_cg_observer import block

HARNESS = r'''
#include <cassert>
#include <array>
#include <vector>
#include "EMI_bddc_prolongation_diagnostics.hh"
struct Entry {
 double v=0;
 Entry()=default; Entry(double x):v(x) {}
 operator double() const { return v; }
 Entry& operator=(double x) { v=x;return *this; }
 double& operator[](int) { return v; }
 double const& operator[](int) const { return v; }
 bool operator==(Entry const&) const=default;
};
struct Vec:std::vector<Entry> {
 using std::vector<Entry>::vector;
 Vec& operator=(double x) { for(auto& v:*this)v=x;return *this; }
 Vec& operator+=(Vec const& x) { for(size_t i=0;i<size();++i)(*this)[i]=(*this)[i].v+x[i].v;return *this; }
};
struct LocalMatrix {
 void smv(double a,Vec const& x,Vec& y) const {
   y[0]=a*(2*x[0].v-x[1].v);y[1]=a*(-x[0].v+2*x[1].v);
 }
};
namespace Kaskade::BDDC {
template<int m,class S,class C,class T> struct Subdomain {
 using Vector=Vec;
 int nx=2;
 Vec rawCorrection,x,rr,z;
 std::vector<int> interfaceDofs,globalIndex_notActive;
 LocalMatrix A;
 void solveDirichlet(Vec& delta,Vec const& rhs) const {
   int i=1-interfaceDofs[0];delta[i]=rhs[i].v/2;
 }
 Vec diagnosticExactProlongation(Vec const&) const;
 Vec getSolution() const { return x; }
 Vec const& getRestrictedResidual() const { return rr; }
 Vec const& diagnosticRawCorrection() const { return rawCorrection; }
 Vec diagnosticProlongatedCorrection() const { return z; }
 Vec diagnosticProlongationWeights() const {
   Vec w(2);w=1;w[interfaceDofs[0]]=0.5;return w;
 }
 // Controlled producer samples test report assembly independently of the codec.
 Subdomain const& transfer() const {return *this;}
 struct Sample {
   Vec contribution;
   size_t groups=1,metadataBytes=8,streamBytes=72,paddingBytes=0,exactBytes=0;
   size_t totalBytes=80,packedBytes=80,messages=1;
   bool codecValid=true;
 };
 Sample diagnosticGroupedProlongation(size_t n,int bits,size_t group) const {
   Sample s;s.contribution=Vec(n);s.contribution=0;
   auto k=interfaceDofs[0];
   s.contribution[k]=(bits==8 && group==0) ? z[k].v/2 : rawCorrection[k].v/2;
   if (group) s.contribution[k].v+=0.05;
   if (bits==16) {s.totalBytes=100;s.streamBytes=92;}
   return s;
 }
};
METHOD
}
using Sub=Kaskade::BDDC::Subdomain<1,double,double,int>;
struct Storage:std::vector<Sub> {
 using std::vector<Sub>::vector;
 bool has(size_t) const { return true; }
};
struct Matrix {
 void mv(Vec const& x,Vec& y) const {
   y[0]=2*x[0].v-x[1].v;y[1]=-x[0].v+4*x[1].v-x[2].v;y[2]=-x[1].v+2*x[2].v;
 }
};
struct TwoChains {
 void mv(Vec const& x,Vec& y) const {
   for(int j: {0,3}) {
     y[j]=2*x[j].v-x[j+1].v;
     y[j+1]=-x[j].v+4*x[j+1].v-x[j+2].v;
     y[j+2]=-x[j+1].v+2*x[j+2].v;
   }
 }
};
struct Data {
 std::vector<std::vector<size_t>> localDofs{{0,1},{1,2}};
 // Deliberately different from actual prolongation weights: projection only.
 std::vector<Vec> weights{Vec{1,0.25},Vec{0.75,1}};
};
struct Event { bool restriction; int solverCall,pass; double rho=0; };
template<class T> struct Sparse {
 std::vector<std::pair<size_t,T>> entries;
 struct It {
   typename std::vector<std::pair<size_t,T>>::const_iterator p;
   size_t index() const {return p->first;}
   T const& operator*() const {return p->second;}
   T const* operator->() const {return &p->second;}
   It& operator++() {++p;return *this;}
   bool operator!=(It const& b) const {return p!=b.p;}
 };
 It begin() const {return {entries.begin()};} It end() const {return {entries.end()};}
 size_t N() const {return entries.size();}
};
struct Interfaces {
 using Block=std::array<std::array<double,1>,1>;
 using Row=Sparse<Block>;
 std::vector<Sparse<Row>> matrices;
 auto coarseConstraints() const {return std::vector<std::pair<std::vector<int>,int>>{{{0,1},1}};}
 template<class T> auto coarseConstraint(size_t s) const {return matrices[s];}
};
int main(int argc,char** argv) {
 assert(argc==2);
 Storage s(2);
 s[0].rawCorrection=Vec{4,2};s[1].rawCorrection=Vec{4,4};
 s[0].interfaceDofs={1};s[1].interfaceDofs={0};
 s[0].x=Vec{0,0};s[1].x=Vec{0,0};
 s[0].rr=Vec{1,-1};s[1].rr=Vec{-1,1};
 s[0].z=Vec{6,6};s[1].z=Vec{6,5};
 auto saved=s;
 Vec rhs{1,-2,1};Matrix A;Data d;int reductions=0;
 Interfaces ifs;
 Interfaces::Row row0,row1;
 row0.entries.push_back({1,Interfaces::Block{{{1}}}});
 row1.entries.push_back({0,Interfaces::Block{{{1}}}});
 Sparse<Interfaces::Row> c0,c1;
 c0.entries.push_back({0,row0});c1.entries.push_back({0,row1});ifs.matrices={c0,c1};
 auto topology=EmiBddc::prolongationConstraintGroups(ifs,d);
 assert(topology.size()==1 && topology[0].dofs==std::vector<size_t>{1} && topology[0].owners==2);
 ifs.matrices[1].entries[0].second.entries[0].second[0][0]=2;
 bool unsupported=false;
 try {EmiBddc::prolongationConstraintGroups(ifs,d);}
 catch(std::runtime_error const&) {unsupported=true;}
 assert(unsupported);
 EmiBddc::ProlongationDiagnosticReport report(A,rhs,s,d,true,
     [&](std::vector<double>&){++reductions;});
 assert(report.open(argv[1]));
 report.constraintProbe(topology,1,std::string(argv[1])+"-probe");
 report.observe(Event{true,1,1});assert(reductions==0);
 report.observe(Event{false,1,1,-1});assert(reductions==4);
 for(int i=0;i<2;++i) {
   assert(s[i].rawCorrection==saved[i].rawCorrection && s[i].x==saved[i].x);
   assert(s[i].rr==saved[i].rr && s[i].z==saved[i].z);
 }
 s[0].z=Vec{4.5,3};s[1].z=Vec{3,3.5};
 report.observe(Event{false,2,1,2});assert(reductions==6);
 double a=3.0000001;
 s[0].z=Vec{4.5+a/2,3+a};s[1].z=Vec{3+a,3.5+a/2};
 report.constraintProbe({{{1},2}},3,std::string(argv[1])+"-round");
 report.observe(Event{false,3,1,2-a});assert(reductions==10);
 bool rejected=false;
 try {report.constraintProbe({{{1},2},{{1},2}},4,"bad");}
 catch(std::runtime_error const&) {rejected=true;}
 assert(rejected);
 assert(report.close());
 Storage grouped=saved;
 EmiBddc::ProlongationDiagnosticReport grouping(A,rhs,grouped,d,true,
     [](std::vector<double>&){});
 assert(grouping.open(std::string(argv[1])+"-group-observer.csv"));
 grouping.groupingProbe(1,8,std::string(argv[1])+"-grouping.csv");
 grouping.observe(Event{false,1,1,-1});assert(grouping.close());
 for(int i=0;i<2;++i) assert(grouped[i].z==saved[i].z && grouped[i].rawCorrection==saved[i].rawCorrection);
 // Inactive interfaces match the production mask and do not leak into the solve.
 s[0].globalIndex_notActive={1};
 auto masked=s[0].diagnosticExactProlongation(Vec{0,99});
 assert(masked[0].v==4 && masked[1].v==0);
 // A multi-node group: remove the mean but retain the within-group error.
 Storage multi(4);
 for(int j: {0,2}) {
   multi[j]=saved[0];multi[j+1]=saved[1];
   double e=j==0?3:1;
   multi[j].z=Vec{4.5+e/2,3+e};multi[j+1].z=Vec{3+e,3.5+e/2};
 }
 auto original=multi;
 Data md;md.localDofs={{0,1},{1,2},{3,4},{4,5}};
 md.weights={Vec{1,0.25},Vec{0.75,1},Vec{1,0.25},Vec{0.75,1}};
 Vec mrhs{1,-2,1,1,-3,1};TwoChains ma;
 EmiBddc::ProlongationDiagnosticReport multireport(ma,mrhs,multi,md,true,
   [](std::vector<double>&){});
 assert(multireport.open(std::string(argv[1])+"-multi.csv"));
 multireport.constraintProbe({{{1,4},2}},1,std::string(argv[1])+"-multi");
 multireport.observe(Event{false,1,1,-4});assert(multireport.close());
 for(int i=0;i<4;++i) {
   assert(multi[i].z==original[i].z && multi[i].rawCorrection==original[i].rawCorrection);
 }
}
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kaskade-root',type=Path,default=Path('/kaskade7'))
    args=p.parse_args()
    source=(args.kaskade_root/'mg/bddc.hpp').read_text()
    method=block(source,'  template <int m, class Scalar_, class CoarseScalar, class Transfer>\n'
                        '  auto Subdomain<m,Scalar_,CoarseScalar,Transfer>::diagnosticExactProlongation(')
    repo=Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix='emi-prolongation-diagnostic-') as tmp:
        root=Path(tmp);cpp=root/'test.cpp';exe=root/'test';out=root/'report.csv'
        cpp.write_text(HARNESS.replace('METHOD',method))
        subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                       ['-std=c++20','-O2','-I'+str(repo),str(cpp),'-o',str(exe)],check=True)
        subprocess.run([str(exe),str(out)],check=True)
        rows=list(csv.DictReader(out.open()))
        assert len(rows)==3
        expected={'rho_before_p':2,'rho_after_p':-1,'rho_change_p':-3,
                  'rho_closure_error':0,'restriction_error_rel':0,'p_error_energy':27,
                  'actual_z_jump':0,'weight_partition_defect':0,'inactive_global_dofs':0,'rho_internal':-1}
        for key,want in expected.items():
            assert math.isclose(float(rows[0][key]),want,abs_tol=1e-12),(key,rows[0][key],want)
        assert math.isclose(float(rows[0]['p_error_norm']),math.sqrt(13.5))
        assert float(rows[1]['p_error_norm'])==0 and float(rows[1]['rho_before_p'])==2
        probe=list(csv.DictReader(Path(str(out)+'-probe-summary.csv').open()))[0]
        for key,want in {'rho_mean_corrected':2,'rho_float32_mean_corrected':2,
                         'attributed_rho_error':-3,'unattributed_rho_error':0,
                         'max_remaining_mean_error':0,'estimated_payload32_per_call':8,
                         'estimated_payload64_per_call':16}.items():
            assert math.isclose(float(probe[key]),want,abs_tol=1e-12),(key,probe[key])
        rounded=list(csv.DictReader(Path(str(out)+'-round-summary.csv').open()))[0]
        assert math.isclose(float(rounded['rho_mean_corrected']),2,abs_tol=1e-12)
        assert math.isclose(float(rounded['rho_float32_mean_corrected']),2-1e-7,abs_tol=1e-12)
        assert math.isclose(float(rounded['max_remaining_float32_mean_error']),1e-7,abs_tol=1e-12)
        grouping=list(csv.DictReader(Path(str(out)+'-grouping.csv').open()))
        assert len(grouping)==9
        for row in grouping:
            bits=int(row['bits']);group=int(row['group_size'])
            expected=-1 if bits==8 and group==0 else 2 if group==0 else 1.9
            assert math.isclose(float(row['rho_candidate']),expected,abs_tol=1e-12),row
            assert int(row['modeled_p_bytes'])==(200 if bits==16 else 160)
            if int(row['current_precision']): assert float(row['relative_difference_vs_actual'])<1e-14
        multi=list(csv.DictReader(Path(str(out)+'-multi-summary.csv').open()))[0]
        for key,want in {'rho_before_p':1,'rho_after_p':-4,'rho_mean_corrected':2,
                         'attributed_rho_error':-5,'unattributed_rho_error':0,
                         'max_remaining_mean_error':0}.items():
            assert math.isclose(float(multi[key]),want,abs_tol=1e-12),(key,multi[key])
    print('Production prolongation shadow checks passed: harmonic extension, rho decomposition, actual weights, state unchanged.')


if __name__=='__main__':main()
