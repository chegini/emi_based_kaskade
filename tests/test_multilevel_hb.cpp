#include "mg/bddcLiftingPlan.hh"
#include <iostream>
#include <random>
#include <numeric>
using Plan=Kaskade::BDDC::IntegerLiftingPlan<double>;
void check(bool ok) { if (!ok) throw std::runtime_error("Multilevel HB regression"); }
template<class F> void rejects(F f)
{
  bool rejected=false; try {f();} catch(std::runtime_error const&) {rejected=true;}
  check(rejected);
}
Plan path(int n)
{
  Plan p; p.nodeCount=n;
  for(int i=0;i<n;i+=2) p.coarse.push_back(i);
  for(int i=1;i<n;i+=2)
  {
    p.fine.push_back(i); p.neighbors.push_back({i-1});
    if(i+1<n) p.neighbors.back().push_back(i+1);
    p.weights.emplace_back(p.neighbors.back().size(),1.0);
  }
  return p;
}
std::vector<std::vector<int>> graph(int n)
{
  std::vector<std::vector<int>> g(n);
  for(int i=1;i<n;++i) {g[i-1].push_back(i);g[i].push_back(i-1);}
  return g;
}
int main()
{
  auto old=path(31); auto bytes=old.pack(); auto single=old;
  single.extendHierarchy(graph(31),1); check(single.pack()==bytes);
  auto multi=old; multi.extendHierarchy(graph(31),0);
  check(multi.levels==4 && multi.coarse.size()==2 && multi.fine.size()==29);
  auto two=old;two.extendHierarchy(graph(31),2);check(two.coarse.size()==8);
  std::mt19937 random(42);
  for(int n:{0,1,2,3,4,7,16,31,32,127})
    for(unsigned depth:{0u,1u,2u,3u,64u})
    {
      auto p=path(n);p.extendHierarchy(graph(n),depth);
      auto restored=Plan::unpack(p.pack(),n);
      for(bool zigzag:{false,true}) for(int trial=0;trial<20;++trial)
      {
        std::vector<std::uint32_t> u(n*2);
        for(auto& v:u) v=random();
        if(n>0) {u[0]=0;u[1]=0xffffffffu;}
        check(restored.inverse<2>(p.forward<2>(u,zigzag),zigzag)==u);
        std::vector<double> x(n*2);
        for(auto& v:x) v=double(int(random()%20000)-10000)/7;
        auto y=restored.floating<2>(p.floating<2>(x,false),true);
        for(std::size_t i=0;i<x.size();++i) check(std::abs(x[i]-y[i])<1e-10);
      }
      auto x=std::vector<double>(n,3.25);auto y=p.floating<1>(x,false);
      for(int i:p.fine) check(y[i]==0);
    }
  auto packed=multi.pack();
  for(std::size_t i=0;i<packed.size();++i)
    rejects([&]{Plan::unpack({packed.begin(),packed.begin()+i},31);});
  auto bad=multi;bad.neighbors.back()[0]=bad.fine.front();rejects([&]{bad.pack();});
  bad=multi;bad.levels=1;rejects([&]{bad.pack();});
  // A disconnected singleton never predicts or is predicted across components.
  auto disconnected=path(7);disconnected.nodeCount=8;disconnected.coarse.push_back(7);
  auto g=graph(7);g.emplace_back();disconnected.extendHierarchy(g,0);
  check(std::find(disconnected.coarse.begin(),disconnected.coarse.end(),7)!=disconnected.coarse.end());
  for(auto const& neighbors:disconnected.neighbors)
    check(std::find(neighbors.begin(),neighbors.end(),7)==neighbors.end());
  // Demonstrate coefficient concentration on a smooth path, without claiming
  // any achieved solver compression ratio or general first-moment property.
  std::vector<double> smooth(31);
  for(int i=0;i<31;++i) smooth[i]=2+0.125*i;
  auto a=old.floating<1>(smooth,false),b=multi.floating<1>(smooth,false);
  auto zeros=[](auto const& x){return std::count(x.begin(),x.end(),0.0);};
  check(zeros(b)>zeros(a));
  std::cout<<"PASS: legacy unchanged; 31 -> 16 -> 8 -> 4 -> 2; exact modular codes; "
           <<"floating roundtrip; components; malformed plans. Smooth zeros "
           <<zeros(a)<<" -> "<<zeros(b)<<" of 31.\n";
}
