#include "EMI_partitioning.hh"

#include <iostream>

using namespace EmiBddc::Partitioning;

void require(bool ok, char const* reason)
{
  if (!ok) throw std::runtime_error(reason);
}

template<class F>
void rejected(F&& f)
{
  bool bad=false;
  try { f(); } catch (std::exception const&) { bad=true; }
  require(bad,"Invalid partition input was accepted");
}

int main()
{
  Problem p;
  for (int i=0; i<8; ++i)
  {
    p.tags.push_back(100+i);
    p.intracellular.push_back(i%2);
    p.localDofs.push_back(10);
  }
  p.edges={{0,1,100},{2,3,100},{4,5,100},{6,7,100},
           {1,2,1},{3,4,1},{5,6,1},{0,7,1}};
  // Synthetic nonconsecutive shared groups test policy decisions, not hardware timings.
  auto nodes=nodeIds({0,1,0,1});
  require(nodes==std::vector<int>({0,1,0,1}),"Node IDs assume adjacent ranks");
  auto b=baseline(8,4), c=coupled(p,nodes);
  validateOwners(p,c,nodes);
  require(b==std::vector<int>({0,1,2,3,0,1,2,3}),"P2 ownership changed");
  require(cost(p,c,nodes)<cost(p,b,nodes),"P1 failed to colocate strong I/E pairs");
  auto cb=loads(p,b,nodes), cc=loads(p,c,nodes);
  require(*std::max_element(cc.rank.begin(),cc.rank.end())<=
          *std::max_element(cb.rank.begin(),cb.rank.end()),"P1 rank load regressed");
  require(*std::max_element(cc.node.begin(),cc.node.end())<=
          *std::max_element(cb.node.begin(),cb.node.end()),"P1 node load regressed");
  require(c==coupled(p,nodes),"P1 is nondeterministic");
  auto pools=separatedPools(p,nodes,"0","1");
  auto s=separated(p,nodes,pools);
  validateSeparated(p,s,nodes,pools);
  for (auto const& e : p.edges)
    if (p.intracellular[e.a]!=p.intracellular[e.b])
      require(nodes[s[e.a]]!=nodes[s[e.b]],"P3 leaves an I/E interface node-local");
  auto complement=separatedPools(p,nodes,"0","");
  require(complement.extra==std::vector<int>({1}),"Node pool complement is wrong");
  auto automatic=separatedPools(p,nodes,"","");
  validateSeparated(p,separated(p,nodes,automatic),nodes,automatic);
  rejected([&] { separatedPools(p,nodes,"0","0"); });
  rejected([&] { separatedPools(p,{0,0,0,0},"",""); });
  rejected([&] { pool("0,0",2); });
  rejected([&] { pool("1,",2); });
  rejected([&] { pool("2",2); });
  rejected([&] { nodeIds({0,3}); });
  rejected([&] { coupled(p,nodes,.9); });
  rejected([&] { validateOwners(p,{0,1},nodes); });
  auto bad=b;bad[0]=4;rejected([&] { validateOwners(p,bad,nodes); });
  rejected([&] { validateSeparated(p,b,nodes,pools); });
  require(mode("P1")=="coupled" && mode("P2")=="baseline" &&
          mode("P3")=="separated","Policy aliases differ");
  rejected([&] { mode("unknown"); });
  // Inputs are immutable, and weights are full local DOFs, not interface sizes.
  require(p.edges[0].dofs==100 && p.localDofs[0]==10,"Planner changed problem data");
  std::cout<<"PASS: P1/P2/P3, nonconsecutive nodes, load bounds, node pools and invalid maps\n";
}
