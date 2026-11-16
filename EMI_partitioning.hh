#ifndef EMI_PARTITIONING_HH
#define EMI_PARTITIONING_HH

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifdef __linux__
#include <sched.h>
#endif

#include "mg/bddcMpi.hh"
#include "mg/bddcMpiRouting.hh"

namespace EmiBddc::Partitioning
{
  struct Edge
  {
    int a, b;
    std::uint64_t dofs;
  };

  /** Material regions and actual algebraic interface incidence, not myocyte IDs.
   * Vertex weights are full overlapped local DOFs, available before owner assembly.
   * Matrix nonzeros are reported later; DOF weights are not measured solve costs.
   */
  struct Problem
  {
    std::vector<int> tags;
    std::vector<bool> intracellular;
    std::vector<std::uint64_t> localDofs;
    std::vector<Edge> edges;
  };

  inline std::string mode(std::string name)
  {
    if (name == "P1") name = "coupled";
    if (name == "P2") name = "baseline";
    if (name == "P3") name = "separated";
    if (name != "baseline" && name != "coupled" && name != "separated")
      throw std::invalid_argument("bddcPartitionMode must be baseline/P2, coupled/P1 or separated/P3");
    return name;
  }

  /** Dense logical node IDs derived from actual MPI shared-memory leaders. */
  inline std::vector<int> nodeIds(std::vector<int> const& leaders)
  {
    if (leaders.empty()) throw std::invalid_argument("Empty rank/node topology");
    auto unique = leaders;
    std::sort(unique.begin(), unique.end());
    unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    std::vector<int> result;
    for (int leader : leaders)
    {
      if (leader < 0 || leader >= static_cast<int>(leaders.size()) || leaders[leader] != leader)
        throw std::invalid_argument("Invalid shared-memory leader vector");
      result.push_back(std::lower_bound(unique.begin(), unique.end(), leader) - unique.begin());
    }
    return result;
  }

  inline std::vector<int> baseline(std::size_t subdomains, int ranks)
  {
    if (ranks <= 0 || subdomains < static_cast<std::size_t>(ranks))
      throw std::invalid_argument("Partition needs at least one subdomain per MPI rank");
    std::vector<int> owners(subdomains);
    for (std::size_t i=0; i<subdomains; ++i) owners[i] = i % ranks;
    return owners;
  }

  inline std::vector<int> pool(std::string const& text, int nodes)
  {
    std::vector<int> result;
    if (text.empty()) return result;
    std::stringstream input(text);
    std::string token;
    if (text.back() == ',') throw std::invalid_argument("Trailing comma in partition node pool");
    while (std::getline(input, token, ','))
    {
      std::size_t end=0;
      int n=std::stoi(token, &end);
      if (end != token.size() || n < 0 || n >= nodes ||
          std::find(result.begin(), result.end(), n) != result.end())
        throw std::invalid_argument("Invalid or repeated partition node index");
      result.push_back(n);
    }
    std::sort(result.begin(), result.end());
    return result;
  }

  inline void validateProblem(Problem const& p)
  {
    if (p.tags.empty() || p.tags.size() != p.intracellular.size() ||
        p.tags.size() != p.localDofs.size() ||
        std::set<int>(p.tags.begin(),p.tags.end()).size() != p.tags.size())
      throw std::invalid_argument("Invalid material partition metadata");
    for (auto w : p.localDofs)
      if (!w) throw std::invalid_argument("Empty material subdomain in partition metadata");
    std::set<std::pair<int,int>> pairs;
    for (auto const& e : p.edges)
      if (e.a < 0 || e.a >= e.b || e.b >= static_cast<int>(p.tags.size()) || !e.dofs ||
          !pairs.emplace(e.a,e.b).second)
        throw std::invalid_argument("Invalid or repeated partition interface edge");
  }

  inline void validateOwners(Problem const& p, std::vector<int> const& owners,
                             std::vector<int> const& nodes)
  {
    validateProblem(p);
    if (owners.size() != p.tags.size())
      throw std::invalid_argument("Ownership map does not cover every material subdomain");
    for (int rank : owners)
      if (rank < 0 || rank >= static_cast<int>(nodes.size()))
        throw std::invalid_argument("Ownership map contains an out-of-range MPI rank");
  }

  struct Loads
  {
    std::vector<std::uint64_t> rank, node;
  };

  inline Loads loads(Problem const& p, std::vector<int> const& owners,
                     std::vector<int> const& nodes)
  {
    Loads result{std::vector<std::uint64_t>(nodes.size()),
                 std::vector<std::uint64_t>(*std::max_element(nodes.begin(),nodes.end())+1)};
    for (std::size_t s=0; s<owners.size(); ++s)
    {
      result.rank.at(owners[s]) += p.localDofs[s];
      result.node.at(nodes.at(owners[s])) += p.localDofs[s];
    }
    return result;
  }

  /// Lexicographic objective: prefer I/E colocation, then minimize all byte-weighted cuts.
  inline auto cost(Problem const& p, std::vector<int> const& owners,
                   std::vector<int> const& nodes)
  {
    std::uint64_t ie=0, total=0;
    for (auto const& e : p.edges)
      if (nodes[owners[e.a]] != nodes[owners[e.b]])
      {
        total += e.dofs;
        if (p.intracellular[e.a] != p.intracellular[e.b]) ie += e.dofs;
      }
    return std::make_pair(ie,total);
  }

  inline std::vector<int> decreasingDofs(Problem const& p)
  {
    std::vector<int> ids(p.tags.size());
    std::iota(ids.begin(),ids.end(),0);
    std::stable_sort(ids.begin(),ids.end(),[&](int a,int b)
    { return p.localDofs[a] > p.localDofs[b]; });
    return ids;
  }

  /** Select a least-loaded node, then a rank on it. Unused eligible ranks are
   * filled first, avoiding a heavy-pair placement on one node just because rank
   * numbers are adjacent. All decisions use the supplied discovered topology.
   */
  inline int leastLoaded(std::vector<int> const& allowed, std::vector<int> const& nodes,
                         Loads const& work)
  {
    bool unused=false;
    for (int r : allowed) unused = unused || work.rank[r] == 0;
    int best=-1;
    for (int r : allowed)
    {
      if (unused && work.rank[r]) continue;
      if (best < 0 || std::tie(work.node[nodes[r]],work.rank[r],r) <
                      std::tie(work.node[nodes[best]],work.rank[best],best)) best=r;
    }
    if (best < 0) throw std::invalid_argument("Partition node pool has no MPI ranks");
    return best;
  }

  inline void assign(Problem const& p, int s, int r, std::vector<int>& owners,
                     std::vector<int> const& nodes, Loads& work)
  {
    owners[s]=r;
    work.rank[r]+=p.localDofs[s];
    work.node[nodes[r]]+=p.localDofs[s];
  }

  /** Deterministic bounded heuristic, not a global graph optimum. Start from
   * strongest I/E bundles and refine by owner swaps. Accept only placements
   * within the requested baseline max-rank/max-node DOF load bounds, with all
   * ranks active, and no worse weighted objective than the baseline.
   */
  inline std::vector<int> coupled(Problem const& p, std::vector<int> const& nodes,
                                  double balanceFactor=1.0)
  {
    if (!std::isfinite(balanceFactor) || balanceFactor < 1)
      throw std::invalid_argument("Partition balance factor must be finite and >=1");
    auto original=baseline(p.tags.size(),nodes.size());
    validateOwners(p,original,nodes);
    auto base=loads(p,original,nodes);
    auto maxRank=*std::max_element(base.rank.begin(),base.rank.end());
    auto maxNode=*std::max_element(base.node.begin(),base.node.end());
    auto admissible=[&](std::vector<int> const& owners)
    {
      auto w=loads(p,owners,nodes);
      return std::find(w.rank.begin(),w.rank.end(),0)==w.rank.end() &&
             *std::max_element(w.rank.begin(),w.rank.end()) <= balanceFactor*maxRank &&
             *std::max_element(w.node.begin(),w.node.end()) <= balanceFactor*maxNode;
    };
    std::vector<int> parent(p.tags.size());
    std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int v) { while (parent[v] != v) v=parent[v]; return v; };
    auto weights=p.localDofs;
    auto edges=p.edges;
    std::stable_sort(edges.begin(),edges.end(),[](Edge const& a,Edge const& b)
    { return a.dofs > b.dofs; });
    std::size_t components=p.tags.size();
    for (auto const& e : edges)
    {
      if (p.intracellular[e.a] == p.intracellular[e.b]) continue;
      int a=root(e.a), b=root(e.b);
      if (a != b && components > base.node.size() &&
          weights[a]+weights[b] <= balanceFactor*maxNode)
      { parent[b]=a; weights[a]+=weights[b]; --components; }
    }
    std::map<int,std::vector<int>> grouped;
    for (int s : decreasingDofs(p)) grouped[root(s)].push_back(s);
    std::vector<std::vector<int>> bundles;
    for (auto const& [r,ids] : grouped) { (void)r; bundles.push_back(ids); }
    std::stable_sort(bundles.begin(),bundles.end(),[&](auto const& a,auto const& b)
    { return weights[root(a.front())] > weights[root(b.front())]; });
    std::vector<int> candidate(p.tags.size(),-1);
    Loads work{std::vector<std::uint64_t>(nodes.size()),
               std::vector<std::uint64_t>(base.node.size())};
    for (auto const& bundle : bundles)
    {
      int n=std::min_element(work.node.begin(),work.node.end())-work.node.begin();
      std::vector<int> ranks;
      for (std::size_t r=0; r<nodes.size(); ++r) if (nodes[r]==n) ranks.push_back(r);
      for (int s : bundle) assign(p,s,leastLoaded(ranks,nodes,work),candidate,nodes,work);
    }
    auto best=original;
    if (admissible(candidate) && cost(p,candidate,nodes) < cost(p,best,nodes)) best=candidate;
    for (unsigned pass=0; pass<16; ++pass)
    {
      auto next=best;
      auto score=cost(p,best,nodes);
      for (std::size_t a=0; a<best.size(); ++a)
        for (std::size_t b=a+1; b<best.size(); ++b)
        {
          if (nodes[best[a]]==nodes[best[b]]) continue;
          auto trial=best;
          std::swap(trial[a],trial[b]);
          auto value=cost(p,trial,nodes);
          if (value < score && admissible(trial)) { next=std::move(trial); score=value; }
        }
      if (next==best) break;
      best=std::move(next);
    }
    return best;
  }

  struct Pools { std::vector<int> intra, extra; };

  inline Pools separatedPools(Problem const& p, std::vector<int> const& nodes,
                               std::string const& intraText, std::string const& extraText)
  {
    int count=*std::max_element(nodes.begin(),nodes.end())+1;
    if (count < 2) throw std::invalid_argument("Separated ownership requires >=2 physical nodes");
    Pools result{pool(intraText,count),pool(extraText,count)};
    auto ni=std::count(p.intracellular.begin(),p.intracellular.end(),true);
    auto ne=p.tags.size()-ni;
    if (!ni || !ne) throw std::invalid_argument("Separated ownership needs both I and E regions");
    if (result.intra.empty() && result.extra.empty())
    {
      std::uint64_t iw=0,total=std::accumulate(p.localDofs.begin(),p.localDofs.end(),0ULL);
      for (std::size_t s=0; s<p.tags.size(); ++s) if (p.intracellular[s]) iw+=p.localDofs[s];
      double error=std::numeric_limits<double>::infinity();
      int selected=-1;
      for (int k=1; k<count; ++k)
      {
        auto ranks=std::count_if(nodes.begin(),nodes.end(),[&](int n) { return n<k; });
        if (ranks > ni || nodes.size()-ranks > ne) continue;
        double e=std::abs(double(ranks)/nodes.size()-double(iw)/total);
        if (e<error) { selected=k; error=e; }
      }
      if (selected<0) throw std::invalid_argument("No fully occupied separated node pools fit region counts");
      for (int n=0; n<count; ++n) (n<selected ? result.intra : result.extra).push_back(n);
    }
    else if (result.intra.empty() || result.extra.empty())
    {
      auto const& present=result.intra.empty() ? result.extra : result.intra;
      auto& missing=result.intra.empty() ? result.intra : result.extra;
      for (int n=0; n<count; ++n)
        if (std::find(present.begin(),present.end(),n)==present.end()) missing.push_back(n);
    }
    std::set<int> covered;
    for (auto const* v : {&result.intra,&result.extra})
      for (int n : *v)
        if (!covered.insert(n).second) throw std::invalid_argument("I/E node pools overlap");
    if (result.intra.empty() || result.extra.empty() || covered.size()!=std::size_t(count))
      throw std::invalid_argument("I/E node pools must be nonempty and cover the allocation");
    return result;
  }

  inline void validateSeparated(Problem const& p, std::vector<int> const& owners,
                                std::vector<int> const& nodes, Pools const& pools)
  {
    validateOwners(p,owners,nodes);
    for (std::size_t s=0; s<owners.size(); ++s)
    {
      auto const& allowed=p.intracellular[s] ? pools.intra : pools.extra;
      if (std::find(allowed.begin(),allowed.end(),nodes[owners[s]])==allowed.end())
        throw std::invalid_argument("Ownership map violates the I/E physical-node pools");
    }
  }

  inline std::vector<int> separated(Problem const& p, std::vector<int> const& nodes,
                                    Pools const& pools)
  {
    std::vector<int> owners(p.tags.size(),-1);
    Loads work{std::vector<std::uint64_t>(nodes.size()),
               std::vector<std::uint64_t>(*std::max_element(nodes.begin(),nodes.end())+1)};
    for (bool intra : {false,true})
    {
      auto const& allowedNodes=intra ? pools.intra : pools.extra;
      std::vector<int> ranks;
      for (std::size_t r=0; r<nodes.size(); ++r)
        if (std::find(allowedNodes.begin(),allowedNodes.end(),nodes[r])!=allowedNodes.end())
          ranks.push_back(r);
      if (ranks.size()>std::size_t(std::count(p.intracellular.begin(),p.intracellular.end(),intra)))
        throw std::invalid_argument("Node pool has more MPI ranks than its material regions");
      for (int s : decreasingDofs(p))
        if (p.intracellular[s]==intra)
          assign(p,s,leastLoaded(ranks,nodes,work),owners,nodes,work);
    }
    validateSeparated(p,owners,nodes,pools);
    return owners;
  }

  inline std::vector<int> readOwners(std::string const& path, Problem const& p)
  {
    std::ifstream in(path);
    std::string line;
    if (!in || !std::getline(in,line))
      throw std::invalid_argument("Cannot read owner CSV");
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line!="subdomain_id,material_tag,owner_rank")
      throw std::invalid_argument("Owner CSV header must be subdomain_id,material_tag,owner_rank");
    std::vector<int> owners(p.tags.size(),-1);
    while (std::getline(in,line))
    {
      if (std::count(line.begin(),line.end(),',')!=2)
        throw std::invalid_argument("Owner CSV row must contain exactly three fields");
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream fields(line);
      int s,tag,rank;
      std::string extra;
      if (!(fields>>s>>tag>>rank) || (fields>>extra) || s<0 ||
          s>=static_cast<int>(p.tags.size()) || p.tags[s]!=tag || owners[s]!=-1 || rank<0)
        throw std::invalid_argument("Malformed, duplicate or mismatched owner CSV row");
      owners[s]=rank;
    }
    return owners;
  }

  struct Fingerprint
  {
    std::uint64_t value=14695981039346656037ULL;
    void add(std::uint64_t word)
    {
      for (unsigned i=0; i<8; ++i) { value^=(word>>(8*i))&255; value*=1099511628211ULL; }
    }
    std::string str() const
    { std::ostringstream out; out<<std::hex<<std::setw(16)<<std::setfill('0')<<value; return out.str(); }
  };

  struct State
  {
    std::string policy;
    Problem problem;
    std::vector<int> leaders, nodes, owners;
    Pools pools;
    double setupSeconds=0;
    std::string problemFingerprint;
  };

  /// Query through one helper so older test/options fixtures keep the default API.
  template<class Options>
  Kaskade::BDDC::MpiOwnerMap ownerMap(Options const& options)
  {
    if constexpr (requires { options.bddcOwnership; }) return options.bddcOwnership;
    return {};
  }

  /** Collectively prepare one immutable map BEFORE material/interface assembly.
   * P2 without reporting/map input returns immediately, with no new collectives.
   * Exceptions are agreed before any rank enters the owner-specific assembly.
   */
  template<class Options,class Data>
  void prepare(Options& options, Data const& data,
               std::vector<int> const& extraTags, std::vector<int> const& intraTags)
  {
    auto policy=mode(options.bddcPartitionMode);
    if (policy=="baseline" && !options.bddcPartitionReport && options.bddcPartitionMap.empty())
      return;
#ifdef KASKADE_HAVE_MPI
    auto start=std::chrono::steady_clock::now();
    int rank=0,ranks=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    Kaskade::BDDC::MpiNodeTopology topology;
    topology.configure();
    auto state=std::make_shared<State>();
    state->policy=policy; state->leaders=topology.leaders(); state->nodes=nodeIds(state->leaders);
    std::string error;
    try
    {
      auto& p=state->problem;
      p.tags=data.tags;
      for (std::size_t s=0; s<data.tags.size(); ++s)
      {
        bool i=std::find(intraTags.begin(),intraTags.end(),data.tags[s])!=intraTags.end();
        bool e=std::find(extraTags.begin(),extraTags.end(),data.tags[s])!=extraTags.end();
        if (i==e) throw std::invalid_argument("Partition material tag must have exactly one I/E class");
        p.intracellular.push_back(i); p.localDofs.push_back(data.localDofs[s].size());
      }
      std::map<std::pair<int,int>,std::uint64_t> edges;
      for (auto const& copies : data.sharedDofs)
        for (std::size_t a=0; a<copies.size(); ++a)
          for (std::size_t b=a+1; b<copies.size(); ++b)
          {
            int x=copies[a].s,y=copies[b].s;
            if (x>y) std::swap(x,y);
            if (x!=y) ++edges[{x,y}];
          }
      for (auto const& [ends,w] : edges) p.edges.push_back({ends.first,ends.second,w});
      validateProblem(p);
      state->owners=baseline(p.tags.size(),ranks);
      if (policy=="separated")
        state->pools=separatedPools(p,state->nodes,options.bddcIntracellularNodes,
                                   options.bddcExtracellularNodes);
      if (!options.bddcPartitionMap.empty()) state->owners=readOwners(options.bddcPartitionMap,p);
      else if (policy=="coupled") state->owners=coupled(p,state->nodes,options.bddcPartitionBalanceFactor);
      else if (policy=="separated") state->owners=separated(p,state->nodes,state->pools);
      validateOwners(p,state->owners,state->nodes);
      if (!options.bddcPartitionMap.empty() && policy!="baseline")
      {
        auto work=loads(p,state->owners,state->nodes);
        if (std::find(work.rank.begin(),work.rank.end(),0)!=work.rank.end())
          throw std::invalid_argument("Controlled owner CSV must keep every allocated rank active");
        if (policy=="coupled")
        {
          auto original=baseline(p.tags.size(),ranks);
          auto base=loads(p,original,state->nodes);
          if (!std::isfinite(options.bddcPartitionBalanceFactor) || options.bddcPartitionBalanceFactor<1 ||
              *std::max_element(work.rank.begin(),work.rank.end())>
                options.bddcPartitionBalanceFactor*(*std::max_element(base.rank.begin(),base.rank.end())) ||
              *std::max_element(work.node.begin(),work.node.end())>
                options.bddcPartitionBalanceFactor*(*std::max_element(base.node.begin(),base.node.end())) ||
              cost(p,state->owners,state->nodes)>cost(p,original,state->nodes))
            throw std::invalid_argument("Coupled owner CSV exceeds the baseline cut/load bounds");
        }
      }
      if (policy=="baseline" && state->owners!=baseline(p.tags.size(),ranks))
        throw std::invalid_argument("Baseline partition cannot override cyclic ownership");
      if (policy=="separated") validateSeparated(p,state->owners,state->nodes,state->pools);
      Fingerprint fp;
      for (std::size_t s=0; s<p.tags.size(); ++s)
      {
        fp.add(p.tags[s]); fp.add(p.intracellular[s]); fp.add(p.localDofs[s]);
        for (auto i : data.localDofs[s]) fp.add(i);
      }
      for (auto const& e : p.edges) { fp.add(e.a); fp.add(e.b); fp.add(e.dofs); }
      state->problemFingerprint=fp.str();
    }
    catch (std::exception const& e) { error=e.what(); }
    int bad=!error.empty(),any=0;
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (any) throw std::runtime_error("Invalid collective EMI ownership setup: "+
                                      (error.empty() ? std::string("another rank rejected the map") : error));
    int count=state->owners.size(), minCount=0, maxCount=0;
    MPI_Allreduce(&count,&minCount,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&count,&maxCount,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (minCount!=maxCount) throw std::runtime_error("Partition metadata lengths differ across ranks");
    // Exact vector agreement, not a collision-prone hash comparison.
    auto rootOwners=state->owners;
    MPI_Bcast(rootOwners.data(),rootOwners.size(),MPI_INT,0,MPI_COMM_WORLD);
    bad=rootOwners!=state->owners;
    auto fp=state->problemFingerprint;
    std::string rootFp=fp;
    MPI_Bcast(rootFp.data(),rootFp.size(),MPI_CHAR,0,MPI_COMM_WORLD);
    bad=bad || rootFp!=fp;
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (any) throw std::runtime_error("EMI ownership/metadata differ across MPI ranks");
    state->setupSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    options.bddcPartitionState=state;
    if (policy!="baseline")
      options.bddcOwnership=std::make_shared<std::vector<int> const>(state->owners);
#else
    (void)options; (void)data; (void)extraTags; (void)intraTags;
    throw std::invalid_argument("Controlled ownership requires an MPI build");
#endif
  }

  inline std::string quoted(std::string s)
  {
    std::string out="\"";
    for (char c : s) { if (c=='\"') out+='\"'; out+=(c=='\n' || c=='\r') ? ' ' : c; }
    return out+'\"';
  }

  /** Export actual placement and full local-matrix work proxies once, outside
   * solver iterations. There is no anatomical pairing inferred from tag numbers.
   */
  template<class Options,class Data,class Matrix>
  void writeMetadata(Options const& options, Data const& data, Matrix const& global)
  {
    if (!options.bddcPartitionReport || !options.bddcPartitionState) return;
#ifdef KASKADE_HAVE_MPI
    auto const& state=*options.bddcPartitionState;
    int rank=0,ranks=state.nodes.size();
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    char host[MPI_MAX_PROCESSOR_NAME]={}; int len=0;
    MPI_Get_processor_name(host,&len);
    std::vector<char> hosts(ranks*MPI_MAX_PROCESSOR_NAME);
    MPI_Allgather(host,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,hosts.data(),MPI_MAX_PROCESSOR_NAME,
                  MPI_CHAR,MPI_COMM_WORLD);
    std::string affinity="unavailable";
#ifdef __linux__
    cpu_set_t mask;
    if (sched_getaffinity(0,sizeof(mask),&mask)==0)
    {
      affinity.clear();
      for (int i=0; i<CPU_SETSIZE; ++i)
        if (CPU_ISSET(i,&mask)) affinity+=(affinity.empty() ? "" : ";")+std::to_string(i);
    }
#endif
    char binding[4096]={};
    if (affinity.size()>=sizeof(binding)) affinity="affinity-list-too-large";
    std::copy(affinity.begin(),affinity.end(),binding);
    std::vector<char> bindings(ranks*sizeof(binding));
    MPI_Allgather(binding,sizeof(binding),MPI_CHAR,bindings.data(),sizeof(binding),MPI_CHAR,MPI_COMM_WORLD);
    int failed=0;
    if (rank==0)
    {
      try
      {
        std::filesystem::create_directories(options.outputDir);
        auto open=[&](char const* name)
        {
          std::ofstream out(std::filesystem::path(options.outputDir)/name);
          out.exceptions(std::ios::failbit|std::ios::badbit);
          out<<std::setprecision(17); return out;
        };
        auto placement=open("partition-ranks.csv");
        placement<<"rank,node_group,shared_leader,hostname,cpu_affinity\n";
        for (int r=0; r<ranks; ++r)
          placement<<r<<','<<state.nodes[r]<<','<<state.leaders[r]<<','
                   <<quoted(hosts.data()+r*MPI_MAX_PROCESSOR_NAME)<<','
                   <<quoted(bindings.data()+r*sizeof(binding))<<'\n';
        char library[MPI_MAX_LIBRARY_VERSION_STRING]={}; int libraryLength=0;
        MPI_Get_library_version(library,&libraryLength);
        auto environment=open("partition-mpi.csv");
        environment<<"mpi_library\n"<<quoted(std::string(library,libraryLength))<<'\n';
        auto owners=open("partition-ownership.csv");
        owners<<"subdomain_id,material_tag,region,owner_rank,node_group,local_dofs,local_matrix_nnz\n";
        Fingerprint matrices,weights,allGlobal;
        std::vector<std::uint64_t> nnz(ranks),rankDofs(ranks);
        auto matrixHash=[](auto const& matrix,Fingerprint& hash)
        {
          std::uint64_t count=0; hash.add(matrix.N()); hash.add(matrix.M());
          for (auto row=matrix.begin(); row!=matrix.end(); ++row)
            for (auto col=row.begin(); col!=row.end(); ++col)
            {
              hash.add(row.index()); hash.add(col.index());
              hash.add(std::bit_cast<std::uint64_t>(double((*col)[0][0]))); ++count;
            }
          return count;
        };
        for (std::size_t s=0; s<data.tags.size(); ++s)
        {
          auto n=matrixHash(data.localMatrices[s],matrices);
          int r=state.owners[s]; nnz[r]+=n; rankDofs[r]+=data.localDofs[s].size();
          owners<<s<<','<<data.tags[s]<<','<<(state.problem.intracellular[s] ? "I" : "E")
                <<','<<r<<','<<state.nodes[r]<<','<<data.localDofs[s].size()<<','<<n<<'\n';
          for (auto const& w : data.weights[s])
            weights.add(std::bit_cast<std::uint64_t>(double(w[0])));
        }
        matrixHash(global,allGlobal);
        auto problem=open("partition-problem.csv");
        problem<<"mode,metadata_fingerprint,global_matrix_fingerprint,local_matrices_fingerprint,weights_fingerprint,planning_s\n"
               <<state.policy<<','<<state.problemFingerprint<<','<<allGlobal.str()<<','
               <<matrices.str()<<','<<weights.str()<<','<<state.setupSeconds<<'\n';
        auto work=open("partition-work.csv");
        work<<"scope,id,local_dofs,local_matrix_nnz,owned_subdomains\n";
        auto totals=loads(state.problem,state.owners,state.nodes);
        for (int r=0; r<ranks; ++r)
          work<<"rank,"<<r<<','<<rankDofs[r]<<','<<nnz[r]<<','
              <<std::count(state.owners.begin(),state.owners.end(),r)<<'\n';
        for (std::size_t n=0; n<totals.node.size(); ++n)
        {
          std::uint64_t nNnz=0, owned=0;
          for (int r=0; r<ranks; ++r) if (state.nodes[r]==int(n))
          { nNnz+=nnz[r]; owned+=std::count(state.owners.begin(),state.owners.end(),r); }
          work<<"node,"<<n<<','<<totals.node[n]<<','<<nNnz<<','<<owned<<'\n';
        }
        auto edges=open("partition-interfaces.csv");
        edges<<"subdomain_a,subdomain_b,interface_dofs,rank_a,rank_b,node_a,node_b,route,I_E\n";
        for (auto const& e : state.problem.edges)
        {
          int a=state.owners[e.a],b=state.owners[e.b];
          int route=Kaskade::BDDC::mpiRouteClass(a,b,state.leaders);
          edges<<e.a<<','<<e.b<<','<<e.dofs<<','<<a<<','<<b<<','<<state.nodes[a]<<','
               <<state.nodes[b]<<','<<route<<','
               <<(state.problem.intracellular[e.a]!=state.problem.intracellular[e.b])<<'\n';
        }
      }
      catch (std::exception const&) { failed=1; }
    }
    MPI_Bcast(&failed,1,MPI_INT,0,MPI_COMM_WORLD);
    if (failed) throw std::runtime_error("Cannot write collective partition metadata");
#else
    (void)data; (void)global;
#endif
  }

  /// Export the actual unchanged coarse rows, outside the solve and on rank 0.
  template<class Options,class Interface>
  void writeCoarseMetadata(Options const& options, Interface const& interface,
                          std::size_t subdomains)
  {
    if constexpr (requires { options.bddcPartitionReport; })
    {
      if (!options.bddcPartitionReport) return;
#ifdef KASKADE_HAVE_MPI
      int rank=0,failed=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (!rank)
      {
        try
        {
          Fingerprint hash;
          auto constraints=interface.coarseConstraints();
          hash.add(constraints.size());
          for (auto const& [owners,count] : constraints)
          {
            hash.add(owners.size()); for (int owner : owners) hash.add(owner);
            hash.add(count);
          }
          for (std::size_t s=0; s<subdomains; ++s)
          {
            auto matrix=interface.template coarseConstraint<double>(s);
            hash.add(matrix.N()); hash.add(matrix.M());
            for (auto row=matrix.begin(); row!=matrix.end(); ++row)
              for (auto col=row->begin(); col!=row->end(); ++col)
              {
                hash.add(row.index()); hash.add(col.index());
                hash.add(std::bit_cast<std::uint64_t>(double((*col)[0][0])));
              }
          }
          std::ofstream out(std::filesystem::path(options.outputDir)/"partition-coarse.csv");
          out.exceptions(std::ios::failbit|std::ios::badbit);
          out<<"constraint_groups,coarse_rows_fingerprint\n"<<constraints.size()<<','<<hash.str()<<'\n';
          out.close();
        }
        catch (std::exception const&) { failed=1; }
      }
      MPI_Bcast(&failed,1,MPI_INT,0,MPI_COMM_WORLD);
      if (failed) throw std::runtime_error("Cannot write partition coarse metadata");
#else
      (void)interface; (void)subdomains;
#endif
    }
  }
}
#endif
