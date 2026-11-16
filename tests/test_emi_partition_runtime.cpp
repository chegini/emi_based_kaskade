#include "EMI_partitioning.hh"

#include <iostream>

using namespace EmiBddc::Partitioning;

struct Block { double v; double const* operator[](int) const { return &v; } };
struct SmallMatrix
{
  struct Row
  {
    std::array<Block,2> a;
    struct It
    {
      Row const* row; int i;
      int index() const { return i; }
      Block const& operator*() const { return row->a[i]; }
      It& operator++() { ++i; return *this; }
      bool operator!=(It const& b) const { return i!=b.i; }
    };
    It begin() const { return {this,0}; }
    It end() const { return {this,2}; }
  };
  std::array<Row,2> rows;
  explicit SmallMatrix(double scale=1):rows{{Row{{Block{4*scale},Block{scale}}},
                                           Row{{Block{scale},Block{3*scale}}}}} {}
  int N() const { return 2; }
  int M() const { return 2; }
  struct It
  {
    SmallMatrix const* matrix; int i;
    int index() const { return i; }
    Row const* operator->() const { return &matrix->rows[i]; }
    Row::It begin() const { return matrix->rows[i].begin(); }
    Row::It end() const { return matrix->rows[i].end(); }
    It& operator++() { ++i; return *this; }
    bool operator!=(It const& b) const { return i!=b.i; }
  };
  It begin() const { return {this,0}; }
  It end() const { return {this,2}; }
};

struct Data
{
  struct Dof { int s, i; };
  std::vector<int> tags;
  std::vector<std::vector<std::size_t>> localDofs;
  std::vector<std::vector<Dof>> sharedDofs;
  std::vector<SmallMatrix> localMatrices;
  std::vector<std::vector<std::array<double,1>>> weights;
};

struct Options
{
  std::string bddcPartitionMode="baseline", bddcPartitionMap;
  std::string bddcIntracellularNodes, bddcExtracellularNodes, outputDir;
  double bddcPartitionBalanceFactor=1;
  int bddcPartitionReport=0;
  Kaskade::BDDC::MpiOwnerMap bddcOwnership;
  std::shared_ptr<State const> bddcPartitionState;
};

struct Constraints
{
  auto coarseConstraints() const { return std::vector<std::pair<std::vector<int>,int>>{{{0,1},1}}; }
  template<class T> SmallMatrix coarseConstraint(int) const { return SmallMatrix{}; }
};

void require(bool ok, char const* reason)
{ if (!ok) throw std::runtime_error(reason); }

template<class F>
void rejectedEverywhere(F&& f, int ranks)
{
  int bad=0,total=0;
  try { f(); } catch (std::exception const&) { bad=1; }
  MPI_Allreduce(&bad,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  require(total==ranks,"Invalid collective ownership did not fail on every rank");
}

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);
  int rank=0,ranks=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try
  {
    require(argc==2 || argc==3,"Expected temporary output directory [physical-node-count]");
    int n=std::max(8,2*ranks);
    Data data; data.sharedDofs.resize(2);
    std::vector<int> intra,extra;
    for (int s=0; s<n; ++s)
    {
      data.tags.push_back(100+s);
      (s%2 ? intra : extra).push_back(100+s);
      data.localDofs.push_back({0,1});
      data.localMatrices.emplace_back(1.0/n);
      data.weights.push_back({{1.0/n},{1.0/n}});
      for (int i=0; i<2; ++i) data.sharedDofs[i].push_back({s,i});
    }
    Options defaultOptions;
    prepare(defaultOptions,data,extra,intra);
    require(!defaultOptions.bddcOwnership && !defaultOptions.bddcPartitionState,
            "P2 default unexpectedly enabled a new setup path");
    Kaskade::BDDC::MpiNodeTopology topology; topology.configure();
    if (argc==3) require(topology.nodes()==std::stoi(argv[2]),"Physical/shared node count differs from requested preflight");
    for (std::string policy : {"baseline","coupled","separated"})
    {
      Options options; options.bddcPartitionMode=policy; options.bddcPartitionReport=1;
      options.outputDir=std::string(argv[1])+"/"+policy;
      if (policy=="separated" && topology.nodes()==1)
      {
        rejectedEverywhere([&] { prepare(options,data,extra,intra); },ranks);
        continue; // No physical separation is claimed on one local machine.
      }
      prepare(options,data,extra,intra);
      require(bool(options.bddcPartitionState),"Missing partition metadata state");
      require((policy=="baseline")==!options.bddcOwnership,"P2 did not keep the legacy map path");
      auto const& saved=options.bddcPartitionState->owners;
      Kaskade::BDDC::MpiDomainLayout layout;
      layout.configure(true,n,false,options.bddcOwnership);
      for (int s=0; s<n; ++s) require(layout.owner(s)==saved[s],"Assembly/solver map mismatch");
      writeMetadata(options,data,SmallMatrix{}); // Instantiate real export against small algebra.
      writeCoarseMetadata(options,Constraints{},n);
      // A small SPD assembly/solution test; this is not an EMI or FCG benchmark.
      for (int step=0; step<3; ++step)
      {
        std::array<double,4> a{};
        for (int s=0; s<n; ++s) if (layout.owns(s))
        {
          require(layout.owner(s)==saved[s],"Owner changed across timesteps");
          a[0]+=4.0/n; a[1]+=1.0/n; a[2]+=1.0/n; a[3]+=3.0/n;
        }
        MPI_Allreduce(MPI_IN_PLACE,a.data(),4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        require(std::abs(a[0]-4)<1e-14 && std::abs(a[1]-1)<1e-14 &&
                std::abs(a[2]-1)<1e-14 && std::abs(a[3]-3)<1e-14,
                "Ownership changed the underlying SPD operator");
        double det=a[0]*a[3]-a[1]*a[2];
        double x=(a[3]-2*a[1])/det, y=(2*a[0]-a[2])/det;
        double r0=1-a[0]*x-a[1]*y,r1=2-a[2]*x-a[3]*y;
        require(std::hypot(r0,r1)/std::sqrt(5.0)<1e-14,"True residual changed under ownership");
        require(std::abs(x-1.0/11)<1e-14 && std::abs(y-7.0/11)<1e-14,
                "Small SPD solution changed under ownership");
      }
    }
    // Core layout rejects malformed explicit maps before owned storage is used.
    for (auto invalid : {std::vector<int>{0}, std::vector<int>(n,ranks), std::vector<int>(n,-1)})
      rejectedEverywhere([&]
      {
        Kaskade::BDDC::MpiDomainLayout layout;
        layout.configure(true,n,false,std::make_shared<std::vector<int> const>(invalid));
      },ranks);
    auto mismatch=data;
    if (rank==ranks-1) mismatch.localDofs[0].push_back(2);
    if (ranks>1)
      rejectedEverywhere([&]
      {
        Options o;o.bddcPartitionMode="coupled";prepare(o,mismatch,extra,intra);
      },ranks);
    if (rank==0)
      std::cout<<"PASS: actual collective preparation/export, P2 default path, immutable owners, "
                   "three-step SPD operator/solution/residual equivalence; physical node groups="
               <<topology.nodes()<<"; no EMI trajectory executed\n";
    MPI_Finalize();return 0;
  }
  catch (std::exception const& e)
  { std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1; }
}
