#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../EMI_phase_profile.hh"

struct FakeClock
{
  using duration = std::chrono::seconds;
  using time_point = std::chrono::time_point<FakeClock,duration>;
  inline static duration elapsed{};
  inline static int calls = 0;
  static time_point now() { ++calls; return time_point(elapsed); }
};

int main(int argc, char** argv)
{
  int rank = 0, ranks = 1;
  bool mpi = false;
#ifdef KASKADE_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  mpi = true;
#else
  (void)argc;
  (void)argv;
#endif
  using Profile = EmiBddc::BasicPhaseProfile<FakeClock>;
  using Phase = EmiBddc::EmiPhase;
  std::ostringstream disabledOutput;
  Profile disabled(false,mpi);
  disabled.next(Phase::meshLoading);
  disabled.printGlobal(disabledOutput);
  assert(FakeClock::calls == 0 && disabledOutput.str().empty());

  // Deterministic rank-dependent intervals verify exclusive accounting,
  // accumulating repeated phases, and the MPI min/mean/max reductions.
  int const scale = rank+1;
  Profile profile(true,mpi);
  profile.next(Phase::threadPool);
  FakeClock::elapsed += FakeClock::duration(2*scale);
  profile.next(Phase::meshLoading);
  FakeClock::elapsed += FakeClock::duration(5*scale);
  profile.next(Phase::other);
  FakeClock::elapsed += FakeClock::duration(3*scale);
  profile.next(Phase::meshLoading);
  FakeClock::elapsed += FakeClock::duration(3*scale);
  std::ostringstream output;
  profile.printGlobal(output);
  assert(FakeClock::calls == 6);
  profile.next(Phase::solver);
  std::ostringstream duplicate;
  profile.printGlobal(duplicate);
  assert(duplicate.str().empty() && FakeClock::calls == 6);
  if (rank == 0)
  {
    std::vector<std::pair<std::string,double>> expected = {
      {"thread-pool initialization",2.0},
      {"mesh loading and grid creation",8.0},
      {"other setup/bookkeeping",3.0},
      {"accounted total",13.0}
    };
    for (auto const& [label,value] : expected)
    {
      auto const pos = output.str().find(label);
      assert(pos != std::string::npos);
      std::istringstream row(output.str().substr(pos+label.size()));
      double minimum, mean, maximum;
      row >> minimum >> mean >> maximum;
      assert(row && std::abs(minimum-value) < 1e-6);
      assert(std::abs(mean-value*(ranks+1)/2.0) < 1e-6);
      assert(std::abs(maximum-value*ranks) < 1e-6);
    }
    std::cout << "EMI phase profiling checks passed (" << ranks << " rank(s))\n";
  }
  else
    assert(output.str().empty());
#ifdef KASKADE_HAVE_MPI
  MPI_Finalize();
#endif
}
