#ifndef EMI_PHASE_PROFILE_HH
#define EMI_PHASE_PROFILE_HH

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ostream>
#include <sstream>

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

namespace EmiBddc
{
  enum class EmiPhase
  {
    threadPool, meshLoading, materialAndRefinement, spacesAndInitialState,
    initialVtk, globalAssembly, partitionMetadata, materialAssembly,
    materialExchange, interfaceDiscovery, interfaceAssembly, interfaceRedistribution,
    localMatrixExtraction, diagonalWeights, diagnostics, solver, vtkOutput,
    other, count
  };

  // Mutually exclusive local intervals, with no synchronization at phase
  // transitions. Waiting is charged to the phase in which it occurs. Only
  // final reporting reduces timings; column maxima need not occur on the same
  // rank and must not be added to estimate total elapsed time.
  template <class Clock>
  class BasicPhaseProfile
  {
  public:
    BasicPhaseProfile(bool enabled, bool mpiEnabled)
      : enabled_(enabled), mpiEnabled_(mpiEnabled),
        start_(enabled ? Clock::now() : typename Clock::time_point{})
    {}

    void next(EmiPhase phase)
    {
      if (!enabled_ || finished_)
        return;
      auto const now = Clock::now();
      times_[static_cast<size_t>(phase_)] +=
        std::chrono::duration<double>(now-start_).count();
      phase_ = phase;
      start_ = now;
    }

    void printGlobal(std::ostream& output = std::cout, bool show = true)
    {
      if (!enabled_ || finished_)
        return;
      next(EmiPhase::other);
      finished_ = true;
      std::array<double,sampleCount> local{};
      for (size_t i = 0; i < phaseCount; ++i)
        local[i] = times_[i];
      local[phaseCount] = std::accumulate(times_.begin(),times_.end(),0.0);
      auto minimum = local, sum = local, maximum = local;
      int rank = 0, ranks = 1;
#ifdef KASKADE_HAVE_MPI
      if (mpiEnabled_)
      {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (initialized)
        {
          MPI_Comm_rank(MPI_COMM_WORLD,&rank);
          MPI_Comm_size(MPI_COMM_WORLD,&ranks);
          MPI_Reduce(local.data(),minimum.data(),sampleCount,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
          MPI_Reduce(local.data(),sum.data(),sampleCount,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(local.data(),maximum.data(),sampleCount,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        }
      }
#else
      (void)mpiEnabled_;
#endif
      if (rank != 0 || !show)
        return;
      std::ostringstream report;
      report << "Profile (wall seconds, EMI end-to-end phases; " << ranks << " rank(s)):\n"
             << std::left << std::setw(44) << "  phase"
             << std::right << std::setw(14) << "min"
             << std::setw(14) << "mean" << std::setw(14) << "max" << '\n'
             << std::fixed << std::setprecision(6);
      for (size_t i = 0; i < local.size(); ++i)
        report << std::left << std::setw(44) << labels_[i]
               << std::right << std::setw(14) << minimum[i]
               << std::setw(14) << sum[i]/ranks << std::setw(14) << maximum[i] << '\n';
      report << "  Local phase intervals are exclusive; waiting is included.\n"
             << "  Phase maxima are not additive; use the accounted-total row.\n"
             << "  Solver wrapper includes local/coarse setup and existing solver profiles.\n";
      output << report.str();
    }

  private:
    static constexpr size_t phaseCount = static_cast<size_t>(EmiPhase::count);
    static constexpr int sampleCount = static_cast<int>(phaseCount)+1;
    inline static constexpr std::array<char const*,sampleCount> labels_ = {
      "  thread-pool initialization", "  mesh loading and grid creation",
      "  material tags and refinement", "  spaces, initial state and DOF support",
      "  initial VTK output", "  global system assembly",
      "  BDDC partition metadata", "  material subdomain assembly",
      "  material matrix exchange",
      "  interface-neighbor discovery", "  interface mass-block assembly",
      "  interface mass redistribution", "  local matrix extraction",
      "  diagonal weights", "  diagnostics and interface splitting",
      "  solver wrapper (including setup)", "  timestep/final VTK output",
      "  other setup/bookkeeping", "  accounted total"
    };
    bool enabled_, mpiEnabled_, finished_ = false;
    EmiPhase phase_ = EmiPhase::other;
    typename Clock::time_point start_;
    std::array<double,phaseCount> times_{};
  };

  using PhaseProfile = BasicPhaseProfile<std::chrono::steady_clock>;
}

#endif
