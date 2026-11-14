#ifndef EMI_BDDC_TIME_STUDY_HH
#define EMI_BDDC_TIME_STUDY_HH

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "mg/bddcTransferTiming.hh"
#include "mg/bddcLiftingPerformance.hh"

namespace EmiBddc
{
/** Fixed-topology FCG observers. No I/O or new collectives in a timestep.
 * Rank-local step walls include RHS assembly, restart and true-residual stopping;
 * solve-call walls exclude caller residual checks and state snapshots. Timings
 * are nested: per-call and per-step profiles must never be added together.
 * Detailed scalar/level observers perturb work; mode 0 is the speed endpoint.
 */
class TimeStudy
{
  using Clock=std::chrono::steady_clock;
  using Counters=Kaskade::BDDC::TransferTiming::Counters;
  using Packet=std::tuple<std::uint64_t,std::uint64_t,std::array<double,5>,std::uint64_t>;
  struct Row
  {
    int step=0, calls=0;
    double wall=0, fcg=0, relative=0;
    Packet packet{};
  };
  struct StageRow
  {
    int step, call;
    unsigned phase, stage, level;
    Kaskade::BDDC::TransferTiming::Sample sample;
  };
  inline static Clock::time_point applicationStart{};
  bool enabled;
  double applicationSetup=0;
  std::uint64_t ownedSubdomains=0, ownedDofs=0, globalDofs=0;
  Clock::time_point setupStart{}, stepStart{}, callStart{};
  double fcgSeconds=0;
  unsigned call=0;
  Packet previousPacket{};
  Counters stepProfile{}, previousCall{};
  std::vector<Row> rows;
  std::vector<StageRow> stages;

  static double elapsed(Clock::time_point start)
  { return std::chrono::duration<double>(Clock::now()-start).count(); }

  void profile(int step, int index, Counters const& data)
  {
    for (unsigned p=0;p<3;++p)
      for (unsigned s=0;s<static_cast<unsigned>(Kaskade::BDDC::TransferTiming::Stage::Count);++s)
        for (unsigned l=0;l<65;++l)
          if (data[p][s][l].calls) stages.push_back({step,index,p,s,l,data[p][s][l]});
  }

  static Packet difference(Packet const& now, Packet const& before)
  {
    auto times=std::get<2>(now);
    for (unsigned i=0;i<times.size();++i) times[i]-=std::get<2>(before)[i];
    return {std::get<0>(now)-std::get<0>(before),
            std::get<1>(now)-std::get<1>(before),times,std::get<3>(now)-std::get<3>(before)};
  }

public:
  static void beginApplication() { applicationStart=Clock::now(); }

  explicit TimeStudy(bool observe):enabled(observe)
  {
    if (!enabled) return;
    Kaskade::BDDC::TransferTiming::reset();
    setupStart=Clock::now();
  }

  void setWork(std::uint64_t subdomains, std::uint64_t dofs, std::uint64_t global)
  { ownedSubdomains=subdomains; ownedDofs=dofs; globalDofs=global; }

  void setupComplete(Packet const& packet)
  {
    if (!enabled) return;
    applicationSetup=elapsed(applicationStart==Clock::time_point{} ? setupStart : applicationStart);
    rows.push_back({0,0,elapsed(setupStart),0,0,packet});
    if (Kaskade::BDDC::TransferTiming::mode())
      profile(0,0,Kaskade::BDDC::TransferTiming::snapshot());
    previousPacket=packet;
  }

  void beginStep()
  {
    if (!enabled) return;
    if (Kaskade::BDDC::TransferTiming::mode())
      stepProfile=Kaskade::BDDC::TransferTiming::snapshot();
    fcgSeconds=0; call=0;
    stepStart=Clock::now();
  }

  void beginCall()
  {
    if (!enabled) return;
    // Exclude RHS restart and outer residual checks from the call-only profile.
    if (Kaskade::BDDC::TransferTiming::mode() && call<2)
      previousCall=Kaskade::BDDC::TransferTiming::snapshot();
    callStart=Clock::now();
  }

  void endCall(int step)
  {
    if (!enabled) return;
    fcgSeconds+=elapsed(callStart);
    ++call;
    if (Kaskade::BDDC::TransferTiming::mode() && call<=2)
    {
      auto now=Kaskade::BDDC::TransferTiming::snapshot();
      profile(step,call,Kaskade::BDDC::TransferTiming::difference(now,previousCall));
      previousCall=now;
    }
  }

  void endStep(int step, int calls, double relative, Packet const& packet)
  {
    if (!enabled) return;
    double wall=elapsed(stepStart);
    if (calls!=static_cast<int>(call)) throw std::runtime_error("Incomplete FCG timing calls");
    rows.push_back({step,calls,wall,fcgSeconds,relative,difference(packet,previousPacket)});
    previousPacket=packet;
    if (Kaskade::BDDC::TransferTiming::mode())
      profile(step,0,Kaskade::BDDC::TransferTiming::difference(
          Kaskade::BDDC::TransferTiming::snapshot(),stepProfile));
  }

  void write(std::string const& directory, int rank) const
  {
    if (!enabled) return;
    double simulation=elapsed(setupStart);
    double application=elapsed(applicationStart==Clock::time_point{} ? setupStart : applicationStart);
    std::filesystem::create_directories(directory);
    std::string stem=directory+"/hb-time-rank-"+std::to_string(rank);
    std::ofstream total(stem+"-total.csv");
    total << "simulation_s,application_compute_s,application_setup_s,owned_subdomains,owned_local_dofs,global_dofs,process_peak_rss_bytes\n"
          << std::setprecision(17) << simulation << ',' << application << ',' << applicationSetup << ',' << ownedSubdomains
          << ',' << ownedDofs << ',' << globalDofs << ','
          << Kaskade::BDDC::LiftingPerformance::peakRss() << '\n';
    total.close();
    if (!total) throw std::runtime_error("Cannot write simulation timing "+stem);
    std::ofstream out(stem+".csv");
    out << "step,calls,wall_s,fcg_s,relative_residual,original_bytes,packed_bytes,"
           "pack_s,unpack_s,exchange_s,verify_s,packet_solve_s,inter_node_packed_bytes\n" << std::setprecision(17);
    for (auto const& row:rows)
    {
      out << row.step << ',' << row.calls << ',' << row.wall << ',' << row.fcg << ','
          << row.relative << ',' << std::get<0>(row.packet) << ',' << std::get<1>(row.packet);
      for (double value:std::get<2>(row.packet)) out << ',' << value;
      out << ',' << std::get<3>(row.packet) << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot write timestep timings "+stem);
    std::ofstream detail(stem+"-stages.csv");
    detail << "step,call,phase,stage,level,calls,seconds\n" << std::setprecision(17);
    char const* phases[]={"restrict","prolongate","other"};
    for (auto const& row:stages)
      detail << row.step << ',' << row.call << ',' << phases[row.phase] << ','
             << Kaskade::BDDC::TransferTiming::name(
                  static_cast<Kaskade::BDDC::TransferTiming::Stage>(row.stage)) << ','
             << row.level << ',' << row.sample.calls << ',' << row.sample.seconds << '\n';
    detail.close();
    if (!detail) throw std::runtime_error("Cannot write timestep stage timings "+stem);
  }
};
}
#endif
