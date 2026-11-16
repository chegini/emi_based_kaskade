#ifndef EMI_PARTITION_REPORT_HH
#define EMI_PARTITION_REPORT_HH

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "mg/bddcMpiRouting.hh"
#include "mg/bddcTransferTiming.hh"

namespace EmiBddc::Partitioning
{
  /** Buffer observations until after the retained solve, with no per-exchange
   * I/O/collectives. Mixed Alltoallv elapsed time has no network-only attribution.
   * Parent codec/phase and child kernel timers are nested, not additive costs.
   */
  class TransportReport
  {
    using Clock=std::chrono::steady_clock;
    using Counters=Kaskade::BDDC::TransferTiming::Counters;
    using Stage=Kaskade::BDDC::TransferTiming::Stage;
    struct Row
    {
      int step, calls, accepted, timingMode;
      std::string status;
      double wall, fcg;
      Kaskade::BDDC::MpiTransportTraffic traffic;
      Counters counters;
    };
    bool enabled_;
    Clock::time_point start;
    Counters before{};
    std::vector<Row> rows;

    static double seconds(Counters const& c, unsigned phase, Stage s)
    { return c[phase][static_cast<unsigned>(s)][0].seconds; }

    static Kaskade::BDDC::MpiMessageBucket combine(
        Kaskade::BDDC::MpiMessageBucket a, Kaskade::BDDC::MpiMessageBucket const& b)
    {
      a.rawReferenceBytes+=b.rawReferenceBytes; a.prepackingBytes+=b.prepackingBytes;
      a.sentBytes+=b.sentBytes; a.messages+=b.messages;
      a.maxRawReferenceBytes=std::max(a.maxRawReferenceBytes,b.maxRawReferenceBytes);
      a.maxSentBytes=std::max(a.maxSentBytes,b.maxSentBytes);
      a.rawReferenceKnown=a.rawReferenceKnown && b.rawReferenceKnown;
      return a;
    }

  public:
    explicit TransportReport(bool enabled):enabled_(enabled) {}
    bool enabled() const { return enabled_; }

    void beginStep()
    {
      if (!enabled_) return;
      start=Clock::now();
      if (Kaskade::BDDC::TransferTiming::mode())
        before=Kaskade::BDDC::TransferTiming::snapshot();
    }

    void endStep(int step, int calls, int accepted, std::string status,
                 Kaskade::BDDC::MpiTransportTraffic traffic, double fcg=-1)
    {
      if (!enabled_) return;
      Counters delta{};
      if (Kaskade::BDDC::TransferTiming::mode())
        delta=Kaskade::BDDC::TransferTiming::difference(
            Kaskade::BDDC::TransferTiming::snapshot(),before);
      rows.push_back({step,calls,accepted,int(Kaskade::BDDC::TransferTiming::mode()),
                     std::move(status),std::chrono::duration<double>(Clock::now()-start).count(),
                     fcg,std::move(traffic),std::move(delta)});
    }

    void write(std::string const& directory, int rank) const
    {
      if (!enabled_) return;
      std::filesystem::create_directories(directory);
      auto file=[&](char const* kind)
      {
        std::ofstream out(std::filesystem::path(directory)/
            (std::string("partition-")+kind+"-rank-"+std::to_string(rank)+".csv"));
        out.exceptions(std::ios::failbit|std::ios::badbit);
        out<<std::setprecision(17); return out;
      };
      auto bytes=file("transport"), timing=file("timing"), steps=file("steps");
      bytes<<"step,phase,scope,phase_applications,raw_reference_bytes,prepacking_bytes,"
             "sent_bytes,mpi_messages,mean_raw_reference_bytes,mean_sent_bytes,"
             "max_raw_reference_bytes,max_sent_bytes\n";
      timing<<"step,phase,timing_mode,compression_s,decompression_s,phase_wall_s,"
              "serialization_s,deserialization_s,mpi_payload_collective_s,mpi_count_collective_s,"
              "mpi_layout_validation_s,mpi_all_calls_s,float_forward_s,integer_forward_s,"
              "float_inverse_s,integer_inverse_s,quantize_s,dequantize_s,huffman_encode_s,"
              "huffman_decode_s,recenter_s,uncenter_s\n";
      steps<<"step,solver_calls,accepted_updates,status,local_step_wall_s,fcg_call_s,coarse_s\n";
      for (auto const& row : rows)
      {
        steps<<row.step<<','<<row.calls<<','<<row.accepted<<','<<row.status<<','<<row.wall<<',';
        if (row.fcg>=0) steps<<row.fcg;
        steps<<',';
        if (row.timingMode) steps<<seconds(row.counters,2,Stage::CoarseSolve);
        steps<<'\n';
        for (unsigned p=0; p<2; ++p)
        {
          auto const& data=row.traffic.phase[p];
          auto phase=p ? "prolongation" : "restriction";
          for (int scope=1; scope<=3; ++scope)
          {
            auto b=scope==3 ? combine(data.route[1],data.route[2]) : data.route[scope];
            bytes<<row.step<<','<<phase<<','
                 <<(scope==1 ? "same_node" : scope==2 ? "inter_node" : "inter_rank")
                 <<','<<data.applications<<',';
            if (b.rawReferenceKnown) bytes<<b.rawReferenceBytes;
            bytes<<','<<b.prepackingBytes<<','<<b.sentBytes<<','<<b.messages<<',';
            if (b.messages && b.rawReferenceKnown) bytes<<double(b.rawReferenceBytes)/b.messages;
            bytes<<',';
            if (b.messages) bytes<<double(b.sentBytes)/b.messages;
            bytes<<',';
            if (b.messages && b.rawReferenceKnown) bytes<<b.maxRawReferenceBytes;
            bytes<<',';
            if (b.messages) bytes<<b.maxSentBytes;
            bytes<<'\n';
          }
          timing<<row.step<<','<<phase<<','<<row.timingMode<<',';
          for (auto stage : {Stage::Compression,Stage::Decompression,
                             p ? Stage::Prolongation : Stage::Restriction})
          { if (row.timingMode) timing<<seconds(row.counters,p,stage); timing<<','; }
          timing<<data.serialization<<','<<data.deserialization<<','<<data.mpi.payload<<','
                <<data.mpi.counts<<','<<data.mpi.validation<<',';
          if (row.timingMode) timing<<seconds(row.counters,p,Stage::Mpi);
          for (auto stage : {Stage::FloatForward,Stage::IntegerForward,Stage::FloatInverse,
                             Stage::IntegerInverse,Stage::Quantize,Stage::Dequantize,
                             Stage::HuffmanEncode,Stage::HuffmanDecode,Stage::Recenter,Stage::Uncenter})
          { timing<<','; if (row.timingMode>=2) timing<<seconds(row.counters,p,stage); }
          timing<<'\n';
        }
      }
      bytes.close(); timing.close(); steps.close();
    }
  };
}
#endif
