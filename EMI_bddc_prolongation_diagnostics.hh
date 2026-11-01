#ifndef EMI_BDDC_PROLONGATION_DIAGNOSTICS_HH
#define EMI_BDDC_PROLONGATION_DIAGNOSTICS_HH

#include "EMI_bddc_cg_diagnostics.hh"
#include <set>
#include <stdexcept>

namespace EmiBddc
{
  struct ProlongationConstraintGroup
  {
    std::vector<size_t> dofs;
    size_t owners = 0;
  };

  // Read the real C rows, not pairwise communication interfaces, which overlap.
  // The EMI scalar InterfaceAverages construction uses disjoint rows of ones.
  template<class Interfaces, class Data>
  std::vector<ProlongationConstraintGroup> prolongationConstraintGroups(
      Interfaces const& interfaces, Data const& data)
  {
    auto const constraints = interfaces.coarseConstraints();
    std::vector<ProlongationConstraintGroup> groups(constraints.size());
    std::vector<std::vector<size_t>> ids(data.localDofs.size());
    for (size_t g=0; g<constraints.size(); ++g)
    {
      groups[g].owners = constraints[g].first.size();
      for (int s: constraints[g].first) ids.at(s).push_back(g);
    }
    for (size_t s=0; s<ids.size(); ++s)
    {
      auto const C = interfaces.template coarseConstraint<double>(s);
      if (C.N()!=ids[s].size()) throw std::runtime_error("Constraint probe: row count mismatch");
      for (auto row=C.begin(); row!=C.end(); ++row)
      {
        std::vector<size_t> dofs;
        for (auto col=row->begin(); col!=row->end(); ++col)
        {
          if ((*col)[0][0]!=1) throw std::runtime_error("Constraint probe needs rows of ones");
          dofs.push_back(data.localDofs[s].at(col.index()));
        }
        std::sort(dofs.begin(),dofs.end());
        auto& group=groups.at(ids[s].at(row.index()));
        if (group.dofs.empty()) group.dofs=dofs;
        else if (group.dofs!=dofs)
          throw std::runtime_error("Constraint probe: owner supports differ");
      }
    }
    return groups;
  }

  /**
   * Separate the actual correction into exact-prolongation shadow plus error.
   * Both use the SAME local/coarse result after the actual restricted residual.
   * The shadow includes exact interface averaging and interior harmonic extension;
   * it is never used by the numerical solver. This does not simulate exact restriction.
   * Two diagnostic global sum reductions and extra local Dirichlet solves are
   * deliberately outside transfer byte counters. Timings are not benchmarks.
   */
  template<class Matrix, class Vector, class Storage, class Data>
  class ProlongationDiagnosticReport
  {
  public:
    using Sum = std::function<void(std::vector<double>&)>;
    ProlongationDiagnosticReport(Matrix const& matrix, Vector const& rhs,
                                 Storage const& storage, Data const& data, bool root, Sum sum)
      : matrix_(matrix), rhs_(rhs), storage_(storage), data_(data), root_(root),
        sum_(std::move(sum)), n_(rhs.size()) {}

    bool open(std::string const& path)
    {
      if (!root_) return true;
      output_.open(path);
      output_ << "solver_call,restriction_pass,true_residual_norm,restriction_error_rel,"
                 "rho_before_p,rho_after_p,rho_change_p,rho_closure_error,"
                 "z_before_norm,z_after_norm,p_error_norm,p_error_rel,p_error_energy,"
                 "actual_z_jump,weight_partition_defect,inactive_global_dofs,rho_internal\n";
      output_ << std::setprecision(17);
      return output_.good();
    }

    // Enabled explicitly for one solver call, on step one only by the caller.
    void constraintProbe(std::vector<ProlongationConstraintGroup> groups, int call,
                         std::string const& prefix)
    {
      std::set<size_t> seen;
      for (auto const& group: groups)
      {
        if (group.dofs.empty() || group.owners<2)
          throw std::runtime_error("Constraint probe: empty or unshared group");
        if (!std::is_sorted(group.dofs.begin(),group.dofs.end()))
          throw std::runtime_error("Constraint probe: supports must be sorted");
        for (size_t g: group.dofs)
          if (g>=n_ || !seen.insert(g).second)
            throw std::runtime_error("Constraint probe requires disjoint valid supports");
      }
      groups_=std::move(groups); probeCall_=call; probePrefix_=prefix;
    }

    void groupingProbe(int call, int currentBits, std::string const& path)
    {
      groupingCall_=call; groupingBits_=currentBits; groupingPath_=path;
    }

    template<class Event>
    void observe(Event const& event)
    {
      if (event.restriction) return;
      std::vector<double> global(5*n_,0);
      for (size_t s=0; s<storage_.size(); ++s)
      {
        if (!storage_.has(s)) continue;
        auto const& sub = storage_[s];
        auto const x = sub.getSolution();
        auto const& raw = sub.diagnosticRawCorrection();
        auto const weights = sub.diagnosticProlongationWeights();
        auto const actual = sub.diagnosticProlongatedCorrection();
        auto const& rr = sub.getRestrictedResidual();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
        {
          size_t const g = data_.localDofs[s][i];
          double const w = data_.weights[s][i][0];
          global[g] += w*x[i][0];
          global[n_+g] += rr[i][0];
          global[2*n_+g] += weights[i][0]*raw[i][0];
          global[3*n_+g] += w*actual[i][0];
          global[4*n_+g] += weights[i][0];
        }
      }
      sum_(global);
      std::vector<double> shadow(n_+1,0);
      for (size_t s=0; s<storage_.size(); ++s)
      {
        if (!storage_.has(s)) continue;
        auto const& sub = storage_[s];
        // Subdomain vectors are scalar-blocked, independently of the application's block type.
        auto averaged = sub.diagnosticRawCorrection();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
          averaged[i][0] = global[2*n_+data_.localDofs[s][i]];
        auto const exact = sub.diagnosticExactProlongation(averaged);
        auto const actual = sub.diagnosticProlongatedCorrection();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
        {
          size_t const g = data_.localDofs[s][i];
          shadow[g] += data_.weights[s][i][0]*exact[i][0];
          double const d = actual[i][0]-global[3*n_+g];
          shadow[n_] += d*d;
        }
      }
      sum_(shadow);
      if (event.solverCall==probeCall_) probe(global,shadow,event.solverCall);
      if (event.solverCall==groupingCall_) probeGrouping(global,shadow,event.solverCall);
      if (!root_) return;
      using namespace CgDiagnostics;
      std::vector<double> x(global.begin(),global.begin()+n_);
      auto Ax = multiply(x);
      std::vector<double> r(n_), re(n_), before(shadow.begin(),shadow.begin()+n_),
                          after(global.begin()+3*n_,global.begin()+4*n_), error(n_);
      double weightDefect = 0;
      size_t inactive = 0;
      for (size_t i=0; i<n_; ++i)
      {
        r[i] = rhs_[i][0]-Ax[i];
        re[i] = global[n_+i]-r[i];
        error[i] = after[i]-before[i];
        double const w = global[4*n_+i];
        if (w==0) ++inactive;
        else weightDefect = std::max(weightDefect,std::abs(w-1));
      }
      double const beforeRho = dot(r,before), afterRho = dot(r,after), deltaRho = dot(r,error);
      output_ << event.solverCall << ',' << event.pass;
      for (double v : {norm(r),ratio(norm(re),norm(r)),beforeRho,afterRho,deltaRho,
                       std::abs(afterRho-beforeRho-deltaRho),norm(before),norm(after),
                       norm(error),ratio(norm(error),norm(before)),dot(error,multiply(error)),
                       std::sqrt(shadow[n_])/std::max(1.0,norm(after)),weightDefect})
      {
        output_ << ',';
        if (std::isfinite(v)) output_ << v;
      }
      output_ << ',' << inactive << ',';
      if (std::isfinite(event.rho)) output_ << event.rho;
      output_ << '\n';
      output_.flush(); // Preserve the rejected proposal even if the solver subsequently aborts.
    }

    bool good() const { return !root_ || (output_.good() && probeGood_); }
    bool close()
    {
      if (!root_) return true;
      output_.flush();
      bool const ok = output_.good();
      output_.close();
      return ok && !output_.fail() && probeGood_;
    }

  private:
    void probeGrouping(std::vector<double> const& global, std::vector<double> const& shadow, int call)
    {
      using namespace CgDiagnostics;
      std::ofstream out;
      if (root_)
      {
        out.open(groupingPath_);
        out << std::setprecision(17)
            << "solver_call,bits,group_size,current_precision,rho_exact,rho_actual,rho_candidate,"
               "relative_error_vs_exact,relative_difference_vs_actual,producer_groups,"
               "metadata_bytes,stream_bytes,padding_bytes,exact_bytes,modeled_p_bytes,"
               "packed_p_bytes,logical_p_messages,codec_valid\n";
      }
      std::vector<std::pair<int,size_t>> candidates{{groupingBits_,0}};
      for (auto candidate: std::vector<std::pair<int,size_t>>{
          {16,0},{8,0},{8,8},{8,16},{8,32},{8,64},{8,128},{8,256},{8,512}})
        if (std::find(candidates.begin(),candidates.end(),candidate)==candidates.end())
          candidates.push_back(candidate);
      for (auto const& [bits,groupSize]: candidates)
      {
        // All candidates use the SAME stored raw correction from the original trajectory.
        std::vector<double> packed(n_+9,0);
        for (size_t s=0; s<storage_.size(); ++s)
        {
          if (!storage_.has(s)) continue;
          auto const& sub=storage_[s];
          if constexpr (requires { sub.transfer().diagnosticGroupedProlongation(
                                      sub.diagnosticRawCorrection().size(),bits,groupSize); })
          {
            auto const sample=sub.transfer().diagnosticGroupedProlongation(
                sub.diagnosticRawCorrection().size(),bits,groupSize);
            for (size_t i=0; i<data_.localDofs[s].size(); ++i)
              packed[data_.localDofs[s][i]]+=sample.contribution[i][0];
            size_t offset=n_;
            for (auto v: {sample.groups,sample.metadataBytes,sample.streamBytes,sample.paddingBytes,
                          sample.exactBytes,sample.totalBytes,sample.packedBytes,sample.messages})
              packed[offset++]+=v;
            if (!sample.codecValid) packed[n_+8]+=1;
          }
          else packed[n_+8]+=1;
        }
        sum_(packed);
        std::vector<double> candidate(n_,0);
        for (size_t s=0; s<storage_.size(); ++s)
        {
          if (!storage_.has(s)) continue;
          auto const& sub=storage_[s];
          auto boundary=sub.diagnosticRawCorrection();
          for (size_t i=0; i<data_.localDofs[s].size(); ++i)
            boundary[i][0]=packed[data_.localDofs[s][i]];
          auto const lifted=sub.diagnosticExactProlongation(boundary);
          for (size_t i=0; i<data_.localDofs[s].size(); ++i)
            candidate[data_.localDofs[s][i]]+=data_.weights[s][i][0]*lifted[i][0];
        }
        sum_(candidate);
        if (!root_) continue;
        std::vector<double> x(global.begin(),global.begin()+n_), r=multiply(x),
                            exact(shadow.begin(),shadow.begin()+n_),
                            actual(global.begin()+3*n_,global.begin()+4*n_), error(n_), difference(n_);
        for (size_t i=0; i<n_; ++i)
        {
          r[i]=rhs_[i][0]-r[i]; error[i]=candidate[i]-exact[i];
          difference[i]=candidate[i]-actual[i];
        }
        out << call << ',' << bits << ',' << groupSize << ','
            << (bits==groupingBits_ && groupSize==0) << ',' << dot(r,exact) << ',' << dot(r,actual)
            << ',' << dot(r,candidate) << ',' << ratio(norm(error),norm(exact))
            << ',' << ratio(norm(difference),norm(actual));
        for (size_t i=0; i<8; ++i) out << ',' << packed[n_+i];
        out << ',' << (packed[n_+8]==0) << '\n';
      }
      if (root_) {out.close(); probeGood_=probeGood_ && !out.fail();}
    }

    void probe(std::vector<double> const& global, std::vector<double> const& shadow, int call)
    {
      using namespace CgDiagnostics;
      std::vector<double> r(n_,0);
      if (root_)
      {
        std::vector<double> x(global.begin(),global.begin()+n_);
        r=multiply(x);
        for (size_t i=0; i<n_; ++i) r[i]=rhs_[i][0]-r[i];
      }
      sum_(r); // Distribute root's original global residual for local attribution.
      size_t const ng=groups_.size();
      std::vector<double> means(ng,0);
      for (size_t g=0; g<ng; ++g)
        for (size_t i: groups_[g].dofs)
          means[g]+=(global[3*n_+i]-global[2*n_+i])/groups_[g].dofs.size();
      // Local harmonic lifts give each group's FULL r.H(e), including interiors.
      // Two temporary solves per incident group; never publish the lifted vectors.
      std::vector<double> sums(2*n_+2*ng,0);
      for (size_t s=0; s<storage_.size(); ++s)
      {
        if (!storage_.has(s)) continue;
        auto const& sub=storage_[s];
        auto averaged=sub.diagnosticRawCorrection();
        for (size_t i=0; i<data_.localDofs[s].size(); ++i)
          averaged[i][0]=global[2*n_+data_.localDofs[s][i]];
        auto const exact=sub.diagnosticExactProlongation(averaged);
        for (size_t g=0; g<ng; ++g)
        {
          auto errorBoundary=averaged, meanBoundary=averaged;
          bool incident=false;
          for (size_t i=0; i<data_.localDofs[s].size(); ++i)
          {
            size_t const k=data_.localDofs[s][i];
            if (!std::binary_search(groups_[g].dofs.begin(),groups_[g].dofs.end(),k)) continue;
            incident=true;
            errorBoundary[i][0]=global[3*n_+k];
            meanBoundary[i][0]+=means[g];
          }
          if (!incident) continue;
          auto const liftedError=sub.diagnosticExactProlongation(errorBoundary);
          auto const liftedMean=sub.diagnosticExactProlongation(meanBoundary);
          double const rounded=static_cast<double>(static_cast<float>(means[g]));
          double const factor=means[g]!=0 ? rounded/means[g] : 0;
          for (size_t i=0; i<data_.localDofs[s].size(); ++i)
          {
            size_t const k=data_.localDofs[s][i];
            double const w=data_.weights[s][i][0];
            double const e=liftedError[i][0]-exact[i][0];
            double const m=liftedMean[i][0]-exact[i][0];
            sums[k]+=w*m;
            sums[n_+k]+=w*m*factor;
            sums[2*n_+g]+=w*r[k]*e;
            sums[2*n_+ng+g]+=w*r[k]*m;
          }
        }
      }
      sum_(sums);
      if (!root_) return;
      std::vector<double> before(shadow.begin(),shadow.begin()+n_),
                          after(global.begin()+3*n_,global.begin()+4*n_), corrected=after,
                          corrected32=after;
      for (size_t i=0; i<n_; ++i)
      {
        corrected[i]-=sums[i]; corrected32[i]-=sums[n_+i];
      }
      std::ofstream detail(probePrefix_+"-groups.csv"), summary(probePrefix_+"-summary.csv"),
                    topology(probePrefix_+"-dofs.csv");
      topology << "group,global_dof\n";
      detail << std::setprecision(17);
      summary << std::setprecision(17);
      detail << "group,dofs,owners,exact_p_mean,actual_p_mean,mean_error,"
                "min_error,max_error,rho_error_with_interior,rho_mean_error_with_interior,"
                "mean_error_after_shadow,mean_error_after_float32_shadow,payload32,payload64\n";
      double attributed=0, meanDefect=0, meanDefect32=0;
      size_t payload32=0, payload64=0, covered=0;
      for (size_t g=0; g<ng; ++g)
      {
        auto const& group=groups_[g];
        double lo=std::numeric_limits<double>::infinity(), hi=-lo;
        double exactMean=0, actualMean=0, defect=0, defect32=0;
        for (size_t k: group.dofs)
        {
          topology << g << ',' << k << '\n';
          double const e=after[k]-global[2*n_+k];
          lo=std::min(lo,e); hi=std::max(hi,e);
          exactMean+=global[2*n_+k]/group.dofs.size();
          actualMean+=after[k]/group.dofs.size();
          defect+=(corrected[k]-before[k])/group.dofs.size();
          defect32+=(corrected32[k]-before[k])/group.dofs.size();
        }
        // Hypothetical scalar payload: one contribution per owner to every other owner.
        // Counts logical deliveries even intra-rank; excludes headers/IDs/new messages.
        size_t const deliveries=group.owners*(group.owners-1);
        payload32+=4*deliveries; payload64+=8*deliveries; covered+=group.dofs.size();
        attributed+=sums[2*n_+g];
        meanDefect=std::max(meanDefect,std::abs(defect));
        meanDefect32=std::max(meanDefect32,std::abs(defect32));
        detail << g << ',' << group.dofs.size() << ',' << group.owners << ',' << exactMean
               << ',' << actualMean << ',' << means[g] << ',' << lo << ',' << hi << ','
               << sums[2*n_+g] << ',' << sums[2*n_+ng+g] << ',' << defect << ',' << defect32
               << ',' << 4*deliveries << ',' << 8*deliveries << '\n';
      }
      double const rhoBefore=dot(r,before), rhoAfter=dot(r,after);
      summary << "solver_call,groups,covered_dofs,rho_before_p,rho_after_p,rho_mean_corrected,"
                 "rho_float32_mean_corrected,attributed_rho_error,unattributed_rho_error,"
                 "max_remaining_mean_error,max_remaining_float32_mean_error,"
                 "estimated_payload32_per_call,estimated_payload64_per_call\n";
      summary << call << ',' << ng << ',' << covered << ',' << rhoBefore << ',' << rhoAfter
              << ',' << dot(r,corrected) << ',' << dot(r,corrected32) << ',' << attributed
              << ',' << rhoAfter-rhoBefore-attributed << ',' << meanDefect << ',' << meanDefect32
              << ',' << payload32 << ',' << payload64 << '\n';
      detail.close(); summary.close(); topology.close();
      probeGood_=!detail.fail() && !summary.fail() && !topology.fail();
    }

    std::vector<double> multiply(std::vector<double> const& x) const
    {
      Vector in(n_), out(n_);
      for (size_t i=0; i<n_; ++i) in[i][0] = x[i];
      matrix_.mv(in,out);
      std::vector<double> result(n_);
      for (size_t i=0; i<n_; ++i) result[i] = out[i][0];
      return result;
    }
    Matrix const& matrix_;
    Vector const& rhs_;
    Storage const& storage_;
    Data const& data_;
    bool root_;
    Sum sum_;
    size_t n_;
    std::ofstream output_;
    std::vector<ProlongationConstraintGroup> groups_;
    int probeCall_=-1;
    std::string probePrefix_;
    bool probeGood_=true;
    int groupingCall_=-1, groupingBits_=0;
    std::string groupingPath_;
  };
}
#endif
