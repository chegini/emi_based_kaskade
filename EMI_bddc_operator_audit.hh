#ifndef EMI_BDDC_OPERATOR_AUDIT_HH
#define EMI_BDDC_OPERATOR_AUDIT_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

#include "EMI_bddc_audit.hh"

namespace EmiBddc
{
  namespace OperatorAuditDetail
  {
    template <class Vector>
    double norm(Vector const& x)
    {
      return std::sqrt(x.dot(x));
    }

    // A cancellation-safe scale for a sampled reciprocity test. Report the
    // signed products as well: small defect on a few probes is not an SPD proof.
    inline double relativeDefect(double a, double b, double scale)
    {
      return std::abs(a-b)/std::max(scale,std::numeric_limits<double>::min());
    }

    template <class Matrix>
    std::array<double,2> matrixSkew(Matrix const& matrix)
    {
      double skew = 0.0, magnitude = 0.0;
      for (auto row = matrix.begin(); row != matrix.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
        {
          double const value = (*col)[0][0];
          auto const& otherRow = matrix[col.index()];
          auto const other = otherRow.find(row.index());
          double const transpose = other != otherRow.end() ? (*other)[0][0] : 0.0;
          if (!std::isfinite(value) || !std::isfinite(transpose))
            return {std::numeric_limits<double>::infinity(),1.0};
          skew = std::max(skew,std::abs(value-transpose));
          magnitude = std::max(magnitude,std::abs(value));
        }
      return {skew,magnitude};
    }

    // Exercise the actual uncompressed transfer kernels, with explicit local
    // message delivery. Restriction's output accumulator starts from zero;
    // postRestrict adds residual updates rather than replacing the accumulator.
    // Only shared DOFs enter the dot products; interiors are not transfer DOFs.
    template <class Matrix, class Vector, class Interfaces>
    void auditTransfers(BddcData<Matrix,Vector> const& data, Interfaces const& interfaces)
    {
      using Transfer = Kaskade::BDDC::SpaceTransfer<1,double,double>;
      std::vector<std::unique_ptr<Transfer>> transfers;
      std::vector<Vector> y, r, projected, restricted;
      for (size_t s = 0; s < data.tags.size(); ++s)
      {
        transfers.push_back(std::make_unique<Transfer>(s,interfaces.interfaces(s)));
        Vector diagonal(data.localDofs[s].size());
        y.emplace_back(diagonal.size());
        r.emplace_back(diagonal.size());
        projected.emplace_back(diagonal.size());
        restricted.emplace_back(diagonal.size());
        for (size_t i = 0; i < diagonal.size(); ++i)
        {
          diagonal[i][0] = data.localMatrices[s][i][i][0][0];
          double const key = static_cast<double>(data.localDofs[s][i]+1);
          y[s][i][0] = std::sin(0.37*key+0.19*s);
          r[s][i][0] = std::cos(0.23*key+0.41*s);
        }
        transfers.back()->setInterfaceAveragingWeights(s,diagonal);
        transfers.back()->setGlobalIndex_active();
      }
      for (size_t to = 0; to < transfers.size(); ++to)
        for (auto const& [from,ids] : interfaces.interfaces(to))
        {
          (void)ids;
          Vector weights;
          transfers[from]->getInterfaceAveragingWeights(to,weights);
          transfers[to]->setInterfaceAveragingWeights(from,weights);
        }

      auto project = [&](std::vector<Vector> const& values)
      {
        for (size_t s = 0; s < transfers.size(); ++s)
          transfers[s]->preProlongate(values[s]);
        for (size_t to = 0; to < transfers.size(); ++to)
          for (auto const& [from,ids] : interfaces.interfaces(to))
          {
            (void)ids;
            typename Transfer::TransmissionType message;
            transfers[from]->getInterfaceValues(to,message);
            transfers[to]->setInterfaceValues(from,message);
          }
        for (size_t s = 0; s < transfers.size(); ++s)
        {
          projected[s] = 0.0;
          transfers[s]->postProlongate(values[s],projected[s]);
          projected[s] += values[s];
        }
      };
      project(y);
      for (size_t s = 0; s < transfers.size(); ++s)
        transfers[s]->preRestrict(r[s]);
      for (size_t to = 0; to < transfers.size(); ++to)
        for (auto const& [from,ids] : interfaces.interfaces(to))
        {
          (void)ids;
          typename Transfer::TransmissionType message;
          transfers[from]->getInterfaceResidualUpdate(to,message);
          transfers[to]->setInterfaceResidualUpdate(from,message);
        }
      for (size_t s = 0; s < transfers.size(); ++s)
      {
        restricted[s] = 0.0;
        transfers[s]->postRestrict(restricted[s]);
      }
      double lhs = 0.0, rhs = 0.0, r2 = 0.0, py2 = 0.0, ptr2 = 0.0, y2 = 0.0;
      for (auto const& copies : data.sharedDofs)
        for (auto const& copy : copies)
        {
          double const a = r[copy.s][copy.i][0];
          double const b = projected[copy.s][copy.i][0];
          double const c = restricted[copy.s][copy.i][0];
          double const d = y[copy.s][copy.i][0];
          lhs += a*b;
          rhs += c*d;
          r2 += a*a;
          py2 += b*b;
          ptr2 += c*c;
          y2 += d*d;
        }
      std::vector<Vector> ones(y);
      for (auto& v : ones)
        v = 1.0;
      project(ones);
      double constantError = 0.0;
      for (auto const& copies : data.sharedDofs)
        for (auto const& copy : copies)
          constantError = std::max(constantError,
                                  std::abs(projected[copy.s][copy.i][0]-1.0));
      std::cout << "  transfer <r,P*y>, <P^T*r,y>: " << lhs << ", " << rhs << '\n'
                << "  transfer adjoint defect (norm-scaled): "
                << relativeDefect(lhs,rhs,std::sqrt(r2*py2)+std::sqrt(ptr2*y2)) << '\n'
                << "  transfer constant preservation max error: " << constantError << '\n';
    }
  }

  // Diagnostic-only: factors are constructed once and used for fresh,
  // uncompressed linear preconditioner applications, never CG iterations or
  // line searches. No physical matrix, boundary term or averaging weight changes.
  // The caller exits after this audit so probe state and traffic cannot affect
  // a production timestep solve.
  template <class Matrix, class Vector, class Options>
  void auditBddcOperator(BddcData<Matrix,Vector> const& data,
                         Matrix const& matrix, Options const& options)
  {
    using namespace OperatorAuditDetail;
    using Subdomain = Kaskade::BDDC::Subdomain<1,double,double,
                                             Kaskade::BDDC::SpaceTransfer<1,double,double>>;
    bool const root = mpiRankZero(options.mpi != 0);
    size_t const n = matrix.N();
    if (n == 0 || n > static_cast<size_t>(std::numeric_limits<int>::max()))
      throw std::runtime_error("BDDC operator audit requires 1 <= N <= INT_MAX");
    Kaskade::BDDC::InterfaceAverages<1,int> interfaces(
      data.sharedDofs,data.subdomainSizes,options.bddcInterfaceTypes,data.interfaceComponents);

    // These matrices and weights are replicated application metadata. Evaluate
    // validity on every rank so all ranks stop before any collective if invalid.
    std::vector<double> sums(n,0.0), diagonalSums(n,0.0);
    size_t negative = 0, zero = 0, invalid = 0, badTransferDenominators = 0;
    double minWeight = std::numeric_limits<double>::infinity(), maxWeight = -minWeight;
    double unityError = 0.0, diagonalError = 0.0, weightMatch = 0.0;
    for (size_t s = 0; s < data.tags.size(); ++s)
      for (size_t i = 0; i < data.localDofs[s].size(); ++i)
      {
        size_t const g = data.localDofs[s][i];
        double const weight = data.weights[s][i][0];
        double const diagonal = data.localMatrices[s][i][i][0][0];
        invalid += !std::isfinite(weight) || !std::isfinite(diagonal);
        negative += weight < 0.0;
        zero += weight == 0.0;
        minWeight = std::min(minWeight,weight);
        maxWeight = std::max(maxWeight,weight);
        sums[g] += weight;
        diagonalSums[g] += diagonal;
      }
    for (size_t g = 0; g < n; ++g)
    {
      unityError = std::max(unityError,std::abs(sums[g]-1.0));
      double const diagonal = matrix[g][g][0][0];
      diagonalError = std::max(diagonalError,std::abs(diagonalSums[g]-diagonal));
    }
    for (auto const& copies : data.sharedDofs)
    {
      auto const& first = copies.front();
      size_t const g = data.localDofs[first.s][first.i];
      if (!std::isfinite(diagonalSums[g]) || diagonalSums[g] <= 0.0)
      {
        ++badTransferDenominators;
        continue;
      }
      for (auto const& copy : copies)
        weightMatch = std::max(weightMatch,std::abs(data.weights[copy.s][copy.i][0]
          - data.localMatrices[copy.s][copy.i][copy.i][0][0]/diagonalSums[g]));
    }
    if (root)
    {
      std::cout << std::setprecision(12)
                << "BDDC operator audit (diagnostic only; no timesteps):\n"
                << "  external weight min/max: " << minWeight << ", " << maxWeight << '\n'
                << "  negative weights / zero weights / nonfinite weights or diagonals: "
                << negative << ", " << zero << ", " << invalid << '\n'
                << "  max partition-of-unity error: " << unityError << '\n'
                << "  max local/global diagonal reconstruction error: " << diagonalError << '\n'
                << "  max external/transfer normalized weight difference: " << weightMatch << '\n'
                << "  shared DOFs with nonpositive/nonfinite transfer denominator: "
                << badTransferDenominators << '\n';
    }
    if (invalid || badTransferDenominators)
      throw std::runtime_error("BDDC operator audit: invalid weights; skipping unsafe probes");
    if (root)
    {
      auto const globalSkew = matrixSkew(matrix);
      double localSkew = 0.0;
      for (auto const& local : data.localMatrices)
      {
        auto const skew = matrixSkew(local);
        localSkew = std::max(localSkew,skew[0]/std::max(skew[1],
                                                    std::numeric_limits<double>::min()));
      }
      std::cout << "  global max-entry relative matrix skew: "
                << globalSkew[0]/std::max(globalSkew[1],std::numeric_limits<double>::min())
                << "\n  maximum local max-entry relative matrix skew: " << localSkew << '\n';
      auditTransfers(data,interfaces);
    }
    auto subdomains = constructOwnedSubdomains<Subdomain>(
      data.localMatrices,interfaces,options.assemblyThreads,options.mpi != 0,
      options.masterWorker != 0);
    using Domain = Kaskade::BDDC::SharedMemoryDomain<Subdomain,decltype(subdomains)>;
    Domain domain(subdomains,options.mpi != 0,options.masterWorker != 0);
    std::vector<int> active(data.tags.size());
    std::iota(active.begin(),active.end(),0);
    domain.setActiveSubdomains(active);
    domain.set_solver_type(false);
    Kaskade::BDDC::KKTSolver<Domain> coarse(domain,interfaces.coarseConstraints());

    // Reconstruct 1 from the local constrained energy-minimizing basis:
    // solve [A_i C_i^T; C_i 0] [z;lambda] = [0; C_i*1].
    // The Robin term means exact reproduction is not imposed as an assertion.
    double constantL2 = 0.0, constantMax = 0.0, constraintDefect = 0.0;
    int worstTag = -1;
    for (size_t s = 0; s < subdomains.size(); ++s)
    {
      if (!subdomains.has(s))
        continue;
      auto const C = interfaces.template coarseConstraint<double>(s);
      typename Subdomain::Vector values(C.N()), ones(data.localDofs[s].size());
      typename Subdomain::XVector reconstructed(ones.size());
      ones = 1.0;
      values = 0.0;
      C.mv(ones,values);
      subdomains[s].resetForNewRhs();
      subdomains[s].setCorrection(values);
      subdomains[s].getCorrection(reconstructed);
      double l2 = 0.0;
      for (size_t i = 0; i < reconstructed.size(); ++i)
      {
        double const e = reconstructed[i][0]-1.0;
        if (!std::isfinite(e))
        {
          constantMax = std::numeric_limits<double>::infinity();
          l2 = std::numeric_limits<double>::infinity();
        }
        else
        {
          constantMax = std::max(constantMax,std::abs(e));
          l2 += e*e;
        }
      }
      l2 = std::sqrt(l2/std::max<size_t>(1,reconstructed.size()));
      if (l2 > constantL2)
      {
        constantL2 = l2;
        worstTag = data.tags[s];
      }
      typename Subdomain::Vector check(C.N());
      C.mv(reconstructed,check);
      for (size_t i = 0; i < check.size(); ++i)
        if (!std::isfinite(check[i][0]))
          constraintDefect = std::numeric_limits<double>::infinity();
        else
          constraintDefect = std::max(constraintDefect,std::abs(check[i][0]-values[i][0])
                                     /std::max(1.0,std::abs(values[i][0])));
    }
    // Rank-labelled reports retain the worst local tag without a second
    // distributed selection rule. Root prints the remaining global probes.
    int rank = 0;
#ifdef KASKADE_HAVE_MPI
    if (options.mpi)
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    std::ostringstream localReport;
    localReport << std::setprecision(12)
                << "  rank " << rank
                << " owned-subdomain coarse constant reconstruction: max relative L2="
                << constantL2 << " (tag " << worstTag << "), max absolute="
                << constantMax << ", max normalized constraint defect=" << constraintDefect << '\n';
    std::cout << localReport.str();

    auto apply = [&](Vector const& rhs)
    {
      auto const localRhs = scatterToSubdomains(data,rhs);
      for (size_t s = 0; s < subdomains.size(); ++s)
        if (subdomains.has(s))
        {
          subdomains[s].resetForNewRhs();
          subdomains[s].setRhs(localRhs[s]);
        }
      domain.updateIteration(0);
      domain.restrict();
      coarse.solve();
      (void)domain.prolongate();
      std::vector<double> local(n,0.0), global(n,0.0);
      for (size_t s = 0; s < subdomains.size(); ++s)
        if (subdomains.has(s))
        {
          typename Subdomain::XVector interior(data.localDofs[s].size()), correction(interior.size());
          subdomains[s].getSolution(interior);
          subdomains[s].getCorrection(correction);
          for (size_t i = 0; i < interior.size(); ++i)
            local[data.localDofs[s][i]] += data.weights[s][i][0]
                                        * (interior[i][0]+correction[i][0]);
        }
#ifdef KASKADE_HAVE_MPI
      if (options.mpi)
        MPI_Allreduce(local.data(),global.data(),static_cast<int>(n),
                      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      else
#endif
        global = local;
      Vector output(n);
      for (size_t i = 0; i < n; ++i)
        output[i][0] = global[i];
      return output;
    };
    Vector x(n), y(n), ones(n);
    ones = 1.0;
    for (size_t i = 0; i < n; ++i)
    {
      x[i][0] = std::sin(0.73*(i+1))+std::cos(0.17*(i+3));
      y[i][0] = std::sin(0.29*(i+2))+std::cos(0.61*(i+1));
    }
    for (auto* probe : {&x,&y})
    {
      probe->axpy(-ones.dot(*probe)/static_cast<double>(n),ones);
      *probe *= 1.0/std::max(norm(*probe),std::numeric_limits<double>::min());
    }
    Vector const bx = apply(x), by = apply(y);
    Vector sum(x);
    sum += y;
    Vector linearity = apply(sum);
    linearity -= bx;
    linearity -= by;
    Vector repeat = apply(x);
    repeat -= bx;
    Vector aOne(n);
    matrix.mv(ones,aOne);
    Vector constant = apply(aOne);
    Vector constantError(constant);
    constantError -= ones;
    if (root)
      std::cout << "  <x,B*y>, <y,B*x>: " << x.dot(by) << ", " << y.dot(bx) << '\n'
                << "  preconditioner reciprocity defect (norm-scaled): "
                << relativeDefect(x.dot(by),y.dot(bx),norm(x)*norm(by)+norm(y)*norm(bx)) << '\n'
                << "  preconditioner reciprocity defect (pair-scaled): "
                << relativeDefect(x.dot(by),y.dot(bx),std::abs(x.dot(by))+std::abs(y.dot(bx)))
                << "\n  mean(B*x), mean(B*y): "
                << ones.dot(bx)/static_cast<double>(n) << ", "
                << ones.dot(by)/static_cast<double>(n) << '\n'
                << "  sampled quadratic forms <x,B*x>, <y,B*y>: "
                << x.dot(bx) << ", " << y.dot(by) << '\n'
                << "  preconditioner additivity defect: "
                << norm(linearity)/std::max(norm(bx)+norm(by),std::numeric_limits<double>::min())
                << "\n  repeat-after-reset relative defect: "
                << norm(repeat)/std::max(norm(bx),std::numeric_limits<double>::min())
                << "\n  ||B*A*1-1||_2 / ||1||_2: "
                << norm(constantError)/std::sqrt(static_cast<double>(n))
                << "\n  mean(B*A*1): " << ones.dot(constant)/static_cast<double>(n)
                << "\n  B is one fresh unscaled preconditioner application including the interior solve."
                << "\n  Sampled symmetry/positivity and constant-mode checks are diagnostics, not proofs."
                << "\nBDDC operator audit completed; no timestep solution or final VTK was produced.\n";
  }
}

#endif
