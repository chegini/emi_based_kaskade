#ifndef EMI_INTERFACE_EXCHANGE_HH
#define EMI_INTERFACE_EXCHANGE_HH

#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

#include "mg/bddcMpi.hh"
#include "EMI_interface_blocks.hh"

namespace EmiBddc
{
  // All ranks must reach the same decision before entering payload exchange.
  // A failure on another rank must not leave a peer waiting in Allgatherv.
  inline void checkInterfaceAssemblyError(std::exception_ptr error,
                                          Kaskade::BDDC::MpiDomainLayout const& layout,
                                          char const* context = "Interface assembly")
  {
    int failed = error ? 1 : 0;
#ifdef KASKADE_HAVE_MPI
    if (layout.enabled() && layout.size() > 1)
    {
      int anyFailed = 0;
      MPI_Allreduce(&failed,&anyFailed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      failed = anyFailed;
    }
#else
    (void)layout;
#endif
    if (error)
      std::rethrow_exception(error);
    if (failed)
      throw std::runtime_error(std::string(context) + " failed on another MPI rank");
  }

  // Each directed face is integrated only by its source material's owner.
  // Current BddcData keeps replicated local matrices for diagonal weights,
  // dumps and audits, so sparse blocks are gathered to all ranks (including
  // both subdomain owners). Fine/coarse solver ownership remains unchanged.
  // Retaining these replicas avoids a second, unrelated BddcData redesign.
  inline InterfaceBlocks exchangeInterfaceBlocks(InterfaceBlocks local,
                                                  Kaskade::BDDC::MpiDomainLayout const& layout,
                                                  size_t subdomains, size_t dofs)
  {
#ifdef KASKADE_HAVE_MPI
    if (layout.enabled() && layout.size() > 1)
    {
      std::vector<std::uint64_t> send;
      std::exception_ptr error;
      try
      {
        send = packInterfaceBlocks(local);
        if (send.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
          throw std::runtime_error("Sparse interface send buffer exceeds MPI int count");
      }
      catch (...) { error = std::current_exception(); }
      checkInterfaceAssemblyError(error,layout);

      int const count = static_cast<int>(send.size());
      std::vector<int> counts(layout.size()), offsets(layout.size());
      MPI_Allgather(&count,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD);
      size_t total = 0;
      for (int rank = 0; rank < layout.size(); ++rank)
      {
        if (counts[rank] < 0 || counts[rank]%5 != 0 ||
            total + static_cast<size_t>(counts[rank]) >
              static_cast<size_t>(std::numeric_limits<int>::max()))
          throw std::runtime_error("Sparse interface receive buffer exceeds MPI int count");
        offsets[rank] = static_cast<int>(total);
        total += counts[rank];
      }
      std::vector<std::uint64_t> received(total);
      MPI_Allgatherv(send.data(),count,MPI_UINT64_T,received.data(),counts.data(),
                     offsets.data(),MPI_UINT64_T,MPI_COMM_WORLD);
      InterfaceBlocks result;
      for (int rank = 0; rank < layout.size(); ++rank)
        unpackInterfaceBlocks(received,offsets[rank],offsets[rank]+counts[rank],rank,
                               subdomains,dofs,[&layout](int s) { return layout.owner(s); },result);
      return result;
    }
#else
    (void)layout;
    (void)subdomains;
    (void)dofs;
#endif
    return local;
  }

}

#endif
