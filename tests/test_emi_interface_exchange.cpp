#include "../EMI_interface_exchange.hh"

#include <algorithm>
#include <iostream>

// Requires MPI, but not Dune. Exercises the production sparse collectives
// with normal and master-worker ownership, including empty sender buffers.
int main(int argc, char** argv)
{
#ifndef KASKADE_HAVE_MPI
  (void)argc;
  (void)argv;
  throw std::runtime_error("Compile this test with -DKASKADE_HAVE_MPI and mpicxx");
#else
  MPI_Init(&argc,&argv);
  int rank = 0, ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try
  {
    int const subdomains = std::max(2,ranks);
    for (bool masterWorker : {false,true})
    {
      if (masterWorker && ranks == 1)
        continue;
      Kaskade::BDDC::MpiDomainLayout layout;
      layout.configure(true,subdomains,masterWorker);
      EmiBddc::InterfaceBlocks expected, local;
      // Only subdomain 0 and 1 have interface mass: other workers send zero
      // words, while all ranks must receive the same replicated operators.
      for (int source = 0; source < 2; ++source)
      {
        auto& block = expected[{source,1-source}];
        block[{0,0}] = 0.125*(source+1);
        block[{0,1}] = -0.0625*(source+1);
        block[{1,0}] = -0.0625*(source+1);
        block[{1,1}] = 1e-30*(source+1);
        if (layout.owns(source))
          local[{source,1-source}] = block;
      }
      auto const actual = EmiBddc::exchangeInterfaceBlocks(std::move(local),layout,subdomains,2);
      if (actual != expected)
        throw std::runtime_error("Sparse MPI exchange changed interface blocks");
      // One failing worker must be reported by every rank before exchange.
      std::exception_ptr error;
      if (rank == ranks-1)
        error = std::make_exception_ptr(std::runtime_error("injected assembly failure"));
      bool failed = false;
      try { EmiBddc::checkInterfaceAssemblyError(error,layout); }
      catch (std::runtime_error const&) { failed = true; }
      if (!failed)
        throw std::runtime_error("Collective assembly error was not propagated");
    }
    if (rank == 0)
      std::cout << "EMI sparse interface MPI exchange checks passed (" << ranks << " ranks)\n";
  }
  catch (std::exception const& error)
  {
    std::cerr << "rank " << rank << ": " << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  MPI_Finalize();
#endif
}
