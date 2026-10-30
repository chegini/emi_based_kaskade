#include "../EMI_material_assembly.hh"
#include "emi_material_test_matrix.hh"

#include <iostream>
#include <string>

void requireMaterial(bool ok, char const* message)
{
  if (!ok)
    throw std::runtime_error(message);
}

template <class F>
void rejectMaterial(F action)
{
  bool rejected = false;
  try { action(); }
  catch (std::runtime_error const&) { rejected = true; }
  requireMaterial(rejected,"Invalid material matrix payload accepted");
}

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
  try
  {
    MaterialTestMatrix matrix(3);
    matrix[0][0][0][0] = 1.25;
    matrix[0][2][0][0] = -0.0625;
    matrix[1][1][0][0] = std::numeric_limits<double>::denorm_min();
    matrix[2][0][0][0] = 1e-200; // do not mirror a full extracted matrix
    auto const packed = EmiBddc::packMaterialMatrix(matrix);
    requireMaterial(packed.size() == 12,"Structural zeros should not be transmitted");
    MaterialTestMatrix restored(3);
    restored = 7.0; // decoding must clear old entries
    EmiBddc::restoreMaterialMatrix(packed,restored);
    requireMaterial(EmiBddc::packMaterialMatrix(restored) == packed,
                    "Nonzero matrix bits changed in round trip");
    EmiBddc::InterfaceBlockDifference exact;
    EmiBddc::compareMaterialMatrices(restored,matrix,exact);
    requireMaterial(exact.maxAbs == 0.0,"Material round trip changed the operator");

    // Missing structural zeros are harmless; a missing/extra NONZERO is not.
    auto sparse = matrix;
    sparse.erase(2,2);
    EmiBddc::InterfaceBlockDifference zero;
    EmiBddc::compareMaterialMatrices(sparse,matrix,zero);
    requireMaterial(zero.maxAbs == 0.0,"Audit treated structural zero as numerical error");
    sparse.erase(0,0);
    EmiBddc::InterfaceBlockDifference missing;
    EmiBddc::compareMaterialMatrices(sparse,matrix,missing);
    requireMaterial(missing.maxAbs == 1.25 && !missing.acceptable(),"Audit missed missing nonzero");
    EmiBddc::InterfaceBlockDifference extra;
    EmiBddc::compareMaterialMatrices(matrix,sparse,extra);
    requireMaterial(extra.maxAbs == 1.25 && !extra.acceptable(),"Audit missed extra nonzero");
    rejectMaterial([&] { MaterialTestMatrix wrong(2); EmiBddc::compareMaterialMatrices(matrix,wrong,exact); });
    rejectMaterial([&] { EmiBddc::restoreMaterialMatrix({1},restored); });
    auto bad = packed;
    bad[0] = 3;
    rejectMaterial([&] { EmiBddc::restoreMaterialMatrix(bad,restored); });
    bad = packed;
    bad.insert(bad.end(),packed.begin(),packed.end());
    rejectMaterial([&] { EmiBddc::restoreMaterialMatrix(bad,restored); });
    bad = packed;
    bad[2] = std::bit_cast<std::uint64_t>(std::numeric_limits<double>::infinity());
    rejectMaterial([&] { EmiBddc::restoreMaterialMatrix(bad,restored); });
    auto incompatible = matrix;
    incompatible.erase(0,0);
    rejectMaterial([&] { EmiBddc::restoreMaterialMatrix(packed,incompatible); });
    auto nonfinite = matrix;
    nonfinite[0][0][0][0] = std::numeric_limits<double>::quiet_NaN();
    rejectMaterial([&] { (void)EmiBddc::packMaterialMatrix(nonfinite); });
    rejectMaterial([&] { EmiBddc::compareMaterialMatrices(nonfinite,matrix,exact); });
    EmiBddc::restoreMaterialMatrix({},restored);
    requireMaterial(EmiBddc::packMaterialMatrix(restored).empty(),"Empty material payload did not clear matrix");

    int const subdomains = std::max(3,ranks);
    for (bool masterWorker : {false,true})
    {
      if (masterWorker && ranks == 1)
        continue;
      Kaskade::BDDC::MpiDomainLayout layout;
      layout.configure(mpi,subdomains,masterWorker);
      std::vector<MaterialTestMatrix> expected(subdomains), local(subdomains);
      for (int source = 0; source < subdomains; ++source)
      {
        expected[source] = MaterialTestMatrix(3);
        if (source != 0) // zero-payload subdomain
        {
          expected[source][0][0][0][0] = 0.125*(source+1);
          expected[source][0][2][0][0] = -0.25*(source+1);
          expected[source][2][0][0][0] = 1e-200*(source+1);
        }
        if (layout.owns(source))
          local[source] = expected[source];
      }
      MaterialTestMatrix pattern(3);
      pattern = 9.0; // replicas must copy the pattern, not retain its values
      // Force multiple broadcasts even for this tiny matrix.
      EmiBddc::replicateMaterialMatrices(local,pattern,layout,7);
      for (int source = 0; source < subdomains; ++source)
      {
        requireMaterial(local[source].N() == 3,"Replica was not initialized");
        requireMaterial(EmiBddc::packMaterialMatrix(local[source]) ==
                        EmiBddc::packMaterialMatrix(expected[source]),
                        "Owner matrix exchange changed nonzero bits");
      }
#ifdef KASKADE_HAVE_MPI
      if (ranks > 1)
      {
        // A corrupt owner must not leave receivers waiting for its payload.
        if (layout.owns(0))
          local[0][0][0][0][0] = std::numeric_limits<double>::infinity();
        rejectMaterial([&] { EmiBddc::replicateMaterialMatrices(local,pattern,layout); });
        // All ranks must also observe an error on just one receiver.
        local = expected;
        auto hole = pattern;
        if (rank == (layout.owner(1)+1)%ranks)
          hole.erase(0,0);
        rejectMaterial([&] { EmiBddc::replicateMaterialMatrices(local,hole,layout); });
      }
#endif
    }
    if (rank == 0)
      std::cout << "EMI material matrix/ownership checks passed (" << ranks << " rank(s))\n";
  }
  catch (std::exception const& error)
  {
    std::cerr << "rank " << rank << ": " << error.what() << '\n';
#ifdef KASKADE_HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#else
    return 1;
#endif
  }
#ifdef KASKADE_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
