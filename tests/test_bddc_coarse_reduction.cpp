// Real-MPI regression for the production coarse packing and reduction helpers.
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#define KASKADE_HAVE_MPI
#include "mg/bddcCoarseReduction.hh"

namespace
{
  int calls = 0;
  int coefficients = 0;
  void require(bool ok, char const* message)
  {
    if (!ok) throw std::runtime_error(message);
  }
  double scalarSum(double local)
  {
    double global = 0;
    PMPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return global;
  }
}

// Count actual production calls, leaving the scalar oracle uninstrumented.
extern "C" int MPI_Allreduce(void const* send, void* receive, int count,
                             MPI_Datatype type, MPI_Op op, MPI_Comm comm)
{
  ++calls;
  coefficients += count;
  return PMPI_Allreduce(send,receive,count,type,op,comm);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  try
  {
    using Block = std::array<std::array<double,1>,1>;
    using Matrix = std::vector<std::vector<Block>>;
    for (int limit : {2,std::numeric_limits<int>::max()})
      for (int repeat = 0; repeat < 3; ++repeat)
      {
        Matrix matrix{std::vector<Block>(3),{},std::vector<Block>(2)};
        Matrix expected = matrix;
        int position = 0;
        for (size_t row = 0; row < matrix.size(); ++row)
          for (size_t col = 0; col < matrix[row].size(); ++col)
          {
            // Rank zero has no owned contribution, as in master/worker mode.
            // Keep a structural zero and different signs to expose reordering.
            double value = rank == 0 || position == 1 ? 0.0
              : (rank+1)*(position%2 ? -0.25 : 0.5)*(repeat+1);
            matrix[row][col][0][0] = value;
            expected[row][col][0][0] = scalarSum(value);
            ++position;
          }
        calls = coefficients = 0;
        Kaskade::BDDC::reduceCoarseMatrix(matrix,limit);
        require(matrix == expected,"Coarse matrix differs from scalar-reduction oracle");
        require(calls == (limit == 2 ? 3 : 1) && coefficients == 5,
                "Matrix reduction was not batched/chunked correctly");

        std::vector<std::array<double,1>> rhs(7), expectedRhs(7);
        for (size_t i = 0; i < rhs.size(); ++i)
        {
          rhs[i][0] = rank == 0 ? 0.0 : (rank+1)*(i+1)*0.125*(repeat+1);
          expectedRhs[i][0] = scalarSum(rhs[i][0]);
        }
        calls = coefficients = 0;
        Kaskade::BDDC::reduceCoarseRhs(rhs,limit);
        require(rhs == expectedRhs,"Coarse RHS differs from scalar-reduction oracle");
        require(calls == (limit == 2 ? 4 : 1) && coefficients == 7,
                "RHS reduction was not batched/chunked correctly");
      }
    Matrix empty;
    std::vector<std::array<double,1>> noRhs;
    calls = coefficients = 0;
    Kaskade::BDDC::reduceCoarseMatrix(empty);
    Kaskade::BDDC::reduceCoarseRhs(noRhs);
    require(calls == 0,"Empty coarse system entered a collective");
    std::vector<double> values{1.0};
    bool rejected = false;
    try { Kaskade::BDDC::allreduceCoarseSum(values,0); }
    catch (std::invalid_argument const&) { rejected = true; }
    require(rejected,"Invalid chunk limit accepted");
    if (rank == 0)
      std::cout << "PASS: coarse matrix/RHS scalar parity, inactive owner, structural zero, "
                   "repeated RHS, chunk boundaries, empty system, collective counts\n";
  }
  catch (std::exception const& e)
  {
    std::cerr << "rank " << rank << ": " << e.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  MPI_Finalize();
}
