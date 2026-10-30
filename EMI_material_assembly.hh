#ifndef EMI_MATERIAL_ASSEMBLY_HH
#define EMI_MATERIAL_ASSEMBLY_HH

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "EMI_interface_exchange.hh"

namespace EmiBddc
{
  // Full scalar matrices have already been extracted from the legacy
  // assembler, including its symmetry convention. Exchange only nonzeros:
  // each entry is (global row, global column, exact double bits). There is
  // no summation, reorientation, or second quadrature evaluation in transit.
  template <class Matrix>
  std::vector<std::uint64_t> packMaterialMatrix(Matrix const& matrix)
  {
    std::vector<std::uint64_t> words;
    for (auto row = matrix.begin(); row != matrix.end(); ++row)
      for (auto col = row.begin(); col != row.end(); ++col)
      {
        double const value = (*col)[0][0];
        if (!std::isfinite(value))
          throw std::runtime_error("Nonfinite material matrix entry");
        if (value != 0.0)
          words.insert(words.end(),{static_cast<std::uint64_t>(row.index()),
                                    static_cast<std::uint64_t>(col.index()),
                                    std::bit_cast<std::uint64_t>(value)});
      }
    return words;
  }

  template <class Matrix>
  void restoreMaterialMatrix(std::vector<std::uint64_t> const& words, Matrix& matrix)
  {
    if (words.size()%3 != 0)
      throw std::runtime_error("Invalid material matrix message length");
    matrix = 0.0;
    std::pair<std::uint64_t,std::uint64_t> previous{};
    bool first = true;
    for (size_t i = 0; i < words.size(); i += 3)
    {
      auto const row = words[i], col = words[i+1];
      double const value = std::bit_cast<double>(words[i+2]);
      std::pair<std::uint64_t,std::uint64_t> const key{row,col};
      if (row >= static_cast<std::uint64_t>(matrix.N()) ||
          col >= static_cast<std::uint64_t>(matrix.M()) || !std::isfinite(value) ||
          (!first && key <= previous))
        throw std::runtime_error("Invalid index, order or value in material matrix message");
      // NumaBCRSMatrix throws on a write outside its existing pattern.
      matrix[static_cast<size_t>(row)][static_cast<size_t>(col)][0][0] = value;
      first = false;
      previous = key;
    }
  }

  // Only owners hold assembled matrices on entry. Every rank receives the
  // same replicas required by existing diagonal weights, interface mass
  // redistribution, local-space extraction and diagnostics. Broadcast one
  // subdomain at a time to bound temporary memory; large payloads use chunks
  // rather than exceeding MPI's int count. No floating-point MPI reduction.
  template <class Matrix>
  void replicateMaterialMatrices(std::vector<Matrix>& matrices, Matrix const& pattern,
                                  Kaskade::BDDC::MpiDomainLayout const& layout,
                                  int chunkWords = std::numeric_limits<int>::max())
  {
    if (chunkWords <= 0)
      throw std::runtime_error("Material matrix broadcast chunk must be positive");
#ifdef KASKADE_HAVE_MPI
    if (layout.enabled() && layout.size() > 1)
    {
      for (size_t subdomain = 0; subdomain < matrices.size(); ++subdomain)
      {
        int const owner = layout.owner(static_cast<int>(subdomain));
        std::vector<std::uint64_t> words;
        std::exception_ptr error;
        try
        {
          if (layout.rank() == owner)
          {
            if (matrices[subdomain].N() != pattern.N() || matrices[subdomain].M() != pattern.M())
              throw std::runtime_error("Owned material matrix dimensions disagree with global pattern");
            words = packMaterialMatrix(matrices[subdomain]);
          }
        }
        catch (...) { error = std::current_exception(); }
        checkInterfaceAssemblyError(error,layout,"Material matrix packing");

        std::uint64_t length = words.size();
        MPI_Bcast(&length,1,MPI_UINT64_T,owner,MPI_COMM_WORLD);
        try
        {
          if (length%3 != 0 || length > words.max_size())
            throw std::runtime_error("Material matrix payload exceeds vector capacity");
          if (layout.rank() != owner)
            words.resize(static_cast<size_t>(length));
        }
        catch (...) { error = std::current_exception(); }
        checkInterfaceAssemblyError(error,layout,"Material matrix receive allocation");
        for (size_t offset = 0; offset < words.size();)
        {
          int const count = static_cast<int>(std::min(words.size()-offset,
            static_cast<size_t>(chunkWords)));
          MPI_Bcast(words.data()+offset,count,MPI_UINT64_T,owner,MPI_COMM_WORLD);
          offset += count;
        }
        try
        {
          if (layout.rank() != owner)
          {
            matrices[subdomain] = pattern;
            restoreMaterialMatrix(words,matrices[subdomain]);
          }
        }
        catch (...) { error = std::current_exception(); }
        checkInterfaceAssemblyError(error,layout,"Material matrix reconstruction");
      }
    }
#else
    (void)matrices;
    (void)pattern;
    (void)layout;
#endif
  }

  // Compare full sparse operators by a sorted merge, including structural
  // zeros and entries missing from either pattern. Avoid building another
  // map of all global entries for the production-size equivalence audit.
  template <class Matrix>
  void compareMaterialMatrices(Matrix const& computed, Matrix const& reference,
                                InterfaceBlockDifference& difference)
  {
    if (computed.N() != reference.N() || computed.M() != reference.M())
      throw std::runtime_error("Material matrix audit dimensions disagree");
    for (size_t row = 0; row < static_cast<size_t>(reference.N()); ++row)
    {
      auto const arow = computed[row], brow = reference[row];
      auto a = arow.begin(), b = brow.begin();
      while (a != arow.end() || b != brow.end())
      {
        double x = 0.0, y = 0.0;
        if (b == brow.end() || (a != arow.end() && a.index() < b.index()))
        {
          x = (*a)[0][0];
          ++a;
        }
        else if (a == arow.end() || b.index() < a.index())
        {
          y = (*b)[0][0];
          ++b;
        }
        else
        {
          x = (*a)[0][0];
          y = (*b)[0][0];
          ++a;
          ++b;
        }
        if (!std::isfinite(x) || !std::isfinite(y))
          throw std::runtime_error("Nonfinite material matrix audit entry");
        double const error = x-y;
        difference.maxAbs = std::max(difference.maxAbs,std::abs(error));
        difference.maxReference = std::max(difference.maxReference,std::abs(y));
        difference.squaredError += error*error;
        difference.squaredReference += y*y;
        ++difference.entries;
      }
    }
  }
}

#endif
