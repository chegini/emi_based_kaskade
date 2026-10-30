#ifndef EMI_INTERFACE_BLOCKS_HH
#define EMI_INTERFACE_BLOCKS_HH

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace EmiBddc
{
  using DirectedInterface = std::pair<size_t,size_t>;
  using InterfaceEntry = std::pair<size_t,size_t>;
  using InterfaceBlock = std::map<InterfaceEntry,double>;
  using InterfaceBlocks = std::map<DirectedInterface,InterfaceBlock>;

  // A directed face contributes its inside/inside and inside/outside blocks.
  // The evaluator indices include the mapper's constraint transformation.
  // The EMI assembler declares this scalar block symmetric: scatter retains
  // ONLY global row >= column, and get<Matrix>(false) mirrors it afterward.
  // Preserve that convention even for an individually directed pair block.
  // Accumulate complete face integrals before changing material operators:
  // the legacy split transfers half of an ASSEMBLED directed pair block.
  template <class Indices>
  void accumulateInterfaceFace(InterfaceBlock& block,
                               Indices const& inside, Indices const& outside,
                               std::vector<double> const& diagonal,
                               std::vector<double> const& coupling)
  {
    size_t const n = inside.size(), m = outside.size();
    if (diagonal.size() != n*n || coupling.size() != n*m)
      throw std::runtime_error("Interface face matrix dimensions disagree with evaluator indices");
    for (size_t i = 0; i < n; ++i)
    {
      for (size_t j = 0; j < n; ++j)
        if (inside[i] >= inside[j] && diagonal[i*n+j] != 0.0)
          block[{inside[i],inside[j]}] += diagonal[i*n+j];
      for (size_t j = 0; j < m; ++j)
        if (inside[i] >= outside[j] && coupling[i*m+j] != 0.0)
          block[{inside[i],outside[j]}] += coupling[i*m+j];
    }
  }

  inline void mirrorInterfaceLowerTriangles(InterfaceBlocks& blocks)
  {
    for (auto& [pair,lower] : blocks)
    {
      (void)pair;
      InterfaceBlock full;
      for (auto const& [entry,value] : lower)
      {
        if (entry.first < entry.second)
          throw std::runtime_error("Interface block expected lower-triangle storage");
        if (!std::isfinite(value))
          throw std::runtime_error("Nonfinite integrated interface mass entry");
        if (value != 0.0)
        {
          full[entry] = value;
          if (entry.first != entry.second)
            full[{entry.second,entry.first}] = value;
        }
      }
      lower = std::move(full);
    }
  }

  // Five uint64 words per nonzero: source subdomain, destination subdomain,
  // global row, global column, and the exact IEEE double bits. No floating
  // reduction is used during exchange. Explicit words avoid struct padding.
  inline std::vector<std::uint64_t> packInterfaceBlocks(InterfaceBlocks const& blocks)
  {
    static_assert(sizeof(double) == sizeof(std::uint64_t));
    static_assert(sizeof(size_t) <= sizeof(std::uint64_t));
    std::vector<std::uint64_t> words;
    for (auto const& [pair,block] : blocks)
      for (auto const& [entry,value] : block)
        if (value != 0.0)
        {
          if (!std::isfinite(value))
            throw std::runtime_error("Nonfinite interface mass entry");
          words.insert(words.end(),{pair.first,pair.second,entry.first,entry.second,
                                    std::bit_cast<std::uint64_t>(value)});
        }
    return words;
  }

  template <class Owner>
  void unpackInterfaceBlocks(std::vector<std::uint64_t> const& words,
                             size_t begin, size_t end, int sender,
                             size_t subdomains, size_t dofs, Owner owner,
                             InterfaceBlocks& blocks)
  {
    if (begin > end || end > words.size() || (end-begin)%5 != 0)
      throw std::runtime_error("Invalid sparse interface message length");
    for (size_t i = begin; i < end; i += 5)
    {
      auto const source = words[i], destination = words[i+1];
      auto const row = words[i+2], col = words[i+3];
      double const value = std::bit_cast<double>(words[i+4]);
      if (source >= subdomains || destination >= subdomains || source == destination ||
          row >= dofs || col >= dofs || !std::isfinite(value) ||
          owner(static_cast<int>(source)) != sender)
        throw std::runtime_error("Invalid source ownership or entry in sparse interface message");
      auto& block = blocks[{static_cast<size_t>(source),static_cast<size_t>(destination)}];
      if (!block.emplace(InterfaceEntry{row,col},value).second)
        throw std::runtime_error("Duplicate sparse interface entry from MPI owners");
    }
  }

  struct InterfaceBlockDifference
  {
    double maxAbs = 0.0, maxReference = 0.0;
    double squaredError = 0.0, squaredReference = 0.0;
    size_t entries = 0;

    void compare(InterfaceBlock const& computed, InterfaceBlock const& reference)
    {
      auto a = computed.begin(), b = reference.begin();
      while (a != computed.end() || b != reference.end())
      {
        double x = 0.0, y = 0.0;
        if (b == reference.end() || (a != computed.end() && a->first < b->first))
          x = (a++)->second;
        else if (a == computed.end() || b->first < a->first)
          y = (b++)->second;
        else
        {
          x = (a++)->second;
          y = (b++)->second;
        }
        if (!std::isfinite(x) || !std::isfinite(y))
          throw std::runtime_error("Nonfinite interface entry in matrix equivalence audit");
        double const difference = x-y;
        maxAbs = std::max(maxAbs,std::abs(difference));
        maxReference = std::max(maxReference,std::abs(y));
        squaredError += difference*difference;
        squaredReference += y*y;
        ++entries;
      }
    }

    bool acceptable() const { return maxAbs <= 1e-14 + 1e-11*maxReference; }
  };

  // Deterministic legacy pair order and exact conservation: -M/2 in the
  // source material and +M/2 in its neighbor, including off-diagonal entries.
  template <class AddEntry>
  void redistributeInterfaceBlocks(InterfaceBlocks const& blocks, AddEntry add)
  {
    // Sort by undirected pair first, then directed source, as in the old loop.
    std::map<std::pair<DirectedInterface,size_t>,InterfaceBlock const*> ordered;
    for (auto const& [pair,block] : blocks)
      ordered[{{std::min(pair.first,pair.second),std::max(pair.first,pair.second)},pair.first}] = &block;
    for (auto const& [key,block] : ordered)
    {
      size_t const source = key.second;
      size_t const destination = source == key.first.first ? key.first.second : key.first.first;
      for (auto const& [entry,value] : *block)
        if (value != 0.0)
        {
          double const half = 0.5*value;
          add(source,entry.first,entry.second,-half);
          add(destination,entry.first,entry.second,half);
        }
    }
  }
}

#endif
