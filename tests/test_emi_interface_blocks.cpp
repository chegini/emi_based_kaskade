#include "../EMI_interface_blocks.hh"

#include <array>
#include <iostream>
#include <string>

using namespace EmiBddc;

void require(bool condition, char const* message)
{
  if (!condition)
    throw std::runtime_error(message);
}

template <class F>
void rejects(F action)
{
  bool failed = false;
  try { action(); }
  catch (std::runtime_error const&) { failed = true; }
  require(failed,"Malformed interface data was accepted");
}

int main()
{
  // P1 triangle face: area=3, C_m=2, integral phi_i phi_j is
  // C_m*area*(1+delta_ij)/12. Shuffle global numbering on both traces.
  std::vector<size_t> inside{7,0,5}, outside{2,6,1};
  std::vector<double> diagonal{1,.5,.5, .5,1,.5, .5,.5,1};
  std::vector<double> coupling{-1,-.5,-.5, -.5,-1,-.5, -.5,-.5,-1};
  InterfaceBlocks blocks;
  accumulateInterfaceFace(blocks[{0,1}],inside,outside,diagonal,coupling);
  accumulateInterfaceFace(blocks[{1,0}],outside,inside,diagonal,coupling);
  // A second face of the same pair must accumulate, not overwrite entries.
  accumulateInterfaceFace(blocks[{0,1}],inside,outside,diagonal,coupling);
  mirrorInterfaceLowerTriangles(blocks);

  // Independent reference: scatter full directed face matrices into dense
  // matrices, then extract/mirror their global lower triangles.
  std::array<std::array<double,64>,2> reference{};
  for (size_t direction = 0; direction < 2; ++direction)
  {
    auto const& rows = direction == 0 ? inside : outside;
    auto const& cols = direction == 0 ? outside : inside;
    double const repeats = direction == 0 ? 2.0 : 1.0;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
      {
        reference[direction][rows[i]*8+rows[j]] += repeats*(i == j ? 1.0 : .5);
        reference[direction][rows[i]*8+cols[j]] -= repeats*(i == j ? 1.0 : .5);
      }
    InterfaceBlock sparse;
    for (size_t i = 0; i < 8; ++i)
      for (size_t j = 0; j <= i; ++j)
        if (reference[direction][i*8+j] != 0.0)
        {
          sparse[{i,j}] = reference[direction][i*8+j];
          sparse[{j,i}] = reference[direction][i*8+j];
        }
    require(blocks.at({direction,1-direction}) == sparse,
            "Shared face scatter differs from symmetric legacy pair extraction");
  }
  require(!blocks.at({0,1}).contains({0,6}),
          "Upper-global-triangle directed coupling must be discarded before mirroring");
  require(blocks.at({0,1}).at({2,7}) == -2.0,
          "Global lower-triangle coupling was not mirrored correctly");

  // Simulate owner ranks, including idle rank 0 (master-worker) and an
  // additional idle worker. Transport must preserve every double bit.
  for (int ranks : {1,2,3,8})
    for (bool master : {false,true})
    {
      if (master && ranks == 1)
        continue;
      auto owner = [=](int source) { return master ? 1+source%(ranks-1) : source%ranks; };
      InterfaceBlocks gathered;
      for (int rank = 0; rank < ranks; ++rank)
      {
        InterfaceBlocks local;
        for (auto const& [pair,block] : blocks)
          if (owner(pair.first) == rank)
            local[pair] = block;
        auto const message = packInterfaceBlocks(local);
        unpackInterfaceBlocks(message,0,message.size(),rank,2,8,owner,gathered);
      }
      require(gathered == blocks,"Owner-wise sparse exchange changed matrix entries");
    }

  // Every transferred entry cancels exactly in the embedded global sum.
  std::array<InterfaceBlock,2> delta;
  redistributeInterfaceBlocks(blocks,[&](size_t s,size_t r,size_t c,double value)
  {
    delta[s][{r,c}] += value;
  });
  InterfaceBlock sum = delta[0];
  for (auto const& [entry,value] : delta[1])
    sum[entry] += value;
  for (auto const& [entry,value] : sum)
  {
    (void)entry;
    require(value == 0.0,"Interface redistribution changed the global operator");
  }
  // Check signs and the exact half, not merely conservation.
  require(delta[0].at({7,7}) == -1.0 && delta[1].at({7,7}) == 1.0,
          "Incorrect source/destination half transfer");

  InterfaceBlockDifference difference;
  difference.compare(blocks.at({0,1}),blocks.at({0,1}));
  require(difference.acceptable() && difference.maxAbs == 0.0,"Exact comparison failed");
  auto corrupted = blocks.at({0,1});
  corrupted[{3,3}] = .25; // extra nonzero absent from reference
  difference.compare(corrupted,blocks.at({0,1}));
  require(!difference.acceptable(),"Audit missed an extra nonzero");
  rejects([&] { accumulateInterfaceFace(blocks[{0,1}],inside,outside,{},coupling); });
  auto message = packInterfaceBlocks(blocks);
  InterfaceBlocks decoded;
  auto owner = [](int source) { return source; };
  rejects([&] { unpackInterfaceBlocks(message,0,1,0,2,8,owner,decoded); });
  auto first = packInterfaceBlocks(InterfaceBlocks{{{0,1},blocks.at({0,1})}});
  rejects([&] { unpackInterfaceBlocks(first,0,first.size(),1,2,8,owner,decoded); });
  unpackInterfaceBlocks(first,0,first.size(),0,2,8,owner,decoded);
  rejects([&] { unpackInterfaceBlocks(first,0,first.size(),0,2,8,owner,decoded); });
  auto invalid = first;
  invalid[2] = 8;
  rejects([&] { InterfaceBlocks out; unpackInterfaceBlocks(invalid,0,invalid.size(),0,2,8,owner,out); });
  invalid = first;
  invalid[4] = std::bit_cast<std::uint64_t>(std::numeric_limits<double>::infinity());
  rejects([&] { InterfaceBlocks out; unpackInterfaceBlocks(invalid,0,invalid.size(),0,2,8,owner,out); });
  require(packInterfaceBlocks({}).empty(),"Empty owner sent data");
  std::cout << "EMI interface block checks passed\n";
}
