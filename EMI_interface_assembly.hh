#ifndef EMI_INTERFACE_ASSEMBLY_HH
#define EMI_INTERFACE_ASSEMBLY_HH

#include <exception>
#include <map>
#include <vector>

#include <boost/fusion/include/at_c.hpp>
#include <dune/grid/common/rangegenerators.hh>

#include "fem/assemblerFlags.hh"
#include "fem/functionspace.hh"
#include "fem/quadrature.hh"

#include "EMI_interface_exchange.hh"

namespace EmiBddc
{
  /**
   * Integrate all directed EMI mass interfaces in a single cell traversal.
   *
   * Uses the same Kaskade face quadrature, mapped/combined basis evaluators,
   * and SemiLinearization inner-boundary d2 as the legacy assembler. Each
   * directed face supplies inside/inside and inside/outside blocks; the
   * opposite direction is visited independently with its own face mapping.
   * The global lower triangle is retained and mirrored, matching the EMI
   * assembler's symmetric storage even for each directed pair separately.
   * This preserves orientation on general grids and higher FE orders.
   *
   * The caller selects pure mass, tau=0, and disables pair filtering. There
   * are no volume or outer-boundary mass terms in EMIModel. Physical BCs and
   * model coefficients are not changed. Face-wise accumulation can change
   * roundoff relative to the legacy cell/pair accumulation order.
   */
  template <class Grid, class Spaces, class Material, class Linearization, class Model>
  InterfaceBlocks assembleInterfaceMassBlocks(Grid const& grid, Spaces const& spaces,
                                              Material const& material,
                                              Linearization const& functional,
                                              Model const& model,
                                              std::map<int,size_t> const& tagToSubdomain,
                                              Kaskade::BDDC::MpiDomainLayout const& layout,
                                              bool distributed)
  {
    InterfaceBlocks blocks;
    std::exception_ptr error;
    try
    {
      constexpr int dim = Grid::dimension;
      using Scalar = typename Linearization::Scalar;
      Kaskade::ShapeFunctionCache<Grid,Scalar> sfCache;
      auto evaluators = Kaskade::getEvaluators(spaces,&sfCache);
      auto neighbors = evaluators;
      auto cache = functional.createInnerBoundaryCache(Kaskade::Assembler::MATRIX);
      Kaskade::QuadratureTraits<Dune::QuadratureRule<typename Grid::ctype,dim-1>> quadrature;
      auto const view = grid.leafGridView();
      Dune::FieldVector<double,dim> zero(0.0);
      for (auto const& cell : Dune::elements(view))
      {
        int const tag = material.value(cell,zero);
        auto const source = tagToSubdomain.find(tag);
        if (source == tagToSubdomain.end() ||
            (distributed && !layout.owns(static_cast<int>(source->second))))
          continue;
        bool moved = false;
        for (auto face = view.ibegin(cell); face != view.iend(cell); ++face)
        {
          if (face->boundary() || !face->neighbor())
            continue;
          int const neighborTag = material.value(face->outside(),zero);
          auto const destination = tagToSubdomain.find(neighborTag);
          if (tag == neighborTag || destination == tagToSubdomain.end() ||
              (model.isExtracellular(tag) && model.isExtracellular(neighborTag)))
            continue;
          if (!functional.considerFace(*face))
            continue;
          if (!moved)
          {
            Kaskade::moveEvaluatorsToCell(evaluators,cell);
            moved = true;
          }
          Kaskade::moveEvaluatorsToCell(neighbors,face->outside());
          cache.moveTo(face);
          auto const& inside = boost::fusion::at_c<0>(evaluators);
          auto const& outside = boost::fusion::at_c<0>(neighbors);
          auto const insideIds = inside.globalIndices(), outsideIds = outside.globalIndices();
          size_t const n = insideIds.size(), m = outsideIds.size();
          std::vector<double> diagonal(n*n,0.0), coupling(n*m,0.0);
          auto const& rule = quadrature.rule(face->type(),
            functional.integrationOrder(cell,Kaskade::maxOrder(evaluators),true));
          for (auto const& point : rule)
          {
            auto const position = point.position();
            auto const weight = point.weight()*face->geometry().integrationElement(position);
            Kaskade::moveEvaluatorsToIntegrationPoint(evaluators,face->geometryInInside().global(position));
            Kaskade::moveEvaluatorsToIntegrationPoint(neighbors,face->geometryInOutside().global(position));
            cache.evaluateAt(position,evaluators,neighbors);
            for (size_t i = 0; i < n; ++i)
            {
              // Match the assembler's lower-triangle-then-mirror convention.
              for (size_t j = 0; j <= i; ++j)
                diagonal[i*n+j] += weight*cache.template d2<0,0>(
                  inside.evalData[i],inside.evalData[j],true)[0][0];
              for (size_t j = 0; j < m; ++j)
                coupling[i*m+j] += weight*cache.template d2<0,0>(
                  inside.evalData[i],outside.evalData[j],false)[0][0];
            }
          }
          for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < i; ++j)
              diagonal[j*n+i] = diagonal[i*n+j];
          accumulateInterfaceFace(blocks[{source->second,destination->second}],
                                   insideIds,outsideIds,diagonal,coupling);
        }
      }
      mirrorInterfaceLowerTriangles(blocks);
    }
    catch (...) { error = std::current_exception(); }
    checkInterfaceAssemblyError(error,layout);
    return blocks;
  }
}

#endif
