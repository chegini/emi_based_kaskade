#ifndef EMI_MODEL_HH
#define EMI_MODEL_HH

#include <algorithm>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"
#include "fem/gridBasics.hh"
#include "fem/variables.hh"
#include "utilities/linalg/scalarproducts.hh"

struct HashPair
{
  template <class T1, class T2>
  size_t operator()(std::pair<T1,T2> const& p) const
  {
    auto const h1 = std::hash<T1>{}(p.first);
    auto const h2 = std::hash<T2>{}(p.second);
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
  }
};

template <class Scalar_, class VarSet, class Material, class Grid_, class Spaces, class MembraneModel>
class EMIModel : public Kaskade::FunctionalBase<Kaskade::WeakFormulation>
{
  using Self = EMIModel<Scalar_,VarSet,Material,Grid_,Spaces,MembraneModel>;

public:
  using Scalar = Scalar_;
  using AnsatzVars = VarSet;
  using OriginVars = VarSet;
  using TestVars = VarSet;
  using GridView = typename AnsatzVars::GridView;
  using Membrane = MembraneModel;

  static constexpr int dim = AnsatzVars::Grid::dimension;
  static constexpr int uIdx = 0;
  static constexpr int uSpaceIdx = Kaskade::spaceIndex<AnsatzVars,uIdx>;

  template <int row>
  using TestComponents = std::integral_constant<size_t,TestVars::template Components<row>::m>;

  template <int row>
  using AnsatzComponents = std::integral_constant<size_t,AnsatzVars::template Components<row>::m>;

  class DomainCache
  {
  public:
    DomainCache(Self const& F_, typename AnsatzVars::VariableSet const& vars_, int = 7)
    : F(F_), vars(vars_)
    {}

    template <class Cell>
    void moveTo(Cell const& cell)
    {
      Dune::FieldVector<double,dim> zero(0);
      cellMaterial = F.material.value(cell,zero);
      extracellular = std::find(F.extraTags.begin(), F.extraTags.end(), cellMaterial) != F.extraTags.end();
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const&, Evaluators const& evaluators)
    {
      using namespace boost::fusion;
      u = component<uIdx>(vars).value(at_c<uSpaceIdx>(evaluators));
      du = component<uIdx>(vars).derivative(at_c<uSpaceIdx>(evaluators));
      f = 0.0;
    }

    Scalar d0() const { return 0.0; }

    template <int row>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m>
    d1(Kaskade::VariationalArg<Scalar,dim,TestComponents<row>::value> const& argT) const
    {
      auto const sigma = extracellular ? F.sigmaE : F.sigmaI;
      return -sigma*sp(du,argT.derivative) + f*argT.value;
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<col>::m>
    d2(Kaskade::VariationalArg<Scalar,dim,TestComponents<row>::value> const& argT,
       Kaskade::VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& argA) const
    {
      if constexpr (row != col)
        return 0.0;

      Scalar p = -1.0;
      if (F.massMode == 1)
        p = 0.0;
      if (F.massMode > 1)
        p = 1.0;

      auto const sigma = extracellular ? F.sigmaE : F.sigmaI;
      return -p*sigma*sp(argT.derivative,argA.derivative);
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<col>::m>
    b2(Kaskade::VariationalArg<Scalar,dim> const&,
       Kaskade::VariationalArg<Scalar,dim> const&) const
    {
      return 0.0;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u, f;
    Dune::FieldMatrix<Scalar,AnsatzComponents<uIdx>::value,dim> du;
    Kaskade::LinAlg::EuclideanScalarProduct sp;
    int cellMaterial = 0;
    bool extracellular = false;
  };

  class BoundaryCache
  {
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

  public:
    BoundaryCache(Self const& F_, typename AnsatzVars::VariableSet const& vars_, int = 7)
    : F(F_), vars(vars_)
    {}

    void moveTo(FaceIterator const& face)
    {
      Dune::FieldVector<double,dim> zero(0);
      cellMaterial = F.material.value((*face).inside(),zero);
      extracellular = std::find(F.extraTags.begin(), F.extraTags.end(), cellMaterial) != F.extraTags.end();
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const&,
                    Evaluators const& evaluators)
    {
      using namespace boost::fusion;
      ue = component<uIdx>(vars).value(at_c<uSpaceIdx>(evaluators));
    }

    template <int row>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m>
    d1(Kaskade::VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if (extracellular)
        return -gamma*ue*arg.value;
      return 0.0;
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<col>::m>
    d2(Kaskade::VariationalArg<Scalar,dim,TestComponents<row>::value> const& argT,
       Kaskade::VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& argA) const
    {
      if (F.massMode == 1 || !extracellular)
        return 0.0;

      Scalar p = -1.0;
      if (F.massMode > 1)
        p = 1.0;

      return -p*gamma*argT.value*argA.value;
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<col>::m>
    b2(Kaskade::VariationalArg<Scalar,dim> const&,
       Kaskade::VariationalArg<Scalar,dim> const&) const
    {
      return 0.0;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    Scalar gamma = 1e-5;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> ue = 0.0;
    int cellMaterial = 0;
    bool extracellular = false;
  };

  class InnerBoundaryCache
  {
  public:
    InnerBoundaryCache(Self const& F_, typename AnsatzVars::VariableSet const& vars_, int = 7)
    : F(F_), vars(vars_), C_m(F_.C_m), R(F_.R)
    {}

    template <class FaceIterator>
    void moveTo(FaceIterator const& f)
    {
      face = *f;
      Dune::FieldVector<double,dim> zero(0);
      cellDomain = F.material.value(face.inside(),zero);
      neighbourDomain = F.material.value(face.outside(),zero);
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators, Evaluators const& neighbourEvaluators)
    {
      if (cellDomain == neighbourDomain)
        return;

      Scalar const cellU = vars.template value<uIdx>(evaluators);
      Scalar const neighbourU = vars.template value<uIdx>(neighbourEvaluators);
      Scalar const v = cellU - neighbourU;

      bool const cellExtra = F.isExtracellular(cellDomain);
      bool const neighbourExtra = F.isExtracellular(neighbourDomain);
      twoExtra = cellExtra && neighbourExtra;
      membraneInterface = cellExtra != neighbourExtra;

      if (twoExtra)
        std::tie(current,dcurrentCell,dcurrentNeighbour) = std::make_tuple(0.0,0.0,0.0);
      else if (cellDomain > neighbourDomain)
        std::tie(current,dcurrentCell,dcurrentNeighbour) = ionCurrent(v);
      else
      {
        std::tie(current,dcurrentCell,dcurrentNeighbour) = ionCurrent(-v);
        current = -current;
      }
    }

    std::tuple<Scalar,Scalar,Scalar> ionCurrent(Scalar v) const
    {
      if (membraneInterface)
      {
        Scalar const ionic = -F.membrane().current(v,0);
        return std::make_tuple(ionic,0.0,0.0);
      }

      return std::make_tuple(v/R,0.0,0.0);
    }

    template <int row>
    Dune::FieldVector<Scalar,1>
    d1(Kaskade::VariationalArg<Scalar,dim> const& argT) const
    {
      // Inner-boundary terms are assembled only on material interfaces.
      // If both sides have the same material tag, this is an ordinary internal
      // face and must not contribute a membrane or gap-junction current.
      if (cellDomain == neighbourDomain)
        return 0.0;
      return -argT.value*current;
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,1,1>
    d2(Kaskade::VariationalArg<Scalar,dim> const& argT,
       Kaskade::VariationalArg<Scalar,dim> const& argA,
       bool centerCell) const
    {
      // Same-material faces are not EMI interfaces. In pure mass assembly
      // mode the stiffness/current Jacobian is intentionally suppressed.
      if (cellDomain == neighbourDomain || F.massMode == 1)
        return 0.0;

      Scalar const dCurrent = centerCell ? dcurrentCell : dcurrentNeighbour;
      Scalar p = -1.0;
      if (F.massMode > 1)
        p = 1.0;

      return -p*argT.value*dCurrent*argA.value;
    }

    template <int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<col>::m>
    b2(Kaskade::VariationalArg<Scalar,dim> const& argT,
       Kaskade::VariationalArg<Scalar,dim> const& argA,
       bool centerCell) const
    {
      // The membrane capacitance couples intracellular/extracellular sides
      // only. Faces between two extracellular regions have interface
      // conductivity but no membrane capacitance term.
      if (cellDomain == neighbourDomain || twoExtra || F.massMode == 0)
        return 0.0;

      Scalar const sign = centerCell ? 1.0 : -1.0;
      return C_m*sign*argT.value*argA.value;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    Scalar C_m;
    Scalar R;
    int cellDomain = 0;
    int neighbourDomain = 0;
    Scalar current = 0.0;
    Scalar dcurrentCell = 0.0;
    Scalar dcurrentNeighbour = 0.0;
    Kaskade::Face<GridView> face;
    bool membraneInterface = false;
    bool twoExtra = false;
  };

  template <int row, int col>
  struct D2: public Kaskade::FunctionalBase<Kaskade::VariationalFunctional>::D2<row,col>
  {
    static constexpr bool present = true;
    static constexpr bool symmetric = true;
  };

  template <int row, int col>
  struct B2: public Kaskade::FunctionalBase<Kaskade::VariationalFunctional>::B2<row,col>
  {
    static constexpr bool present = true;
    static constexpr bool symmetric = true;
    static constexpr bool constant = false;
  };

  EMIModel(Material const& material_,
           Grid_ const& grid_,
           Spaces const& spaces_,
           double penalty_ = 1e6,
           double sigmaI_ = 3.0,
           double sigmaE_ = 20.0,
           double C_m_ = 1.0,
           double R_ = 0.1,
           double RExtra_ = 1e-7)
  : penalty(penalty_),
    sigmaI(sigmaI_),
    sigmaE(sigmaE_),
    C_m(C_m_),
    R(R_),
    RExtra(RExtra_),
    material(material_),
    grid(grid_),
    spaces(spaces_)
  {}

  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& us, typename AnsatzVars::VariableSet& u) const
  {
    Kaskade::interpolateGloballyWeak<Kaskade::Volume>(boost::fusion::at_c<row>(u.data),
                                                      ScaledFunction<WeakFunctionView>(row == uIdx,us,*this));
  }

  template <class WeakFunctionView>
  struct ScaledFunction
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static constexpr int components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& us_,
                   Self const& f_)
    : doScaling(doScaling_), us(us_), f(f_)
    {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
    {
      ValueType value = 0.0;
      if (doScaling)
        value = us.value(cell,localCoordinate);
      return value;
    }

  private:
    bool doScaling;
    WeakFunctionView const& us;
    Self const& f;
  };

  bool considerFace(Kaskade::Face<GridView> const& face) const
  {
    Dune::FieldVector<Scalar,dim> zero(0.0);
    return material.value(face.inside(),zero) != material.value(face.outside(),zero);
  }

  bool isExtracellular(int tag) const
  {
    return std::find(extraTags.begin(), extraTags.end(), tag) != extraTags.end();
  }

  Scalar time() const { return t; }
  void time(Scalar tnew) { t = tnew; }

  void Mass_stiff(int mass) { massMode = mass; }

  Membrane const& membrane() const { return memb; }

  void extracellularMaterials(std::vector<int> const& tags)
  {
    extraTags = tags;
  }

  template <class Cell>
  int integrationOrder(Cell const&, int shapeFunctionOrder, bool boundary) const
  {
    return boundary ? 2*shapeFunctionOrder : 2*shapeFunctionOrder - 2;
  }

private:
  double penalty;
  double sigmaI;
  double sigmaE;
  double C_m;
  double R;
  double RExtra;
  double t = 0.0;
  Material const& material;
  Grid_ const& grid;
  Spaces const& spaces;
  MembraneModel memb;
  int massMode = 10;
  std::vector<int> extraTags;
};

#endif
