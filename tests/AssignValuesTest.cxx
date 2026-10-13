/*!
 * \file   tests/AssignValuesTest.cxx
 * \brief  This test checks the `assign_values` function: the errors it
 * reports, the copies between views defined on the same mesh, and the copies
 * between views defined on the same location of two different meshes (the
 * whole mesh and a submesh, or two submeshes). All data layouts are tested:
 * contiguous views and views with a stride.
 * \date   12/10/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <array>
#include <memory>
#include <vector>
#include <optional>
#include <filesystem>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

struct {
  const char* mesh_file = nullptr;
  int parallel = 0;
  int order = 1;
} parameters;

/*!
 * \return the integration rule used by most partial quadrature spaces
 * \param[in] e: finite element
 * \param[in] tr: element transformation
 *
 * \note this rule has more than one integration point for all the elements
 * of the test meshes.
 */
static const mfem::IntegrationRule& getIntegrationRule(
    const mfem::FiniteElement& e,
    const mfem::ElementTransformation& tr) noexcept {
  return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e) + 2);
}  // end of getIntegrationRule

/*!
 * \return an integration rule with less integration points than the one
 * returned by `getIntegrationRule`
 * \param[in] e: finite element
 * \param[in] tr: element transformation
 */
static const mfem::IntegrationRule& getOtherIntegrationRule(
    const mfem::FiniteElement& e,
    const mfem::ElementTransformation& tr) noexcept {
  return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e));
}  // end of getOtherIntegrationRule

struct AssignValuesTest final : public tfel::tests::TestCase {
  AssignValuesTest() : tfel::tests::TestCase("MFEMMGIS", "AssignValuesTest") {}
  tfel::tests::TestResult execute() override {
    if (parameters.parallel) {
#ifdef MFEM_USE_MPI
      this->test1<true>();
      this->test2();
      this->test3<true>();
#else  /* MFEM_USE_MPI */
      mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
      this->test2();
      this->test3<false>();
    }
    return this->result;
  }  // end of execute

 private:
  //! \brief data layout of a view
  struct Layout {
    //! \brief position of the first component in a block of values
    mfem_mgis::size_type begin;
    //! \brief size of a block of values
    mfem_mgis::size_type stride;
  };
  //! \brief value of the data which shall not be modified
  static constexpr auto sentinel = mfem_mgis::real{-1000};
  //! \brief tolerance
  static constexpr auto eps = mfem_mgis::real{1e-12};
  /*!
   * \return a finite element discretization built on the given mesh file
   * \param[in, out] ctx: execution context
   */
  static mfem_mgis::FiniteElementDiscretization makeFiniteElementDiscretization(
      mfem_mgis::Context& ctx) {
    return mfem_mgis::FiniteElementDiscretization{
        ctx,
        {{"MeshFileName", parameters.mesh_file},
         {"FiniteElementFamily", "H1"},
         {"FiniteElementOrder", parameters.order},
         {"UnknownsSize", 1},
         {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
         {"Parallel", bool(parameters.parallel)}}};
  }  // end of makeFiniteElementDiscretization
  //! \brief errors reported by `assign_values`
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    const auto m = MaterialIdentifier{.id = 1};
    const auto qspace = std::make_shared<const PartialQuadratureSpace>(
        fed, m, &getIntegrationRule);
    auto f = PartialQuadratureFunction(qspace, 1);
    // the same material of another mesh discretization
    auto fed2 = makeFiniteElementDiscretization(ctx);
    auto v1 = PartialQuadratureFunction(
        std::make_shared<const PartialQuadratureSpace>(fed2, m,
                                                       &getIntegrationRule),
        1);
    TFEL_TESTS_CHECK(!assign_values(ctx, f, v1));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "unmatched mesh discretization");
    // another material
    auto v2 = PartialQuadratureFunction(
        std::make_shared<const PartialQuadratureSpace>(
            fed, MaterialIdentifier{.id = 2}, &getIntegrationRule),
        1);
    TFEL_TESTS_CHECK(!assign_values(ctx, f, v2));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "unmatched location");
    // another number of components
    auto v3 = PartialQuadratureFunction(qspace, 2);
    TFEL_TESTS_CHECK(!assign_values(ctx, f, v3));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "unmatched number of components");
    // same mesh, other integration rule. In parallel, a process on which the
    // material is not defined has nothing to copy.
    auto v4 = PartialQuadratureFunction(
        std::make_shared<const PartialQuadratureSpace>(
            fed, m, &getOtherIntegrationRule),
        1);
    TFEL_TESTS_CHECK(!assign_values(ctx, f, v4));
    // submesh of the material, other integration rule
    const auto submesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(submesh));
    const auto fespace =
        fed.getFiniteElementSpacesManager().getFiniteElementSpace<parallel>(
            ctx, *submesh, 1);
    TFEL_TESTS_ASSERT(isValid(fespace));
    auto v5 = PartialQuadratureFunction(
        std::make_shared<const PartialQuadratureSpace>(
            fed, *fespace, m, &getOtherIntegrationRule),
        1);
    TFEL_TESTS_CHECK(!assign_values(ctx, f, v5));
  }  // end of test1
  //! \brief copies between views defined on the same mesh
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    const auto m = MaterialIdentifier{.id = 1};
    const auto b = BoundaryIdentifier{.id = 1};
    const auto mqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, m, &getIntegrationRule);
    // another partial quadrature space on the same mesh
    const auto mqspace2 = std::make_shared<const PartialQuadratureSpace>(
        fed, m, &getIntegrationRule);
    const auto bqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, b, &getIntegrationRule);
    this->checkAssignments(ctx, mqspace, mqspace, true);
    this->checkAssignments(ctx, mqspace, mqspace2, true);
    this->checkAssignments(ctx, bqspace, bqspace, true);
  }  // end of test2
  //! \brief copies between views defined on different meshes
  template <bool parallel>
  void test3() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    // submeshes and finite element spaces
    const auto msubmesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    const auto bsubmesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_BOUNDARIES);
    const auto bsubmesh2 =
        fed.getSubMesh<parallel>(ctx, Parameter{std::vector<Parameter>{1, 3}},
                                 MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(msubmesh));
    TFEL_TESTS_ASSERT(isValid(bsubmesh));
    TFEL_TESTS_ASSERT(isValid(bsubmesh2));
    auto manager = fed.getFiniteElementSpacesManager();
    const auto mfespace =
        manager.getFiniteElementSpace<parallel>(ctx, *msubmesh, 1);
    const auto bfespace =
        manager.getFiniteElementSpace<parallel>(ctx, *bsubmesh, 1);
    const auto bfespace2 =
        manager.getFiniteElementSpace<parallel>(ctx, *bsubmesh2, 1);
    TFEL_TESTS_ASSERT(isValid(mfespace));
    TFEL_TESTS_ASSERT(isValid(bfespace));
    TFEL_TESTS_ASSERT(isValid(bfespace2));
    // partial quadrature spaces
    const auto m = MaterialIdentifier{.id = 1};
    const auto b = BoundaryIdentifier{.id = 1};
    // material on the whole mesh and on its submesh
    const auto mqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, m, &getIntegrationRule);
    const auto mqspace2 = std::make_shared<const PartialQuadratureSpace>(
        fed, *mfespace, m, &getIntegrationRule);
    // boundary on the whole mesh, on its submesh, on the submesh of the
    // material and on the submesh of the boundaries 1 and 3
    const auto bqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, b, &getIntegrationRule);
    const auto bqspace2 = std::make_shared<const PartialQuadratureSpace>(
        fed, *bfespace, b, &getIntegrationRule);
    const auto bqspace3 = std::make_shared<const PartialQuadratureSpace>(
        fed, *mfespace, b, &getIntegrationRule);
    const auto bqspace4 = std::make_shared<const PartialQuadratureSpace>(
        fed, *bfespace2, b, &getIntegrationRule);
    this->checkAssignments(ctx, mqspace, mqspace2, true);
    this->checkAssignments(ctx, bqspace, bqspace2, true);
    this->checkAssignments(ctx, bqspace, bqspace4, true);
    this->checkAssignments(ctx, bqspace2, bqspace4, true);
    // boundary seen through the submesh of a material. In 3D, some boundary
    // elements of this submesh are not oriented as in the whole mesh, so the
    // copy shall either fail or be exact. In 2D, the boundary 1 of the test
    // mesh is exterior and inside the material 1, so the copy shall be exact.
    const auto exact = fed.getMesh<parallel>().Dimension() == 2;
    this->checkAssignments(ctx, bqspace, bqspace3, exact);
    this->checkAssignments(ctx, bqspace2, bqspace3, exact);
  }  // end of test3
  /*!
   * \return the value of a linear field
   * \param[in] x: position
   */
  static mfem_mgis::real getLinearField(const mfem::Vector& x) noexcept {
    return 1 + x[0] + 2 * x[1] + (x.Size() > 2 ? 3 * x[2] : mfem_mgis::real{0});
  }  // end of getLinearField
  /*!
   * \return the values of the linear field at the integration points of a
   * partial quadrature space
   * \param[in] s: partial quadrature space
   */
  static std::shared_ptr<mfem_mgis::PartialQuadratureFunction>
  evaluateLinearField(
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& s) {
    return mfem_mgis::PartialQuadratureFunction::evaluate(
        s, [](const mfem::FiniteElement&, mfem::ElementTransformation& tr) {
          auto x = mfem::Vector{};
          tr.Transform(tr.GetIntPoint(), x);
          return getLinearField(x);
        });
  }  // end of evaluateLinearField
  /*!
   * \return the maximum error of the copy of the values of a view into
   * another view, or nothing if the copy fails
   *
   * The values of the source view are the values of the linear field shifted
   * by ten times the index of the component. The values of the blocks which
   * are not part of the views are set to `sentinel` and shall not be
   * modified.
   *
   * \param[in, out] ctx: execution context
   * \param[in] fs: partial quadrature space of the destination view
   * \param[in] fl: layout of the destination view
   * \param[in] vs: partial quadrature space of the source view
   * \param[in] vl: layout of the source view
   * \param[in] nc: number of components of the views
   */
  static std::optional<mfem_mgis::real> getAssignmentError(
      mfem_mgis::Context& ctx,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& fs,
      const Layout& fl,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& vs,
      const Layout& vl,
      const mfem_mgis::size_type nc) {
    using namespace mfem_mgis;
    auto fvalues = std::vector<real>(getSpaceSize(*fs) * fl.stride, sentinel);
    auto vvalues = std::vector<real>(getSpaceSize(*vs) * vl.stride, sentinel);
    auto f = PartialQuadratureFunctionView(
        fs, fvalues,
        {.data_begin = fl.begin, .data_size = nc, .data_stride = fl.stride});
    auto v = PartialQuadratureFunctionView(
        vs, vvalues,
        {.data_begin = vl.begin, .data_size = nc, .data_stride = vl.stride});
    const auto vref = evaluateLinearField(vs);
    for (const auto& [i, o] : vs->getOffsets()) {
      static_cast<void>(o);
      for (size_type g = 0; g != vs->getNumberOfQuadraturePoints(i); ++g) {
        auto* const values = v.data(i, g);
        for (size_type c = 0; c != nc; ++c) {
          values[c] = vref->getIntegrationPointValue(i, g) + 10 * c;
        }
      }
    }
    if (!assign_values(ctx, f, v)) {
      return {};
    }
    // the stored values are read directly to check the layout
    const auto fref = evaluateLinearField(fs);
    auto e = real{};
    for (const auto& [i, o] : fs->getOffsets()) {
      for (size_type g = 0; g != fs->getNumberOfQuadraturePoints(i); ++g) {
        const auto* const values = fvalues.data() + (o + g) * fl.stride;
        for (size_type c = 0; c != fl.stride; ++c) {
          if ((c < fl.begin) || (c >= fl.begin + nc)) {
            e = std::max(e, std::abs(values[c] - sentinel));
          } else {
            const auto r =
                fref->getIntegrationPointValue(i, g) + 10 * (c - fl.begin);
            e = std::max(e, std::abs(values[c] - r));
          }
        }
      }
    }
    return e;
  }  // end of getAssignmentError
  /*!
   * \brief check the copies between views defined on two partial quadrature
   * spaces, in both directions, for one and two components, for contiguous
   * views and for views with a stride
   * \param[in, out] ctx: execution context
   * \param[in] s1: partial quadrature space
   * \param[in] s2: partial quadrature space
   * \param[in] exact: if false, the copies may fail, but shall be exact if
   * they succeed
   */
  void checkAssignments(
      mfem_mgis::Context& ctx,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& s1,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& s2,
      const bool exact) {
    using namespace mfem_mgis;
    for (const size_type nc : {1, 2}) {
      const auto contiguous = Layout{.begin = 0, .stride = nc};
      const auto strided = Layout{.begin = 1, .stride = 3};
      const auto strided2 = Layout{.begin = 0, .stride = 3};
      const auto layouts =
          std::array<std::array<Layout, 2>, 4>{{{contiguous, contiguous},
                                                {contiguous, strided},
                                                {strided, contiguous},
                                                {strided, strided2}}};
      for (const auto& [fl, vl] : layouts) {
        for (const auto& [fs, vs] : {std::array{s1, s2}, std::array{s2, s1}}) {
          const auto oe = getAssignmentError(ctx, fs, fl, vs, vl, nc);
          if (isInvalid(oe)) {
            TFEL_TESTS_CHECK(!exact);
            ctx.clearErrorMessages();
            continue;
          }
          TFEL_TESTS_CHECK(*oe < eps);
        }
      }
    }
  }  // end of checkAssignments
};

TFEL_TESTS_GENERATE_PROXY(AssignValuesTest, "AssignValuesTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parameters.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (parameters.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  const auto mesh = std::filesystem::path(parameters.mesh_file).stem().string();
  if (parameters.parallel) {
    m.addXMLTestOutput(
        "ParallelAssignValuesTest-" + mesh + "-" +
        std::to_string(mfem_mgis::getMPISizeForDefaultCommunicator()) + "-" +
        std::to_string(mfem_mgis::getMPIRankForDefaultCommunicator()) + ".xml");
  } else {
    m.addXMLTestOutput("AssignValuesTest-" + mesh + ".xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
