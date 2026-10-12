/*!
 * \file   tests/GridFunctionOnSubMeshTest.cxx
 * \brief  This test checks that a grid function can only be built from partial
 * quadrature functions defined on all the materials or on all the boundaries
 * of its mesh, and that partial quadrature functions are correctly updated
 * from a linear field defined on another mesh.
 * \date   18/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <memory>
#include <optional>
#include <filesystem>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "mfem/fem/coefficient.hpp"
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
 * \return the integration rule used by the partial quadrature spaces
 * \param[in] e: finite element
 * \param[in] tr: element transformation
 */
static const mfem::IntegrationRule& getIntegrationRule(
    const mfem::FiniteElement& e,
    const mfem::ElementTransformation& tr) noexcept {
  return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e));
}  // end of getIntegrationRule

struct GridFunctionOnSubMeshTest final : public tfel::tests::TestCase {
  GridFunctionOnSubMeshTest()
      : tfel::tests::TestCase("MFEMMGIS", "GridFunctionOnSubMeshTest") {
  }  // end of GridFunctionOnSubMeshTest
  tfel::tests::TestResult execute() override {
    if (parameters.parallel) {
#ifdef MFEM_USE_MPI
      this->test1<true>();
      this->test2<true>();
      this->test3<true>();
#else  /* MFEM_USE_MPI */
      mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
      this->test2<false>();
      this->test3<false>();
    }
    return this->result;
  }  // end of execute

 private:
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
  //! \brief function defined on a material
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    auto f = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
    auto c = mfem::ConstantCoefficient(1);
    f.ProjectCoefficient(c);
    // a function defined on the first material only
    auto qspace = std::make_shared<PartialQuadratureSpace>(
        fed, MaterialIdentifier{.id = 1}, &getIntegrationRule);
    auto qf = PartialQuadratureFunction(qspace, 1);
    TFEL_TESTS_ASSERT(update(ctx, qf, f));
    // the whole mesh contains a material on which no function is defined
    auto ctx2 = Context{};
    TFEL_TESTS_CHECK(isInvalid(makeGridFunction<parallel>(ctx2, {qf})));
    TFEL_TESTS_CHECK(isInvalid(
        makeGridFunction<parallel>(ctx2, {qf}, fed.getMesh<parallel>())));
    // submesh of another material
    const auto submesh2 = fed.getSubMesh<parallel>(
        ctx, 2, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(submesh2));
    TFEL_TESTS_CHECK(
        isInvalid(makeGridFunction<parallel>(ctx2, {qf}, *submesh2)));
    // submesh of a boundary having the same attribute than the material
    const auto bsubmesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(bsubmesh));
    TFEL_TESTS_CHECK(
        isInvalid(makeGridFunction<parallel>(ctx2, {qf}, *bsubmesh)));
    // submesh of the material
    const auto submesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(submesh));
    auto og = makeGridFunction<parallel>(ctx, {qf}, *submesh);
    TFEL_TESTS_ASSERT(isValid(og));
    TFEL_TESTS_ASSERT(updateGridFunction<parallel>(ctx, *og, {qf}));
    // no function given
    TFEL_TESTS_CHECK(!updateGridFunction<parallel>(ctx2, *og, {}));
    *og -= 1;
    TFEL_TESTS_CHECK(og->Normlinf() < 1e-10);
  }  // end of test1
  //! \brief function defined on a boundary
  template <bool parallel>
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    // a function defined on the second boundary
    auto qspace = std::make_shared<PartialQuadratureSpace>(
        fed, BoundaryIdentifier{.id = 2}, &getIntegrationRule);
    auto qf = PartialQuadratureFunction(qspace, 1);
    // the quadrature space is built on the boundary elements of the whole mesh
    const auto omesh = qspace->getMesh<parallel>(ctx);
    TFEL_TESTS_ASSERT(isValid(omesh));
    TFEL_TESTS_CHECK(&(*omesh) == &(fed.getMesh<parallel>()));
    const auto bsubmesh = fed.getSubMesh<parallel>(
        ctx, 2, MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(bsubmesh));
    // a grid function on the whole mesh can be used to update the function
    auto c = mfem::ConstantCoefficient(1);
    auto ctx2 = Context{};
    auto f0 = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
    f0.ProjectCoefficient(c);
    TFEL_TESTS_CHECK(update(ctx2, qf, f0));
    for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
      TFEL_TESTS_CHECK(std::abs(*(qf.data(i)) - 1) < 1e-14);
    }
    // grid function on the submesh of the boundary
    const auto fespace =
        fed.getFiniteElementSpacesManager().getFiniteElementSpace<parallel>(
            ctx, *bsubmesh, 1);
    TFEL_TESTS_ASSERT(isValid(fespace));
    auto f = GridFunction<parallel>{fespace.get()};
    f.ProjectCoefficient(c);
    for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
      *(qf.data(i)) = 0;
    }
    TFEL_TESTS_ASSERT(update(ctx, qf, f));
    for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
      TFEL_TESTS_CHECK(std::abs(*(qf.data(i)) - 1) < 1e-14);
    }
    // the function can only be projected on the submesh of the boundary
    TFEL_TESTS_CHECK(isInvalid(makeGridFunction<parallel>(ctx2, {qf})));
    TFEL_TESTS_CHECK(isInvalid(
        makeGridFunction<parallel>(ctx2, {qf}, fed.getMesh<parallel>())));
    // submesh of a material having the same attribute than the boundary
    const auto submesh = fed.getSubMesh<parallel>(
        ctx, 2, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(submesh));
    TFEL_TESTS_CHECK(
        isInvalid(makeGridFunction<parallel>(ctx2, {qf}, *submesh)));
    // functions defined on materials and on boundaries can not be mixed
    auto qspace2 = std::make_shared<PartialQuadratureSpace>(
        fed, MaterialIdentifier{.id = 2}, &getIntegrationRule);
    auto qf2 = PartialQuadratureFunction(qspace2, 1);
    TFEL_TESTS_CHECK(isInvalid(
        makeGridFunction<parallel>(ctx2, {qf, qf2}, fed.getMesh<parallel>())));
    // submesh of the boundary
    auto og = makeGridFunction<parallel>(ctx, {qf}, *bsubmesh);
    TFEL_TESTS_ASSERT(isValid(og));
    TFEL_TESTS_ASSERT(updateGridFunction<parallel>(ctx, *og, {qf}));
    *og -= 1;
    TFEL_TESTS_CHECK(og->Normlinf() < 1e-10);
  }  // end of test2
  //! \brief update of functions from a linear field defined on another mesh
  template <bool parallel>
  void test3() {
    using namespace mfem_mgis;
    constexpr auto eps = real{1e-12};
    auto ctx = Context{};
    auto fed = makeFiniteElementDiscretization(ctx);
    auto& fespace = fed.getFiniteElementSpace<parallel>();
    // submeshes and finite element spaces
    const auto bsubmesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_BOUNDARIES);
    const auto bsubmesh2 =
        fed.getSubMesh<parallel>(ctx, Parameter{std::vector<Parameter>{1, 3}},
                                 MeshDiscretization::Location::ON_BOUNDARIES);
    const auto msubmesh = fed.getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(bsubmesh));
    TFEL_TESTS_ASSERT(isValid(bsubmesh2));
    TFEL_TESTS_ASSERT(isValid(msubmesh));
    auto manager = fed.getFiniteElementSpacesManager();
    const auto bfespace =
        manager.getFiniteElementSpace<parallel>(ctx, *bsubmesh, 1);
    const auto bfespace2 =
        manager.getFiniteElementSpace<parallel>(ctx, *bsubmesh2, 1);
    const auto mfespace =
        manager.getFiniteElementSpace<parallel>(ctx, *msubmesh, 1);
    TFEL_TESTS_ASSERT(isValid(bfespace));
    TFEL_TESTS_ASSERT(isValid(bfespace2));
    TFEL_TESTS_ASSERT(isValid(mfespace));
    // partial quadrature spaces
    const auto b = BoundaryIdentifier{.id = 1};
    const auto m = MaterialIdentifier{.id = 1};
    const auto bqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, b, &getExactIntegrationRule);
    const auto bqspace2 = std::make_shared<const PartialQuadratureSpace>(
        fed, *bfespace, b, &getExactIntegrationRule);
    const auto bqspace3 = std::make_shared<const PartialQuadratureSpace>(
        fed, *mfespace, b, &getExactIntegrationRule);
    const auto mqspace = std::make_shared<const PartialQuadratureSpace>(
        fed, m, &getExactIntegrationRule);
    const auto mqspace2 = std::make_shared<const PartialQuadratureSpace>(
        fed, *mfespace, m, &getExactIntegrationRule);
    // boundary space on the whole mesh, grid function on the submesh of the
    // boundary, and conversely
    checkUpdate<parallel>(ctx, bqspace, *bfespace, eps);
    checkUpdate<parallel>(ctx, bqspace2, fespace, eps);
    // material space on the whole mesh, grid function on the submesh of the
    // material, and conversely
    checkUpdate<parallel>(ctx, mqspace, *mfespace, eps);
    checkUpdate<parallel>(ctx, mqspace2, fespace, eps);
    // two different submeshes
    checkUpdate<parallel>(ctx, bqspace2, *bfespace2, eps);
    // boundary seen through the submesh of a material. In 3D, some boundary
    // elements of this submesh are not oriented as in the whole mesh, so the
    // update shall either fail or be exact. In 2D, the boundary elements of
    // an exterior boundary are oriented as in the whole mesh. The boundary 1
    // of the 2D test mesh is exterior and inside the material 1, so the
    // update shall be exact.
    const auto exact = fed.getMesh<parallel>().Dimension() == 2;
    const auto check = [exact](const std::optional<real>& oe) {
      if (isInvalid(oe)) {
        return !exact;
      }
      return *oe < eps;
    };
    TFEL_TESTS_CHECK(check(getUpdateError<parallel>(ctx, bqspace, *mfespace)));
    TFEL_TESTS_CHECK(check(getUpdateError<parallel>(ctx, bqspace3, fespace)));
    // extrapolation at nodes on the submesh of the boundary
    const auto f = evaluateLinearField(bqspace);
    auto og = makeGridFunction<parallel>(ctx, {*f}, *bsubmesh);
    TFEL_TESTS_ASSERT(isValid(og));
    TFEL_TESTS_ASSERT(updateGridFunction<parallel>(ctx, *og, {*f}));
    auto c = mfem::FunctionCoefficient(getLinearField);
    TFEL_TESTS_CHECK(og->ComputeMaxError(c) < eps);
  }  // end of test3
  /*!
   * \return the value of a linear field
   * \param[in] x: position
   */
  static mfem_mgis::real getLinearField(const mfem::Vector& x) noexcept {
    return 1 + x[0] + 2 * x[1] + (x.Size() > 2 ? 3 * x[2] : mfem_mgis::real{0});
  }  // end of getLinearField
  /*!
   * \return an integration rule which integrates exactly the local
   * projections of a linear field
   * \param[in] e: finite element
   * \param[in] tr: element transformation
   */
  static const mfem::IntegrationRule& getExactIntegrationRule(
      const mfem::FiniteElement& e,
      const mfem::ElementTransformation& tr) noexcept {
    return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e) + 2);
  }  // end of getExactIntegrationRule
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
   * \return the maximum error of the update of a partial quadrature function
   * from the interpolation of the linear field on the given finite element
   * space, or nothing if the update fails
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   * \param[in] fespace: finite element space
   */
  template <bool parallel>
  static std::optional<mfem_mgis::real> getUpdateError(
      mfem_mgis::Context& ctx,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& s,
      mfem_mgis::FiniteElementSpace<parallel>& fespace) {
    using namespace mfem_mgis;
    auto c = mfem::FunctionCoefficient(getLinearField);
    auto f = GridFunction<parallel>{&fespace};
    f.ProjectCoefficient(c);
    auto qf = PartialQuadratureFunction(s, 1);
    if (!update(ctx, qf, f)) {
      return {};
    }
    const auto ref = evaluateLinearField(s);
    auto e = real{};
    for (const auto& [i, o] : s->getOffsets()) {
      static_cast<void>(o);
      for (size_type g = 0; g != s->getNumberOfQuadraturePoints(i); ++g) {
        e = std::max(e, std::abs(qf.getIntegrationPointValue(i, g) -
                                 ref->getIntegrationPointValue(i, g)));
      }
    }
    return e;
  }  // end of getUpdateError
  /*!
   * \brief check that the update of a partial quadrature function from the
   * interpolation of the linear field on the given finite element space
   * succeeds and is exact
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   * \param[in] fespace: finite element space
   * \param[in] eps: tolerance
   */
  template <bool parallel>
  void checkUpdate(
      mfem_mgis::Context& ctx,
      const std::shared_ptr<const mfem_mgis::PartialQuadratureSpace>& s,
      mfem_mgis::FiniteElementSpace<parallel>& fespace,
      const mfem_mgis::real eps) {
    using namespace mfem_mgis;
    const auto oe = getUpdateError<parallel>(ctx, s, fespace);
    TFEL_TESTS_ASSERT(isValid(oe));
    TFEL_TESTS_CHECK(*oe < eps);
  }  // end of checkUpdate
};

TFEL_TESTS_GENERATE_PROXY(GridFunctionOnSubMeshTest,
                          "GridFunctionOnSubMeshTest");

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
        "ParallelGridFunctionOnSubMeshTest-" + mesh + "-" +
        std::to_string(mfem_mgis::getMPISizeForDefaultCommunicator()) + "-" +
        std::to_string(mfem_mgis::getMPIRankForDefaultCommunicator()) + ".xml");
  } else {
    m.addXMLTestOutput("GridFunctionOnSubMeshTest-" + mesh + ".xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
