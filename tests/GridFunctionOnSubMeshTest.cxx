/*!
 * \file   GridFunctionOnSubMeshTest.cxx
 * \brief  This test checks that a grid function can only be built from partial
 * quadrature functions defined on all the materials or on all the boundaries
 * of its mesh.
 * \date   18/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <cstdlib>
#include <iostream>
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

//! \return the integration rule used by the partial quadrature spaces
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
#else  /* MFEM_USE_MPI */
      mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
      this->test2<false>();
    }
    return this->result;
  }  // end of execute

 private:
  //! \return a finite element discretization built on the given mesh file
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
    auto qspace =
        std::make_shared<PartialQuadratureSpace>(fed, 1, &getIntegrationRule);
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
    updateGridFunction<parallel>(*og, {qf}, *submesh);
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
        fed,
        LocationIdentifier{.material_identifier = {},
                           .boundary_identifier = BoundaryIdentifier{.id = 2}},
        &getIntegrationRule);
    auto qf = PartialQuadratureFunction(qspace, 1);
    // the quadrature space is built on the submesh of the boundary
    const auto bsubmesh = fed.getSubMesh<parallel>(
        ctx, 2, MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(bsubmesh));
    const auto omesh = qspace->getMesh<parallel>(ctx);
    TFEL_TESTS_ASSERT(isValid(omesh));
    TFEL_TESTS_CHECK(&(*omesh) == &(*bsubmesh));
    // a grid function on the whole mesh can not be used to update the function
    auto c = mfem::ConstantCoefficient(1);
    auto ctx2 = Context{};
    auto f0 = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
    f0.ProjectCoefficient(c);
    TFEL_TESTS_CHECK(!update(ctx2, qf, f0));
    // grid function on the submesh of the boundary
    const auto fespace =
        fed.getFiniteElementSpacesManager().getFiniteElementSpace<parallel>(
            ctx, *bsubmesh, 1);
    TFEL_TESTS_ASSERT(isValid(fespace));
    auto f = GridFunction<parallel>{fespace.get()};
    f.ProjectCoefficient(c);
    TFEL_TESTS_ASSERT(update(ctx, qf, f));
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
    auto qspace2 =
        std::make_shared<PartialQuadratureSpace>(fed, 2, &getIntegrationRule);
    auto qf2 = PartialQuadratureFunction(qspace2, 1);
    TFEL_TESTS_CHECK(
        isInvalid(makeGridFunction<parallel>(ctx2, {qf, qf2}, *bsubmesh)));
    // submesh of the boundary
    auto og = makeGridFunction<parallel>(ctx, {qf}, *bsubmesh);
    TFEL_TESTS_ASSERT(isValid(og));
    updateGridFunction<parallel>(*og, {qf}, *bsubmesh);
    *og -= 1;
    TFEL_TESTS_CHECK(og->Normlinf() < 1e-10);
  }  // end of test2
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
  if (parameters.parallel) {
    m.addXMLTestOutput("ParallelGridFunctionOnSubMeshTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("GridFunctionOnSubMeshTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
