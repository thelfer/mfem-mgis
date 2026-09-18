/*!
 * \file   GridFunctionOnSubMeshTest.cxx
 * \brief  This test checks that a grid function can only be built from partial
 * quadrature functions defined on all the materials of its mesh.
 * \date   18/09/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "UnitTestingUtilities.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int parallel = 0;
};  // end of struct TestParameters

static void parseCommandLineOptions(TestParameters& params,
                                    int argc,
                                    char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&params.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&params.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (params.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

template <bool parallel>
bool test(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      ctx,
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", 1},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto f = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
  auto c = mfem::ConstantCoefficient(1);
  f.ProjectCoefficient(c);
  // a function defined on the first material only
  auto qspace = std::make_shared<PartialQuadratureSpace>(
      fed, 1,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto qf = PartialQuadratureFunction(qspace, 1);
  if (!update(ctx, qf, f)) {
    return false;
  }
  // the whole mesh contains a material on which no function is defined
  auto ctx2 = Context{};
  if (!isInvalid(makeGridFunction<parallel>(ctx2, {qf}))) {
    return ctx.registerErrorMessage(
        "building a grid function on the whole mesh shall fail");
  }
  if (!isInvalid(
          makeGridFunction<parallel>(ctx2, {qf}, fed.getMesh<parallel>()))) {
    return ctx.registerErrorMessage(
        "building a grid function on the whole mesh shall fail");
  }
  // submesh of another material
  const auto submesh2 = fed.getSubMesh<parallel>(
      ctx, 2, MeshDiscretization::Location::ON_MATERIALS);
  if (isInvalid(submesh2)) {
    return false;
  }
  if (!isInvalid(makeGridFunction<parallel>(ctx2, {qf}, *submesh2))) {
    return ctx.registerErrorMessage(
        "building a grid function on the submesh of another material shall "
        "fail");
  }
  // submesh of a boundary having the same attribute than the material
  const auto bsubmesh = fed.getSubMesh<parallel>(
      ctx, 1, MeshDiscretization::Location::ON_BOUNDARIES);
  if (isInvalid(bsubmesh)) {
    return false;
  }
  if (!isInvalid(makeGridFunction<parallel>(ctx2, {qf}, *bsubmesh))) {
    return ctx.registerErrorMessage(
        "building a grid function on the submesh of a boundary shall fail");
  }
  // submesh of the material
  const auto submesh = fed.getSubMesh<parallel>(
      ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
  if (isInvalid(submesh)) {
    return false;
  }
  auto og = makeGridFunction<parallel>(ctx, {qf}, *submesh);
  if (isInvalid(og)) {
    return false;
  }
  updateGridFunction<parallel>(*og, {qf}, *submesh);
  *og -= 1;
  const auto e = og->Normlinf();
  if (e > 1e-10) {
    return ctx.registerErrorMessage("invalid values at nodes (maximum error " +
                                    std::to_string(e) + ")");
  }
  return true;
}

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // options treatment
  auto params = TestParameters{};
  initialize(argc, argv);
  parseCommandLineOptions(params, argc, argv);
  const auto success = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test<false>(ctx, params);
  }();
  if (!success) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
