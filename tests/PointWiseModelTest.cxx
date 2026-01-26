/*!
 * \file   PointWiseModelTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   26/01/2026
 */

#include <cstdlib>
#include "MGIS/Model/Model.hxx"
#include "MFEMMGIS/Material.hxx"
#include "UnitTestingUtilities.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int linearsolver = 0;
  int order = 1;
  int parallel = 0;
};  // end of struct TestParameters

static void parseCommandLineOptions(TestParameters& params,
                                    int argc,
                                    char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&params.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&params.linearsolver, "-ls", "--linearsolver",
                 "identifier of the linear solver: 0 -> CG, 1 -> GMRES, 2 -> "
                 "UMFPack (serial), 3-> MUMPS(serial), 2 -> HypreFGMRES "
                 "(//), 3 -> HyprePCG (//), 4 -> HypreGMRES (//)");
  args.AddOption(&params.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&params.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (params.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

// single component test
template <bool parallel>
bool test(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto qspace = std::make_shared<PartialQuadratureSpace>(
      fed, 5,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation&) noexcept
      -> const mfem::IntegrationRule& {
        return mfem::IntRules.Get(e.GetGeomType(), 0);
      });
  const char* library = "src/libBehaviour.so";
  auto omodel = mgis::model::load(ctx, library, "UO2_Shrinkage_RAPHAEL2008",
                                  mgis::behaviour::Hypothesis::PLANESTRAIN);
  if (isInvalid(omodel)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  auto m = Material(qspace, std::make_unique<Behaviour>(*omodel));
  return true;
}  // end of test

int main(int argc, char**argv) {
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