/*!
 * \file   PointWiseModelTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   26/01/2026
 */

#include <cstdlib>
#include "MGIS/Model/Model.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
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
  const char* library = "./libBehaviourTest.so";
  auto omodel = mgis::model::load(ctx, library, "UO2Shrinkage_RAPHAEL2008",
                                  mgis::behaviour::Hypothesis::PLANESTRAIN);
  if (isInvalid(omodel)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  auto m = Material(qspace, std::make_unique<Behaviour>(*omodel));
  mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 893.15);
  mgis::behaviour::setExternalStateVariable(m.s0, "BurnUp_AtPercent", 0);
  mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 893.15);
  mgis::behaviour::setExternalStateVariable(m.s1, "BurnUp_AtPercent", 5);
  const auto r = mgis::behaviour::integrate(
      m, mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, 1,
      0, m.n);
  if (!((r == 1) || (r == 0))) {
    return false;
  }
  const auto s = getInternalStateVariable(m, "Shrinkage", ets);
  const auto Ta = real{750};
  const auto Tm = real{893.15};
  const auto A = std::max(real{5.e-3}, real{-1.26e-2 + 1.8e-5 * Tm});
  const auto expTm = real{1} + exp(-Ta / Tm);
  const auto BUp = real{5};
  const auto s_ref = A * (exp(-3.e-2 * BUp * expTm) - 1);
  for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
    const auto value = s.getIntegrationPointValue(i);
    if (std::abs(value - s_ref) > 1e-8) {
      return ctx.registerErrorMessage("invalid shrinkage value (" +
                                      std::to_string(value) + ", expected " +
                                      std::to_string(s_ref) + ")");
    }
  }
  return true;
}  // end of test

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
