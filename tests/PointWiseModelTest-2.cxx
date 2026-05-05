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
#include "MFEMMGIS/PointWiseModel.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/TimeStep.hxx"
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
  auto or_die = ctx.getFatalFailureHandler();
  auto fed =
      make_shared<FiniteElementDiscretization>(
          ctx, Parameters{{"MeshFileName", params.mesh_file},
                          {"FiniteElementFamily", "H1"},
                          {"FiniteElementOrder", params.order},
                          {"UnknownsSize", 1},
                          {"NumberOfUniformRefinements", parallel ? 1 : 0},
                          {"Parallel", parallel}}) |
      or_die;
  auto qspace = make_shared<const PartialQuadratureSpace>(
                    ctx, *fed, 5,
                    [](const mfem::FiniteElement& e,
                       const mfem::ElementTransformation&) noexcept
                    -> const mfem::IntegrationRule& {
                      return mfem::IntRules.Get(e.GetGeomType(), 0);
                    }) |
                or_die;
  auto model = make_shared<PointWiseModel>(
                   ctx, qspace,
                   Parameters{{"Library", "./libBehaviourTest.so"},
                              {"Model", "UO2Shrinkage_RAPHAEL2008"},
                              {"Hypothesis", "PlaneStrain"}}) |
               or_die;
  auto& m = model->getMaterial();
  mgis::behaviour::setExternalStateVariable(ctx, m.s0, "Temperature", 893.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m.s0, "BurnUp_AtPercent", 0) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m.s1, "Temperature", 893.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m.s1, "BurnUp_AtPercent", 5) |
      or_die;
  //
  auto ps = construct<PhysicalSystem>(ctx, *fed) | or_die;
  ps.setModel(ctx, model) | or_die;
  const auto r = ps.computeNextState(ctx, {.begin = 0, .end = 1, .dt = 1});
  if (!r.first.shallContinue()) {
    std::cerr << ctx.getErrorMessage() << '\n';
    std::exit(-1);
  }
  //
  const auto s = getInternalStateVariable(ctx, m, "Shrinkage", ets) | or_die;
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
