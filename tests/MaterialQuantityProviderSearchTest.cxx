/*!
 * \file   tests/MaterialQuantityProviderSearchTest.cxx
 * \brief  This test checks the search of the behaviour integrator providing a
 * gradient, a thermodynamic force or an internal state variable when several
 * behaviour integrators are defined on the same material, and the export of
 * the integration points results at nodes in this case.
 * \date   23/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <string>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/MaterialQuantityProviderSearch.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int parallel = 0;
  int order = 1;
} parameters;

struct MaterialQuantityProviderSearchTest final : public tfel::tests::TestCase {
  MaterialQuantityProviderSearchTest()
      : tfel::tests::TestCase("MFEMMGIS",
                              "MaterialQuantityProviderSearchTest") {
  }  // end of MaterialQuantityProviderSearchTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }

 private:
  using Status = mfem_mgis::MaterialQuantityProviderSearchResult::Status;
  //! \return the parameters of the non linear evolution problem
  static mfem_mgis::Parameters getProblemParameters() {
    return {{"MeshFileName", parameters.mesh_file},
            {"FiniteElementFamily", "H1"},
            {"FiniteElementOrder", parameters.order},
            {"UnknownsSize", 3},
            {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
            {"Hypothesis", "Tridimensional"},
            {"Parallel", bool(parameters.parallel)}};
  }
  /*!
   * \brief add two behaviour integrators on the material `1`:
   *
   * - the first one, based on the `Elasticity` behaviour, has no internal
   *   state variable.
   * - the second one, based on the `Plasticity` behaviour, has two internal
   *   state variables: `ElasticStrain` and `EquivalentPlasticStrain`.
   *
   * Both behaviours define the `Strain` gradient and the `Stress`
   * thermodynamic force.
   *
   * \param[in] p: non linear evolution problem
   */
  void addBehaviourIntegrators(mfem_mgis::NonLinearEvolutionProblem& p) {
    auto ctx = mfem_mgis::Context{};
    for (const auto& b : {"Elasticity", "Plasticity"}) {
      TFEL_TESTS_ASSERT(mfem_mgis::isValid(p.addBehaviourIntegrator(
          ctx, "Mechanics", 1, parameters.library, b)));
    }
  }
  //! \brief check the search functions
  void test1() {
    using namespace mfem_mgis;
    auto ctx0 = Context{};
    auto oproblem =
        construct<NonLinearEvolutionProblem>(ctx0, getProblemParameters());
    TFEL_TESTS_ASSERT(isValid(oproblem));
    this->addBehaviourIntegrators(*oproblem);
    const auto& p = *oproblem;
    auto ctx = Context{};
    const auto material =
        LocationIdentifier{.material_identifier = MaterialIdentifier{.id = 1},
                           .boundary_identifier = {}};
    // a unique provider
    const auto r1 = hasInternalStateVariableProvider(ctx, p, material,
                                                     "EquivalentPlasticStrain");
    TFEL_TESTS_ASSERT(isValid(r1));
    TFEL_TESTS_CHECK(r1->status == Status::SUCCESS);
    TFEL_TESTS_ASSERT(isValid(*r1));
    const auto obi = p.getBehaviourIntegrator(ctx, 1, 1);
    TFEL_TESTS_ASSERT(isValid(obi));
    TFEL_TESTS_CHECK(&(*(r1->behaviour_integrator)) == &(*obi));
    // multiple providers
    const auto r2 = hasGradientProvider(ctx, p, material, "Strain");
    TFEL_TESTS_ASSERT(isValid(r2));
    TFEL_TESTS_CHECK(r2->status == Status::MULTIPLE_PROVIDERS);
    TFEL_TESTS_CHECK(isInvalid(*r2));
    const auto r3 = hasThermodynamicForceProvider(ctx, p, material, "Stress");
    TFEL_TESTS_ASSERT(isValid(r3));
    TFEL_TESTS_CHECK(r3->status == Status::MULTIPLE_PROVIDERS);
    TFEL_TESTS_CHECK(isInvalid(*r3));
    // no provider
    const auto r4 =
        hasGradientProvider(ctx, p, material, "EquivalentPlasticStrain");
    TFEL_TESTS_ASSERT(isValid(r4));
    TFEL_TESTS_CHECK(r4->status == Status::NO_PROVIDER);
    TFEL_TESTS_CHECK(isInvalid(*r4));
    // behaviour integrators are not defined on boundaries
    const auto boundary =
        LocationIdentifier{.material_identifier = {},
                           .boundary_identifier = BoundaryIdentifier{.id = 1}};
    const auto r5 = hasInternalStateVariableProvider(ctx, p, boundary,
                                                     "EquivalentPlasticStrain");
    TFEL_TESTS_ASSERT(isValid(r5));
    TFEL_TESTS_CHECK(r5->status == Status::NO_PROVIDER);
    // invalid location
    auto ctx2 = Context{};
    const auto r6 =
        hasGradientProvider(ctx2, p, LocationIdentifier{}, "Strain");
    TFEL_TESTS_CHECK(isInvalid(r6));
    TFEL_TESTS_CHECK(ctx2.getRawErrorMessage() ==
                     "invalid location identifier");
  }
  /*!
   * \brief check the export of integration points results at nodes when
   * several behaviour integrators are defined on the same material
   */
  void test2() {
    using namespace mfem_mgis;
    auto ctx0 = Context{};
    auto oproblem =
        construct<NonLinearEvolutionProblem>(ctx0, getProblemParameters());
    TFEL_TESTS_ASSERT(isValid(oproblem));
    this->addBehaviourIntegrators(*oproblem);
    auto& p = *oproblem;
    const auto prefix = std::string{"MaterialQuantityProviderSearchTest-"} +
                        (parameters.parallel ? "parallel" : "sequential");
    auto add_post_processing = [&p, &prefix](Context& ctx, std::string_view r) {
      return p.addPostProcessing(
          ctx, "ParaviewExportIntegrationPointResultsAtNodes",
          {{"OutputFileName", prefix + "-" + std::string{r}},
           {"Results", std::string{r}}});
    };
    auto contains = [](const std::string& s, std::string_view s2) {
      return s.find(s2) != std::string::npos;
    };
    // the internal state variable is only provided by the second behaviour
    // integrator
    auto ctx = Context{};
    TFEL_TESTS_CHECK(add_post_processing(ctx, "EquivalentPlasticStrain"));
    TFEL_TESTS_CHECK(p.executeInitialPostProcessings(ctx, 0));
    // the stress is provided by both behaviour integrators
    auto ctx2 = Context{};
    TFEL_TESTS_CHECK(!add_post_processing(ctx2, "Stress"));
    TFEL_TESTS_CHECK(
        contains(ctx2.getErrorMessage(),
                 "multiple behaviour integrators provide the result 'Stress' "
                 "on material '1'"));
    // no behaviour integrator provides this result
    auto ctx3 = Context{};
    TFEL_TESTS_CHECK(!add_post_processing(ctx3, "UnknownResult"));
    TFEL_TESTS_CHECK(contains(ctx3.getErrorMessage(),
                              "no behaviour integrator provides the result "
                              "'UnknownResult' on material '1'"));
  }
};

TFEL_TESTS_GENERATE_PROXY(MaterialQuantityProviderSearchTest,
                          "MaterialQuantityProviderSearchTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.library, "-l", "--library", "Material library.");
  args.AddOption(&parameters.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parameters.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (parameters.mesh_file == nullptr) ||
      (parameters.library == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (parameters.parallel) {
    m.addXMLTestOutput("ParallelMaterialQuantityProviderSearchTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("MaterialQuantityProviderSearchTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
