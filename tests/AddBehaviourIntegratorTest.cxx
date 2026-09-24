/*!
 * \file   tests/AddBehaviourIntegratorTest.cxx
 * \brief  This test checks that the non throwing version of the
 * `addBehaviourIntegrator` method reports errors in the execution context.
 * \date   24/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct {
  const char* mesh_file = nullptr;
  int parallel = 0;
} parameters;

struct AddBehaviourIntegratorTest final : public tfel::tests::TestCase {
  AddBehaviourIntegratorTest()
      : tfel::tests::TestCase("MFEMMGIS", "AddBehaviourIntegratorTest") {
  }  // end of AddBehaviourIntegratorTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }

 private:
  //! \return the parameters of the non linear evolution problem
  static mfem_mgis::Parameters getProblemParameters() {
    return {{"MeshFileName", parameters.mesh_file},
            {"FiniteElementFamily", "H1"},
            {"FiniteElementOrder", 1},
            {"UnknownsSize", 3},
            {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
            {"Hypothesis", "Tridimensional"},
            {"Parallel", bool(parameters.parallel)}};
  }
  //! \brief unknown material
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto oproblem =
        construct<NonLinearEvolutionProblem>(ctx, getProblemParameters());
    TFEL_TESTS_ASSERT(isValid(oproblem));
    const auto r = oproblem->addBehaviourIntegrator(
        ctx, "Mechanics", "UnknownMaterial", "libBehaviour.so", "Elasticity");
    TFEL_TESTS_CHECK(isInvalid(r));
    TFEL_TESTS_CHECK(!ctx.getErrorMessage().empty());
  }
  //! \brief multi-material integrator disabled
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto params = getProblemParameters();
    params.insert(throwing, "UseMultiMaterialNonLinearIntegrator", false);
    auto oproblem = construct<NonLinearEvolutionProblem>(ctx, params);
    TFEL_TESTS_ASSERT(isValid(oproblem));
    const auto r = oproblem->addBehaviourIntegrator(
        ctx, "Mechanics", 1, "libBehaviour.so", "Elasticity");
    TFEL_TESTS_CHECK(isInvalid(r));
    TFEL_TESTS_CHECK(ctx.getRawErrorMessage() ==
                     "NonLinearEvolutionProblemImplementationBase::"
                     "addBehaviourIntegrator: multi material support has been "
                     "disabled");
  }
};

TFEL_TESTS_GENERATE_PROXY(AddBehaviourIntegratorTest,
                          "AddBehaviourIntegratorTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
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
    m.addXMLTestOutput("ParallelAddBehaviourIntegratorTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("AddBehaviourIntegratorTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
