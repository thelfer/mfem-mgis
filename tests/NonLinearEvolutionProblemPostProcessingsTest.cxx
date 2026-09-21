/*!
 * \file   tests/NonLinearEvolutionProblemPostProcessingsTest.cxx
 * \brief  This test checks that the post-processings of a non linear evolution
 * problem are executed at the initial time and at the end of the time steps.
 * \date   21/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct {
  const char* mesh_file = nullptr;
  int parallel = 0;
  int order = 1;
} parameters;

/*!
 * \brief a post-processing which records the times at which it is executed
 * \tparam parallel: boolean stating if the computations are performed in
 * parallel.
 */
template <bool parallel>
struct TestPostProcessing final
    : mfem_mgis::AbstractNonLinearEvolutionProblemPostProcessing<parallel> {
  /*!
   * \brief constructor
   * \param[in] ts: recorded times
   * \param[in] b: boolean stating if the initial post-processing succeeds
   */
  TestPostProcessing(std::vector<mfem_mgis::real>& ts, const bool b)
      : times(ts), success(b) {}  // end of TestPostProcessing
  //
  bool executeInitialPostProcessing(
      mfem_mgis::Context& ctx,
      mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>&,
      const mfem_mgis::real t) noexcept override {
    if (!this->success) {
      return ctx.registerErrorMessage("initial post-processing failed");
    }
    this->times.push_back(t);
    return true;
  }  // end of executeInitialPostProcessing
  void execute(mfem_mgis::Context&,
               mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>&,
               const mfem_mgis::real t,
               const mfem_mgis::real dt) override {
    this->times.push_back(t + dt);
  }  // end of execute
  //! \brief destructor
  ~TestPostProcessing() override = default;

 private:
  //! \brief recorded times
  std::vector<mfem_mgis::real>& times;
  //! \brief boolean stating if the initial post-processing succeeds
  const bool success;
};  // end of struct TestPostProcessing

struct NonLinearEvolutionProblemPostProcessingsTest final
    : public tfel::tests::TestCase {
  NonLinearEvolutionProblemPostProcessingsTest()
      : tfel::tests::TestCase("MFEMMGIS",
                              "NonLinearEvolutionProblemPostProcessingsTest") {
  }  // end of NonLinearEvolutionProblemPostProcessingsTest
  tfel::tests::TestResult execute() override {
    if (parameters.parallel) {
#ifdef MFEM_USE_MPI
      this->test1<true>();
#else  /* MFEM_USE_MPI */
      mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
    }
    return this->result;
  }

 private:
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto oproblem = construct<NonLinearEvolutionProblem>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", parameters.order},
                  {"UnknownsSize", 3},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Hypothesis", "Tridimensional"},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(oproblem));
    auto& p = oproblem->getImplementation<parallel>();
    // a post-processing which succeeds
    auto times = std::vector<real>{};
    TFEL_TESTS_ASSERT(p.addPostProcessing(
        ctx, std::make_unique<TestPostProcessing<parallel>>(times, true)));
    TFEL_TESTS_CHECK(oproblem->executeInitialPostProcessings(ctx, 1));
    oproblem->executePostProcessings(ctx, 1, 2);
    TFEL_TESTS_CHECK((times == std::vector<real>{1, 3}));
    // a post-processing whose initial post-processing fails
    auto times2 = std::vector<real>{};
    TFEL_TESTS_ASSERT(p.addPostProcessing(
        ctx, std::make_unique<TestPostProcessing<parallel>>(times2, false)));
    TFEL_TESTS_CHECK(!oproblem->executeInitialPostProcessings(ctx, 3));
    TFEL_TESTS_CHECK(ctx.getRawErrorMessage() ==
                     "initial post-processing failed");
    TFEL_TESTS_CHECK((times == std::vector<real>{1, 3, 3}));
    TFEL_TESTS_CHECK(times2.empty());
  }
};

TFEL_TESTS_GENERATE_PROXY(NonLinearEvolutionProblemPostProcessingsTest,
                          "NonLinearEvolutionProblemPostProcessingsTest");

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
    m.addXMLTestOutput("ParallelNonLinearEvolutionProblemPostProcessingsTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("NonLinearEvolutionProblemPostProcessingsTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
