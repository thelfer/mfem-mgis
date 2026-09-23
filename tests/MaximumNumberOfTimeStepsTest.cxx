/*!
 * \file   tests/MaximumNumberOfTimeStepsTest.cxx
 * \author Thomas Helfer
 * \date   06/04/2023
 */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LoopCouplingScheme.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/Simulation.hxx"

const char* mesh_file = nullptr;

struct MaximumNumberOfTimeStepsTest final : public tfel::tests::TestCase {
  MaximumNumberOfTimeStepsTest()
      : tfel::tests::TestCase("MFEMMGIS", "MaximumNumberOfTimeStepsTest") {
  }  // end of MaximumNumberOfTimeStepsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mfem_mgis;
    const auto mesh_parameters =
        Parameters{{"generator", "Grid"}, {"spaceDimension", 3}};
    auto ctx = Context{};
    auto mesh = construct<MeshDiscretization>(
        ctx, Parameters{{"MeshFileName", mesh_file}, {"Parallel", false}});
    TFEL_TESTS_ASSERT(isValid(mesh));
    auto ps = construct<PhysicalSystem>(ctx, *mesh);
    TFEL_TESTS_ASSERT(isValid(ps));
    // creating the coupling scheme
    auto c = make_shared<LoopCouplingScheme>(
        ctx, *mesh, Parameters{{"NumberOfIterations", 1}});
    TFEL_TESTS_ASSERT(isValid(c));
    TFEL_TESTS_ASSERT(ps->setCouplingScheme(ctx, c));
    constexpr auto nSteps = size_type{10};
    auto s = construct<Simulation>(ctx, *ps,
                                   Simulation::TimesDescription{0, 1, nSteps});
    TFEL_TESTS_ASSERT(isValid(s));
    TFEL_TESTS_CHECK(isInvalid(s->setMaximumNumberOfTimeSteps(ctx, 0)));
    ctx.clearErrorMessages();
    TFEL_TESTS_CHECK(isInvalid(s->setMaximumNumberOfTimeSteps(ctx, -1)));
    ctx.clearErrorMessages();
    TFEL_TESTS_CHECK(s->setMaximumNumberOfTimeSteps(ctx, 2));
    TFEL_TESTS_CHECK(s->run(ctx).first.shallContinue());
    TFEL_TESTS_ASSERT(s->getTimes().size() == 9);
    TFEL_TESTS_CHECK(std::abs(*(s->getTimes().cbegin()) - 0.2) < 1e-14);
    TFEL_TESTS_CHECK(s->setMaximumNumberOfTimeSteps(ctx, 4));
    TFEL_TESTS_CHECK(s->run(ctx).first.shallContinue());
    TFEL_TESTS_ASSERT(s->getTimes().size() == 5);
    TFEL_TESTS_CHECK(std::abs(*(s->getTimes().cbegin()) - 0.6) < 1e-14);
    s->unsetMaximumNumberOfTimeSteps();
    TFEL_TESTS_CHECK(s->run(ctx).first.shallContinue());
    TFEL_TESTS_ASSERT(s->getTimes().size() == 1);
    TFEL_TESTS_CHECK(s->run(ctx).first.shallStop());
  }  // en of test1
};

TFEL_TESTS_GENERATE_PROXY(MaximumNumberOfTimeStepsTest,
                          "MaximumNumberOfTimeStepsTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  //
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  if (mesh_file == nullptr) {
    mfem_mgis::getOutputStream() << "no mesh file specified\n";
    return EXIT_FAILURE;
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("MaximumNumberOfTimeStepsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
