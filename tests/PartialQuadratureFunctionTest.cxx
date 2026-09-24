/*!
 * \file   tests/PartialQuadratureFunctionTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <cstdlib>
#include <iostream>
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

auto parameters = mfem_mgis::unit_tests::TestParameters{};
auto ctx = mgis::Context{};

struct PartialQuadratureFunctionTest final : public tfel::tests::TestCase {
  PartialQuadratureFunctionTest()
      : tfel::tests::TestCase("MFEMMGIS", "PartialQuadratureFunctionTest") {
  }  // end of PartialQuadratureFunctionTest
  tfel::tests::TestResult execute() override {
    this->setup();
    this->test1<true>();
    this->test1<false>();
    this->test2();
    return this->result;
  }

 private:
  void setup() {
    using namespace mfem_mgis;
    auto or_raise = ctx.getThrowingFailureHandler();
    this->problem = make_unique<NonLinearEvolutionProblem>(
        ctx, ctx,
        Parameters{{"MeshFileName", parameters.mesh_file},
                   {"FiniteElementFamily", "H1"},
                   {"FiniteElementOrder", parameters.order},
                   {"UnknownsSize", 3},
                   {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                   {"Hypothesis", "Tridimensional"},
                   {"Parallel", bool(parameters.parallel)}});
    if (isInvalid(this->problem)) {
      return;
    }
    problem->addBehaviourIntegrator(ctx, "Mechanics", 1, parameters.library,
                                    parameters.behaviour) |
        or_raise;
    auto& m1 = problem->getMaterial(ctx, 1, 0) | or_raise;
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  }  // end of setup

  template <bool mutable_version>
  void test1() {
    using namespace mfem_mgis;
    using MaterialReference =
        std::conditional_t<mutable_version, Material&, const Material&>;
    using ExpectedResultType =
        std::conditional_t<mutable_version, PartialQuadratureFunction,
                           ImmutablePartialQuadratureFunctionView>;
    if (isInvalid(problem)) {
      return;
    }
    auto or_raise = ctx.getThrowingFailureHandler();
    auto& b = this->problem->getBehaviourIntegrator(ctx, 1, 0) | or_raise;
    const auto& qspace = b.getPartialQuadratureSpace();
    auto& mutable_material = b.getMaterial(ctx) | or_raise;
    auto& m1 = static_cast<MaterialReference&>(mutable_material);
    auto strain = getGradient(ctx, m1, "Strain") | or_raise;
    auto stress = getThermodynamicForce(ctx, m1, "Stress") | or_raise;
    auto p = getInternalStateVariable(ctx, m1, "EquivalentStrain") | or_raise;
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<decltype(strain), ExpectedResultType>));
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<decltype(stress), ExpectedResultType>));
    TFEL_TESTS_STATIC_ASSERT((std::same_as<decltype(p), ExpectedResultType>));
    TFEL_TESTS_ASSERT(strain.getNumberOfComponents() == 6);
    TFEL_TESTS_ASSERT(stress.getNumberOfComponents() == 6);
    TFEL_TESTS_ASSERT(p.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(strain.getDataStride() == 6);
    TFEL_TESTS_ASSERT(stress.getDataStride() == 6);
    TFEL_TESTS_ASSERT(p.getDataStride() == 1);
    if (getSpaceSize(qspace) != 0) {
      if constexpr (!mutable_version) {
        TFEL_TESTS_ASSERT(strain.getValues().data() == m1.s1.gradients.data());
        TFEL_TESTS_ASSERT(stress.getValues().data() ==
                          m1.s1.thermodynamic_forces.data());
        TFEL_TESTS_ASSERT(p.getValues().data() ==
                          m1.s1.internal_state_variables.data());
      }
    }
    // requesting an unknown internal state variable
    auto oeel = getInternalStateVariable(ctx, m1, "ElasticStrain");
    TFEL_TESTS_ASSERT(isInvalid(oeel));
  }  // end of test1
  //! \brief views on external data described by `ViewSpecifications`
  void test2() {
    using namespace mfem_mgis;
    if (isInvalid(problem)) {
      return;
    }
    auto or_raise = ctx.getThrowingFailureHandler();
    auto& b = this->problem->getBehaviourIntegrator(ctx, 1, 0) | or_raise;
    const auto om = std::as_const(b).getMaterial(ctx);
    TFEL_TESTS_ASSERT(isValid(om));
    const auto qspace = om->getPartialQuadratureSpacePointer();
    const auto n = getSpaceSize(*qspace);
    auto values = std::vector<real>(3 * n, real{0});
    for (size_type i = 0; i != n; ++i) {
      values[3 * i + 1] = static_cast<real>(i);
    }
    // the second component of an array of stride 3
    auto of = PartialQuadratureFunction::borrow(
        ctx, qspace, values,
        {.data_begin = 1, .data_size = 1, .data_stride = 3});
    TFEL_TESTS_ASSERT(isValid(of));
    TFEL_TESTS_CHECK(of->getDataStride() == 3);
    TFEL_TESTS_CHECK(of->getNumberOfComponents() == 1);
    for (size_type i = 0; i != n; ++i) {
      TFEL_TESTS_CHECK(std::abs(of->getIntegrationPointValue(i) -
                                static_cast<real>(i)) < 1e-14);
    }
    // the data range is outside the stride
    auto ctx2 = Context{};
    TFEL_TESTS_CHECK(isInvalid(PartialQuadratureFunction::borrow(
        ctx2, qspace, values,
        {.data_begin = 2, .data_size = 2, .data_stride = 3})));
    // the size of the values does not match the stride
    if (n != 0) {
      auto ctx3 = Context{};
      TFEL_TESTS_CHECK(isInvalid(PartialQuadratureFunction::borrow(
          ctx3, qspace, values,
          {.data_begin = 0, .data_size = 1, .data_stride = 2})));
    }
  }  // end of test2
  //
  std::unique_ptr<mfem_mgis::NonLinearEvolutionProblem> problem;
};

TFEL_TESTS_GENERATE_PROXY(PartialQuadratureFunctionTest,
                          "PartialQuadratureFunctionTest");

int main(int argc, char** argv) {
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);  // mfem_mgis::getDefaultLogStream());
  m.addXMLTestOutput("PartialQuadratureFunctionTest-" +
                     std::string(parameters.behaviour) + "-" +
                     std::to_string(mfem_mgis::getMPIsize()) + "-" +
                     std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
