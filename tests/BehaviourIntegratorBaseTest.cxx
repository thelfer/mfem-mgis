/*!
 * \file   tests/BehaviourIntegratorBaseTest.cxx
 * \brief  This test checks that the behaviour integrators reject behaviours
 * whose type, kinematic, symmetry or hypothesis are not the expected ones.
 * \date   24/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>
#include <exception>
#include <string_view>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/TransientHeatTransferBehaviourIntegrator.hxx"

struct {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int parallel = 0;
} parameters;

/*!
 * \return the error message reported when building the given behaviour
 * integrator, or an empty string if the construction succeeded
 * \param[in] fed: finite element discretization
 * \param[in] b: behaviour
 */
template <typename BehaviourIntegrator>
static std::string getConstructionError(
    const mfem_mgis::FiniteElementDiscretization& fed,
    const mgis::behaviour::Behaviour& b) {
  try {
    auto bi = BehaviourIntegrator(
        fed, 1, std::make_unique<const mgis::behaviour::Behaviour>(b));
    static_cast<void>(bi);
  } catch (std::exception& e) {
    return e.what();
  }
  return {};
}

struct BehaviourIntegratorBaseTest final : public tfel::tests::TestCase {
  BehaviourIntegratorBaseTest()
      : tfel::tests::TestCase("MFEMMGIS", "BehaviourIntegratorBaseTest") {
  }  // end of BehaviourIntegratorBaseTest
  tfel::tests::TestResult execute() override {
    using namespace mfem_mgis;
    auto ctx = Context{};
    const auto ofed = construct<FiniteElementDiscretization>(
        ctx, getFiniteElementDiscretizationParameters(3));
    TFEL_TESTS_ASSERT(isValid(ofed));
    this->test1(*ofed);
    this->test2(*ofed);
    // scalar unknowns for the heat transfer integrator
    const auto ofed2 = construct<FiniteElementDiscretization>(
        ctx, getFiniteElementDiscretizationParameters(1));
    TFEL_TESTS_ASSERT(isValid(ofed2));
    this->test3(*ofed2);
    return this->result;
  }

 private:
  /*!
   * \return the parameters of the finite element discretization
   * \param[in] n: number of components of the unknowns
   */
  static mfem_mgis::Parameters getFiniteElementDiscretizationParameters(
      const int n) {
    return {{"MeshFileName", parameters.mesh_file},
            {"FiniteElementFamily", "H1"},
            {"FiniteElementOrder", 1},
            {"UnknownsSize", n},
            {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
            {"Parallel", bool(parameters.parallel)}};
  }
  //! \brief a finite strain integrator requires a PK1 behaviour
  void test1(const mfem_mgis::FiniteElementDiscretization& fed) {
    using namespace mgis::behaviour;
    using Integrator = mfem_mgis::
        IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator;
    const auto h = Hypothesis::TRIDIMENSIONAL;
    // default options: Cauchy stress and its derivative with respect to F
    const auto b1 = load(FiniteStrainBehaviourOptions{}, parameters.library,
                         "SaintVenantKirchhoffElasticity", h);
    TFEL_TESTS_CHECK(contains(getConstructionError<Integrator>(fed, b1),
                              "first Piola-Kirchhoff stress"));
    // options used by mfem-mgis
    auto opts = FiniteStrainBehaviourOptions{};
    opts.stress_measure = FiniteStrainBehaviourOptions::PK1;
    opts.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
    const auto b2 =
        load(opts, parameters.library, "SaintVenantKirchhoffElasticity", h);
    TFEL_TESTS_CHECK(getConstructionError<Integrator>(fed, b2).empty());
    // small strain behaviour
    const auto b3 = load(parameters.library, "Elasticity", h);
    TFEL_TESTS_CHECK(contains(getConstructionError<Integrator>(fed, b3),
                              "invalid behaviour type"));
  }
  //! \brief the hypothesis of the behaviour must match the integrator one
  void test2(const mfem_mgis::FiniteElementDiscretization& fed) {
    using namespace mgis::behaviour;
    using Integrator = mfem_mgis::
        IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator;
    const auto b =
        load(parameters.library, "Elasticity", Hypothesis::PLANESTRAIN);
    TFEL_TESTS_CHECK(
        contains(getConstructionError<Integrator>(fed, b), "does not match"));
  }
  //! \brief the transient heat transfer integrator requires an isotropic
  //! behaviour
  void test3(const mfem_mgis::FiniteElementDiscretization& fed) {
    using namespace mgis::behaviour;
    using Integrator = mfem_mgis::TransientHeatTransferBehaviourIntegrator;
    const auto b = load(parameters.library, "OrthotropicElasticity",
                        Hypothesis::TRIDIMENSIONAL);
    TFEL_TESTS_CHECK(contains(getConstructionError<Integrator>(fed, b),
                              "invalid behaviour symmetry"));
  }
  //! \return if the string s contains the string s2
  static bool contains(const std::string& s, std::string_view s2) {
    return s.find(s2) != std::string::npos;
  }
};

TFEL_TESTS_GENERATE_PROXY(BehaviourIntegratorBaseTest,
                          "BehaviourIntegratorBaseTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.library, "-l", "--library", "Material library.");
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
    m.addXMLTestOutput("ParallelBehaviourIntegratorBaseTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("BehaviourIntegratorBaseTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
