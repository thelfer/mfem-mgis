/*!
 * \file   tests/PartialQuadratureSpaceTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

auto parameters = mfem_mgis::unit_tests::TestParameters{};

struct PartialQuadratureSpaceTest final : public tfel::tests::TestCase {
  PartialQuadratureSpaceTest()
      : tfel::tests::TestCase("MFEMMGIS", "PartialQuadratureSpaceTest") {
  }  // end of PartialQuadratureSpaceTest

  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }

 private:
  void test1() {
    using namespace mfem_mgis;
    constexpr const auto dim = size_type{3};
    auto ctx = Context{};
    // building the non linear problem
    auto oproblem = construct<NonLinearEvolutionProblem>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", parameters.order},
                  {"UnknownsSize", dim},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Hypothesis", "Tridimensional"},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(oproblem));
    // materials
    const auto ok = oproblem->addBehaviourIntegrator(
        ctx, "Mechanics", 1, parameters.library, parameters.behaviour);
    TFEL_TESTS_ASSERT(isValid(ok));
    const auto& obi = oproblem->getBehaviourIntegrator(ctx, 1, 0);
    TFEL_TESTS_ASSERT(isValid(obi));
    const auto& qspace = obi->getPartialQuadratureSpace();
    TFEL_TESTS_CHECK(qspace.isDefinedOnAMaterial());
    TFEL_TESTS_CHECK(!qspace.isDefinedOnABoundary());
    const auto oqinfo = getInformation(ctx, qspace);
    TFEL_TESTS_ASSERT(isValid(oqinfo));
    if (parameters.parallel) {
      TFEL_TESTS_CHECK(oqinfo->identifier == 1);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells == 8);
      TFEL_TESTS_CHECK(oqinfo->number_of_quadrature_points == 8 * 27);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->first ==
          mfem::Geometry::CUBE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->second == 8);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->first == mfem::Geometry::CUBE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->second == 27);
    } else {
      TFEL_TESTS_CHECK(oqinfo->identifier == 1);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells == 1);
      TFEL_TESTS_CHECK(oqinfo->number_of_quadrature_points == 27);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->first ==
          mfem::Geometry::CUBE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->second == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->first == mfem::Geometry::CUBE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->second == 27);
    }
    TFEL_TESTS_CHECK(info(ctx, ctx.log(), *oqinfo));
  }  // end of test1

  void test2() {
    using namespace mfem_mgis;
    constexpr const auto dim = size_type{3};
    auto ctx = Context{};
    auto ofed = construct<FiniteElementDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", parameters.order},
                  {"UnknownsSize", dim},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Hypothesis", "Tridimensional"},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(ofed));
    auto oqspace = construct<PartialQuadratureSpace>(
        ctx, *ofed,
        LocationIdentifier{.material_identifier = {},
                           .boundary_identifier = BoundaryIdentifier{.id = 2}},
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation&) -> const mfem::IntegrationRule& {
          return mfem::IntRules.Get(e.GetGeomType(), 4);
        });
    TFEL_TESTS_ASSERT(isValid(oqspace));
    TFEL_TESTS_CHECK(oqspace->isDefinedOnABoundary());
    TFEL_TESTS_CHECK(!oqspace->isDefinedOnAMaterial());
    const auto oqinfo = getInformation(ctx, *oqspace);
    TFEL_TESTS_ASSERT(isValid(oqinfo));
    TFEL_TESTS_CHECK_EQUAL(oqinfo->identifier, 2);
    TFEL_TESTS_CHECK_EQUAL(oqinfo->name, "boundary (2)");
    if (parameters.parallel) {
      TFEL_TESTS_CHECK(oqinfo->number_of_cells == 4);
      TFEL_TESTS_CHECK(oqinfo->number_of_quadrature_points == 4 * 9);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->first ==
          mfem::Geometry::SQUARE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->second == 4);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->first == mfem::Geometry::SQUARE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->second == 9);
    } else {
      TFEL_TESTS_CHECK(oqinfo->number_of_cells == 1);
      TFEL_TESTS_CHECK(oqinfo->number_of_quadrature_points == 9);
      TFEL_TESTS_CHECK(oqinfo->number_of_cells_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->first ==
          mfem::Geometry::SQUARE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_cells_by_geometric_type.begin()->second == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->first == mfem::Geometry::SQUARE);
      TFEL_TESTS_CHECK(
          oqinfo->number_of_quadrature_points_by_geometric_type.begin()
              ->second == 9);
    }
  }  // end of test2
};

TFEL_TESTS_GENERATE_PROXY(PartialQuadratureSpaceTest,
                          "PartialQuadratureSpaceTest");

int main(int argc, char** argv) {
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);  // mfem_mgis::getDefaultLogStream());
  if (parameters.parallel) {
    m.addXMLTestOutput("ParallelPartialQuadratureSpaceTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("PartialQuadratureSpaceTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
