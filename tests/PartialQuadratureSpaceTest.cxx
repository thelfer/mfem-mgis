/*!
 * \file   tests/UniaxialTensileTest.cxx
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
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  auto ctx = mfem_mgis::Context{};
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", parameters.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", parameters.order},
       {"UnknownsSize", dim},
       {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", bool(parameters.parallel)}});
  // materials
  problem.addBehaviourIntegrator("Mechanics", 1, parameters.library,
                                 parameters.behaviour);
  const auto& qspace =
      problem.getBehaviourIntegrator(1).getPartialQuadratureSpace();
  const auto oqinfo = getInformation(ctx, qspace);
  if (isInvalid(oqinfo)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return false;
  }
  if (parameters.parallel) {
    std::cerr
        << "nqpoints: "
        << oqinfo->number_of_quadrature_points_by_geometric_type.begin()->second
        << '\n';
    assert(oqinfo->identifier == 1);
    assert(oqinfo->number_of_cells == 8);
    assert(oqinfo->number_of_quadrature_points == 8 * 27);
    assert(oqinfo->number_of_cells_by_geometric_type.size() == 1);
    assert(oqinfo->number_of_cells_by_geometric_type.begin()->first ==
           mfem::Geometry::CUBE);
    assert(oqinfo->number_of_cells_by_geometric_type.begin()->second == 8);
    assert(oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
    assert(
        oqinfo->number_of_quadrature_points_by_geometric_type.begin()->first ==
        mfem::Geometry::CUBE);
    assert(
        oqinfo->number_of_quadrature_points_by_geometric_type.begin()->second ==
        27);
  } else {
    assert(oqinfo->identifier == 1);
    assert(oqinfo->number_of_cells == 1);
    assert(oqinfo->number_of_quadrature_points == 27);
    assert(oqinfo->number_of_cells_by_geometric_type.size() == 1);
    assert(oqinfo->number_of_cells_by_geometric_type.begin()->first ==
           mfem::Geometry::CUBE);
    assert(oqinfo->number_of_cells_by_geometric_type.begin()->second == 1);
    assert(oqinfo->number_of_quadrature_points_by_geometric_type.size() == 1);
    assert(
        oqinfo->number_of_quadrature_points_by_geometric_type.begin()->first ==
        mfem::Geometry::CUBE);
    assert(
        oqinfo->number_of_quadrature_points_by_geometric_type.begin()->second ==
        27);
  }
  const auto success = mfem_mgis::info(ctx, std::cout, *oqinfo);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
