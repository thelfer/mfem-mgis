/*!
 * \file   SatohTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2022
 */

#include <memory>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include <MFEMMGIS/Profiler.hxx>


template <bool parallel>
static void dumpPartialQuadratureFunction(
    std::ostream& os,
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& f) {
  const auto& s = f.getPartialQuadratureSpace();
  const auto& fed = s.getFiniteElementDiscretization();
  const auto& fespace = fed.getFiniteElementSpace<parallel>();
  const auto m = s.getId();
  for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
    if (fespace.GetAttribute(i) != m) {
      continue;
    }
    const auto& fe = *(fespace.GetFE(i));
    auto& tr = *(fespace.GetElementTransformation(i));
    const auto& ir = s.getIntegrationRule(fe, tr);
    for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
      // get the gradients of the shape functions
      mfem::Vector p;
      const auto& ip = ir.IntPoint(g);
      tr.SetIntPoint(&ip);
      tr.Transform(tr.GetIntPoint(), p);
      os << p[0] << " " << p[1] << " " << f.getIntegrationPointValue(i, g)
         << '\n';
    }
  }
}  // end of dumpPartialQuadratureFunction


struct TestParameters {
	const char* mesh_file = "inclusion.msh";
	const char* behaviour = "ElasticitySatohTest";
	const char* library = "src/libBehaviour.so";
	int order = 1;
	bool parallel = true;
	int refinement = 0;
	int verbosity_level = 0; // default value : lower level
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
	args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 0");
	args.AddOption(&p.verbosity_level, "-p", "--parallel", "choose between serial (-p 0 and parallel -p 1");

	args.Parse();

	if (!args.Good()) {
		if (mfem_mgis::getMPIrank() == 0)
			args.PrintUsage(std::cout);
		mfem_mgis::finalize();
		exit(0);
	}
	if (p.mesh_file == nullptr) {
		if (mfem_mgis::getMPIrank() == 0)
			std::cout << "ERROR: Mesh file missing" << std::endl;
		args.PrintUsage(std::cout);
		mfem_mgis::abort(EXIT_FAILURE);
	}
	if (mfem_mgis::getMPIrank() == 0)
		args.PrintOptions(std::cout);
}


/*
 * This test models a 2D plate of lenght 1 in plane strain clamped on the left
 * and right boundaries and submitted to a parabolic thermal gradient along the
 * x-axis:
 *
 * - the temperature profile is minimal on the left and right boundaries
 * - the temperature profile is maximal for x = 0.5
 *
 * This example shows how to define an external state variable using an
 * analytical profile.
 */
int main(int argc, char** argv) {
  // mpi initialization here 
  mfem_mgis::initialize(argc, argv);
  // init timers
  mfem_mgis::Profiler::timers::init_timers();

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);


  auto success = true;
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", p.mesh_file},
       {"Materials", mfem_mgis::Parameters{{"plate", 1}}},
       {"Boundaries", mfem_mgis::Parameters{{"left", 2}, {"right", 4}}},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", p.order},
       {"UnknownsSize", 2},
       {"NumberOfUniformRefinements", 2},
       {"Hypothesis", "PlaneStrain"},
       {"Parallel", p.parallel}});
  // materials
  problem.addBehaviourIntegrator("Mechanics", "plate", p.library, p.behaviour);
  auto& m1 = problem.getMaterial("plate");
  // material properties at the beginning and the end of the time step
  for (auto& s : {&m1.s0, &m1.s1}) {
    mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 150e9);
    mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.3);
  }
  // temperature
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  auto Tg = mfem_mgis::PartialQuadratureFunction::evaluate(
      m1.getPartialQuadratureSpacePointer(),
      [](const mfem_mgis::real x, const mfem_mgis::real) {
        constexpr auto Tref = mfem_mgis::real{293.15};
        constexpr auto dT = mfem_mgis::real{2000} - Tref;
        return 4 * dT * x * (1 - x) + Tref;
      });
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature",
                                            Tg->getValues());
  // boundary conditions
  for (const auto boundary : {"left", "right"}) {
    for (const auto dof : {0, 1}) {
      problem.addBoundaryCondition(
          std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
              problem.getFiniteElementDiscretizationPointer(), boundary, dof));
    }
  }
  // set the solver parameters
  problem.setLinearSolver("HyprePCG", {});
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "SatohTestOutput"}});
  auto results = std::vector<mfem_mgis::Parameter>{
      "Stress", "ImposedTemperature", "HydrostaticPressure"};
  problem.addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "SatohTestIntegrationPointOutput"},
       {"Materials", {"plate"}},
       {"Results", results}});
  // solving the problem on 1 time step
  auto r = problem.solve(0, 1);
  problem.executePostProcessings(0, 1);
  //
  std::ofstream output("HydrostaticPressure.txt");
  const auto pr = getInternalStateVariable(
      static_cast<const mfem_mgis::Material&>(m1), "HydrostaticPressure");
  if(p.parallel == 0)
  {
    dumpPartialQuadratureFunction<false>(output, pr);
  }
  else
  {
    dumpPartialQuadratureFunction<true>(output, pr);
  }
  //
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
