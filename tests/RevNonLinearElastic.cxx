#include <memory>
#include <cstdlib>
#include <iostream>

#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include <MFEMMGIS/Profiler.hxx>
#include <functional>


// We need this class for test case sources
struct TestParameters {
	const char* mesh_file = "cube_2mat_per.mesh";
	const char* behaviour = "SaintVenantKirchhoffElasticity";
	const char* library = "src/libBehaviour.so";
	int order = 1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	bool parallel = true;
	int refinement = 0;
	int post_processing = 1; // default value : disabled
	int verbosity_level = 1; // default value : lower level
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
	args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 0");
	args.AddOption(&p.post_processing, "-p", "--post-processing", "run post processing step");
	args.AddOption(&p.verbosity_level, "-v", "--verbosity-level", "choose the verbosity level");

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

	template<typename Problem>
void add_post_processings(Problem& p, std::string msg)
{
	p.addPostProcessing(
			"ParaviewExportResults",
			{{"OutputFileName", msg}}
			);
} // end timer add_postprocessing_and_outputs

	template<typename Problem>
void execute_post_processings(Problem& p, double start, double end)
{
	CatchTimeSection("common::post_processing_step");
	p.executePostProcessings(start, end);
}

void setup_properties(const TestParameters& p, mfem_mgis::PeriodicNonLinearEvolutionProblem& problem)
{
	using namespace mgis::behaviour;
	using real=mfem_mgis::real;

	CatchTimeSection("set_mgis_stuff");
	problem.addBehaviourIntegrator("Mechanics", 1, p.library, p.behaviour);
	problem.addBehaviourIntegrator("Mechanics", 2, p.library, p.behaviour);
	// materials
	auto& m1 = problem.getMaterial(1);
	auto& m2 = problem.getMaterial(2);
	auto set_properties = [](auto& m, const double yo, const double po) {
		setMaterialProperty(m.s0, "YoungModulus", yo);
		setMaterialProperty(m.s0, "PoissonRatio", po);
		setMaterialProperty(m.s1, "YoungModulus", yo);
		setMaterialProperty(m.s1, "PoissonRatio", po);
	};

	set_properties(m1, 2.0e11       , 0.3);
	set_properties(m2, 8.0e11       , 0.3);

	//
	auto set_temperature = [](auto& m) {
		setExternalStateVariable(m.s0, "Temperature", 293.15);
		setExternalStateVariable(m.s1, "Temperature", 293.15);
	};
	set_temperature(m1);
	set_temperature(m2);



	// macroscopic strain
	std::vector<real> e(9, real{0});
	e[0] = 1.1;
	e[1] = 1.0;
	e[2] = 1.0;
	problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
} 


	template<typename Problem>		
static void setLinearSolver(Problem& p,
		const int verbosity = 0,
		const mfem_mgis::real Tol = 1e-12
		)
{
	CatchTimeSection("set_linear_solver");
	// pilote
	constexpr int defaultMaxNumOfIt	 	= 5000; 		// MaximumNumberOfIterations
	constexpr int adjustMaxNumOfIt 		= 500000; 		// MaximumNumberOfIterations
	auto solverParameters = mfem_mgis::Parameters{};
	solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
	solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
	solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});


	// preconditionner hypreBoomerAMG
	auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	p.setLinearSolver("HyprePCG", solverParameters);
}

	template<typename Problem>
void run_solve(Problem& p, double start, double end)
{
	CatchTimeSection("Solve");
	// solving the problem
	auto statistics = p.solve(0, 1);

	// check status
	if (statistics.status) {
		mfem_mgis::Profiler::Utils::Message("INFO: FAILED");
	}
}

int main(int argc, char* argv[]) 
{
	// mpi initialization here 
	mfem_mgis::initialize(argc, argv);

	// init timers
	mfem_mgis::Profiler::timers::init_timers();

	// get parameters
	TestParameters p;
	mfem::OptionsParser args(argc, argv);
	common_parameters(args, p);

	// add post processing
	const bool use_post_processing = (p.post_processing == 1);

	// 3D
	constexpr const auto dim = mfem_mgis::size_type{3};

	// creating the finite element workspace
	auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
			mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
			{"FiniteElementFamily", "H1"},
			{"FiniteElementOrder", p.order},
			{"UnknownsSize", dim},
			{"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
			{"Parallel", p.parallel}});
	mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

	// set problem
	setup_properties(p, problem);
	setLinearSolver(problem, p.verbosity_level);

	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-rev-non-linear-elastic");

	// main function here
	run_solve(problem, 0, 1);

	if(use_post_processing)	execute_post_processings(problem, 0,1);

	// print and write timetable
	mfem_mgis::Profiler::timers::print_and_write_timers();
	return(EXIT_SUCCESS);
}
