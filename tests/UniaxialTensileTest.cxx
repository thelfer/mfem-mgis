/*!
 * \file   tests/UniaxialTensileTest.cxx
 * \brief
 * This example models a cyclic tension-compression test on a unit cube. The
 * displacement imposed on the face x = 1, equal to the axial strain, goes up
 * to 0.9 %, down to -2.1 % and back up to 1.9 %. The faces x = 0, y = 0 and
 * z = 0 are symmetry planes.
 *
 * At each time step, the first two components of the gradients, the first
 * component of the thermodynamic forces and an internal state variable are
 * extracted at the first integration point. They are saved in the file
 * UniaxialTensileTest-<behaviour>.txt and compared to the values of the
 * reference file, if any.
 *
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <array>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string_view>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct TestParameters {
  const char* mesh_file = "cube.mesh";
  const char* library = "src/libBehaviour.so";
  const char* behaviour = "IsotropicLinearHardeningPlasticity";
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* isv_name = "EquivalentPlasticStrain";
  const char* reference_file = "";
  const char* linearsolver = "CGSolver";
  int order = 1;
  int refinement = 0;
  int nbsteps = 100;
#ifdef MFEM_USE_MPI
  bool parallel = true;
#else
  bool parallel = false;
#endif
  bool post_processing = true;
};

//! \brief values extracted at the first integration point at each time step
struct UniaxialTestResults {
  //! \brief first component of the gradients in the material frame
  std::vector<mfem_mgis::real> g0;
  //! \brief second component of the gradients in the material frame
  std::vector<mfem_mgis::real> g1;
  //! \brief first component of the thermodynamic forces in the material frame
  std::vector<mfem_mgis::real> tf0;
  //! \brief selected internal state variable, if any
  std::vector<mfem_mgis::real> v;
};

static TestParameters parseCommandLineOptions(int argc, char** argv) {
  auto p = TestParameters{};
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.behaviour, "-b", "--behaviour", "Name of the behaviour.");
  args.AddOption(&p.isv_name, "-isv", "--internal-state-variable",
                 "Internal state variable saved and compared to the "
                 "reference file, none if empty.");
  args.AddOption(&p.reference_file, "-rf", "--reference-file",
                 "Reference file, no comparison if empty.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&p.nbsteps, "-ns", "--nbsteps",
                 "Number of time steps. The reference file must have been "
                 "computed with the same number of time steps.");
  args.AddOption(&p.linearsolver, "-ls", "--linearsolver",
                 "Linear solver: GMRESSolver, CGSolver, UMFPackSolver "
                 "(sequential only), MUMPSSolver, HypreFGMRES with the "
                 "HypreILU preconditioner, HyprePCG with the HypreDiagScale "
                 "preconditioner or HypreGMRES (parallel only).");
  args.AddOption(&p.parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.AddOption(&p.post_processing, "-pp", "--post-processing", "-no-pp",
                 "--no-post-processing",
                 "Export or not the results to Paraview.");
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    std::exit(EXIT_SUCCESS);
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
  if (p.nbsteps <= 0) {
    mfem_mgis::getErrorStream() << "Invalid number of time steps\n";
    mfem_mgis::abort(EXIT_FAILURE);
  }
  return p;
}

[[nodiscard]] static bool setLinearSolver(
    mfem_mgis::Context& ctx,
    mfem_mgis::AbstractNonLinearEvolutionProblem& p,
    const std::string_view s) noexcept {
  const auto ilu = mfem_mgis::Parameters{
      {"Name", "HypreILU"},
      {"Options", mfem_mgis::Parameters{{"HypreILULevelOfFill", 1}}}};
  const auto diagscale = mfem_mgis::Parameters{
      {"Name", "HypreDiagScale"},
      {"Options", mfem_mgis::Parameters{{"VerbosityLevel", 0}}}};
  if (s == "GMRESSolver") {
    return p.setLinearSolver(ctx, "GMRESSolver",
                             {{"VerbosityLevel", 1},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "CGSolver") {
    return p.setLinearSolver(ctx, "CGSolver",
                             {{"VerbosityLevel", 1},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
#ifdef MFEM_USE_SUITESPARSE
  } else if (s == "UMFPackSolver") {
    return p.setLinearSolver(ctx, "UMFPackSolver", {});
#endif
#ifdef MFEM_USE_MUMPS
  } else if (s == "MUMPSSolver") {
    return p.setLinearSolver(ctx, "MUMPSSolver", {{"Symmetric", true}});
#endif
  } else if (s == "HypreFGMRES") {
    return p.setLinearSolver(ctx, "HypreFGMRES",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"Preconditioner", ilu},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "HyprePCG") {
    return p.setLinearSolver(ctx, "HyprePCG",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"Preconditioner", diagscale},
                              {"MaximumNumberOfIterations", 5000}});
  } else if (s == "HypreGMRES") {
    return p.setLinearSolver(ctx, "HypreGMRES",
                             {{"VerbosityLevel", 1},
                              {"Tolerance", 1e-12},
                              {"MaximumNumberOfIterations", 5000}});
  }
  return ctx.registerErrorMessage("unsupported linear solver '" +
                                  std::string{s} + "'");
}

/*!
 * \brief add the values of the given state at the first integration point, if
 * any
 * \param[in, out] r: results
 * \param[in] m: material
 * \param[in] s: state of the material
 * \param[in] p: parameters
 */
static void extractResults(UniaxialTestResults& r,
                           const mfem_mgis::Material& m,
                           const mgis::behaviour::MaterialStateManager& s,
                           const TestParameters& p) {
  if (m.n == 0) {
    return;
  }
  r.g0.push_back(s.gradients[0]);
  r.g1.push_back(s.gradients[1]);
  r.tf0.push_back(s.thermodynamic_forces[0]);
  if (std::string_view{p.isv_name}.empty()) {
    r.v.push_back(0);
  } else {
    const auto o = mgis::behaviour::getVariableOffset(m.b.isvs, p.isv_name,
                                                      m.b.hypothesis);
    r.v.push_back(s.internal_state_variables[o]);
  }
}

/*!
 * \return the values of a file of results, or empty results if the file
 * can't be read
 * \param[in] f: file name
 */
static UniaxialTestResults readResults(const std::string& f) {
  auto r = UniaxialTestResults{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    auto values = std::array<mfem_mgis::real, 4>{};
    if (!(std::istringstream(line) >> values[0] >> values[1] >> values[2] >>
          values[3])) {
      mfem_mgis::getErrorStream()
          << "invalid line '" << line << "' in '" << f << "'\n";
      return {};
    }
    r.g0.push_back(values[0]);
    r.g1.push_back(values[1]);
    r.tf0.push_back(values[2]);
    r.v.push_back(values[3]);
  }
  return r;
}

/*!
 * \return if the results match the reference values
 * \param[in] r: results
 * \param[in] m: material
 * \param[in] p: parameters
 */
static bool checkResults(const UniaxialTestResults& r,
                         const mfem_mgis::Material& m,
                         const TestParameters& p) {
  // absolute tolerances on the gradients, the internal state variable and the
  // thermodynamic forces
  constexpr auto eeps = mfem_mgis::real(1.e-10);
  constexpr auto seps = mfem_mgis::real(70.e9) * eeps;
  auto success = true;
  if ((m.n != 0) && (!std::string_view{p.reference_file}.empty())) {
    const auto references = readResults(p.reference_file);
    if (references.g0.size() != r.g0.size()) {
      mfem_mgis::getErrorStream()
          << "test failed ('" << p.reference_file << "' has "
          << references.g0.size() << " values instead of " << r.g0.size()
          << ")\n";
      success = false;
    } else {
      auto check = [&success](const auto i, const auto cv, const auto rv,
                              const auto eps, const auto msg) {
        const auto e = std::abs(cv - rv);
        if (e > eps) {
          mfem_mgis::getErrorStream()
              << "test failed (" << msg << " at time step " << i << ", " << cv
              << " vs " << rv << ", error " << e << ")\n";
          success = false;
        }
      };
      for (std::size_t i = 0; i != r.g0.size(); ++i) {
        check(i, r.g1[i], references.g1[i], eeps,
              "invalid transverse gradients");
        check(i, r.tf0[i], references.tf0[i], seps,
              "invalid thermodynamic force value");
        check(i, r.v[i], references.v[i], eeps,
              "invalid internal state variable");
      }
    }
  }
#ifdef MFEM_USE_MPI
  const auto& fed =
      m.getPartialQuadratureSpace().getFiniteElementDiscretization();
  if (fed.describesAParallelComputation()) {
    MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_CXX_BOOL, MPI_LAND,
                  mfem_mgis::getMPICommunicator(fed));
  }
#endif /* MFEM_USE_MPI */
  return success;
}

/*!
 * \brief save the results. Only the first process holding an integration
 * point writes the file.
 * \param[in] f: file name
 * \param[in] r: results
 * \param[in] m: material
 */
static void saveResults(const std::string& f,
                        const UniaxialTestResults& r,
                        const mfem_mgis::Material& m) {
  auto writer = (m.n != 0);
#ifdef MFEM_USE_MPI
  const auto& fed =
      m.getPartialQuadratureSpace().getFiniteElementDiscretization();
  if (fed.describesAParallelComputation()) {
    const auto comm = mfem_mgis::getMPICommunicator(fed);
    auto rank = int{};
    auto size = int{};
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    auto first = writer ? rank : size;
    MPI_Allreduce(MPI_IN_PLACE, &first, 1, MPI_INT, MPI_MIN, comm);
    writer = (rank == first);
  }
#endif /* MFEM_USE_MPI */
  if (!writer) {
    return;
  }
  std::ofstream out(f);
  out.precision(14);
  for (std::size_t i = 0; i != r.g0.size(); ++i) {
    out << r.g0[i] << " " << r.g1[i] << " " << r.tf0[i] << " " << r.v[i]
        << '\n';
  }
}

int main(int argc, char** argv) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  mfem_mgis::initialize(argc, argv);
  const auto p = parseCommandLineOptions(argc, argv);
  constexpr const auto dim = mfem_mgis::size_type{3};
  // building the non linear problem
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", p.mesh_file},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", p.order},
              {"UnknownsSize", dim},
              {"NumberOfUniformRefinements", p.refinement},
              {"Materials", mfem_mgis::Parameters{{"Cube", 1}}},
              {"Boundaries",
               mfem_mgis::Parameters{
                   {"Ymin", 1}, {"Zmin", 2}, {"Xmax", 3}, {"Xmin", 5}}},
              {"Hypothesis", "Tridimensional"},
              {"Parallel", p.parallel}}) |
      or_die;
  // material
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", p.library,
                                 p.behaviour) |
      or_die;
  auto& m1 = problem.getMaterial(ctx, "Cube", 0) | or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15) |
      or_die;
  if (m1.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
    std::array<mfem_mgis::real, 9u> r = {0, 1, 0,  //
                                         1, 0, 0,  //
                                         0, 0, 1};
    m1.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
  }
  // boundary conditions: symmetry planes and imposed axial displacement
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Xmin"}, {"Component", 0}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Ymin"}, {"Component", 1}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Zmin"}, {"Component", 2}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "Xmax"},
            {"Component", 0},
            {"LoadingEvolution",
             [](const auto t) {
               if (t < 0.3) {
                 return 3e-2 * t;
               } else if (t < 0.6) {
                 return 0.009 - 0.1 * (t - 0.3);
               }
               return -0.021 + 0.1 * (t - 0.6);
             }}}) |
      or_die;
  // solver parameters
  setLinearSolver(ctx, problem, p.linearsolver) | or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // post-processings
  if (p.post_processing) {
    problem.addPostProcessing(
        ctx, "ParaviewExportResults",
        {{"OutputFileName",
          "UniaxialTensileTestOutput-" + std::string(p.behaviour)}}) |
        or_die;
    const auto& b = m1.b;
    if ((b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) &&
        (b.kinematic == mgis::behaviour::Behaviour::SMALLSTRAINKINEMATIC)) {
      problem.addPostProcessing(
          ctx, "ParaviewExportIntegrationPointResultsAtNodes",
          {{"OutputFileName", "UniaxialTensileTestIntegrationPointOutput-" +
                                  std::string(p.behaviour)},
           {"Materials", {"Cube"}},
           {"Results", {"Strain"}}}) |
          or_die;
    }
  }
  // loop over the time steps
  const auto dt = mfem_mgis::real{1} / p.nbsteps;
  auto t = mfem_mgis::real{0};
  auto r = UniaxialTestResults{};
  extractResults(r, m1, m1.s0, p);
  for (int i = 0; i != p.nbsteps; ++i) {
    mfem_mgis::getOutputStream()
        << "time step " << i + 1 << " from " << t << " to " << t + dt << '\n';
    problem.solve(ctx, t, dt) | or_die;
    problem.executePostProcessings(ctx, t, dt) | or_die;
    problem.update(ctx) | or_die;
    t += dt;
    extractResults(r, m1, m1.s1, p);
  }
  // save and compare the results
  saveResults("UniaxialTensileTest-" + std::string(p.behaviour) + ".txt", r,
              m1);
  const auto success = checkResults(r, m1, p);
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
