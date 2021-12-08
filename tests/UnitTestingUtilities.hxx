/*!
 * \file   UnitTestingUtilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   08/04/2021
 */

#ifndef LIB_MFEM_MGIS_UNITTESTINGUTILITIES_HXX
#define LIB_MFEM_MGIS_UNITTESTINGUTILITIES_HXX

#include <vector>
#include <fstream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis::unit_tests {

  struct TestParameters {
    const char* mesh_file = nullptr;
    const char* behaviour = nullptr;
    const char* library = nullptr;
    const char* reference_file = nullptr;
    const char* isv_name = nullptr;
    int linearsolver = 0;
    int order = 1;
  };  // end of struct TestParameters

  struct UniaxialTestResults {
    /*!
     * \brief values of the first component of the gradients in the material
     * frame.
     */
    std::vector<mfem_mgis::real> g0;
    /*!
     * \brief values of the second component of the gradients in the material
     * frame.
     */
    std::vector<mfem_mgis::real> g1;
    /*!
     * \brief values of the first component of the thermodynamic forces in the
     * material frame.
     */
    std::vector<mfem_mgis::real> tf0;
    //! \brief values of the selected internal state variable, if any
    std::vector<mfem_mgis::real> v;
  };  // end of struct UniaxialTestResults

  [[maybe_unused]] static void parseCommandLineOptions(TestParameters& params,
                                                       int argc,
                                                       char** argv) {
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&params.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&params.reference_file, "-r", "--reference-file",
                   "Reference file.");
    args.AddOption(&params.behaviour, "-b", "--behaviour",
                   "Name of the behaviour.");
    args.AddOption(&params.isv_name, "-v", "--internal-state-variable",
                   "Internal variable name to be post-processed.");
    args.AddOption(&params.library, "-l", "--library", "Material library.");
    args.AddOption(
        &params.linearsolver, "-ls", "--linearsolver",
        "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
    args.AddOption(&params.order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.Parse();
    if ((!args.Good()) || (params.mesh_file == nullptr) ||
        (params.library == nullptr) || (params.behaviour == nullptr)) {
      args.PrintUsage(mfem_mgis::getOutputStream());
      mfem_mgis::abort(EXIT_FAILURE);
    }
    // args.PrintOptions(mfem_mgis::getOutputStream());
  }  // end of parseCommandLineOptions

  [[maybe_unused]] static void setLinearSolver(
      mfem_mgis::NonLinearEvolutionProblem& problem,
      const TestParameters& parameters) {
    if (parameters.linearsolver == 0) {
      problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                           {"AbsoluteTolerance", 1e-12},
                                           {"RelativeTolerance", 1e-12},
                                           {"MaximumNumberOfIterations", 300}});
    } else if (parameters.linearsolver == 1) {
      problem.setLinearSolver("GMRESSolver",
                              {{"VerbosityLevel", 1},
                               {"AbsoluteTolerance", 1e-12},
                               {"RelativeTolerance", 1e-12},
                               {"MaximumNumberOfIterations", 300}});
#ifdef MFEM_USE_SUITESPARSE
    } else if (parameters.linearsolver == 2) {
      problem.setLinearSolver("UMFPackSolver", {});
#endif
#ifdef MFEM_USE_MUMPS
    } else if (parameters.linearsolver == 3) {
      problem.setLinearSolver("MUMPSSolver", {{"Symmetric", true}});
#endif
    } else {
      mfem_mgis::getErrorStream() << "unsupported linear solver\n";
      mfem_mgis::abort(EXIT_FAILURE);
    }
  }  // end of setLinearSolver

  [[maybe_unused]] static void extractResults(
      UniaxialTestResults& r,
      const mfem_mgis::Material& m,
      const mgis::behaviour::MaterialStateManager& s,
      const TestParameters& parameters) {
    if (m.n != 0) {
      r.g0.push_back(s.gradients[0]);
      r.g1.push_back(s.gradients[1]);
      r.tf0.push_back(s.thermodynamic_forces[0]);
      if (parameters.isv_name != nullptr) {
        const auto vo = mgis::behaviour::getVariableOffset(
            m.b.isvs, parameters.isv_name, m.b.hypothesis);
        r.v.push_back(s.internal_state_variables[vo]);
      } else {
        r.v.push_back(0);
      }
    }
  }  // end of extractResults

  [[maybe_unused]] static void extractInitialResults(
      UniaxialTestResults& r,
      const mfem_mgis::Material& m,
      const TestParameters& parameters) {
    extractResults(r, m, m.s0, parameters);
  }  // end of extractInitialResults

  [[maybe_unused]] static void extractResults(
      UniaxialTestResults& r,
      const mfem_mgis::Material& m,
      const TestParameters& parameters) {
    extractResults(r, m, m.s1, parameters);
  }  // end of extractResults

  [[maybe_unused]] static bool checkResults(UniaxialTestResults& r,
                                            const mfem_mgis::Material& m,
                                            const TestParameters& parameters,
                                            const real eeps,
                                            const real seps) {
    // comparison to reference results
    bool success = true;
    if ((m.n != 0) && (parameters.reference_file != nullptr)) {
      std::ifstream in(parameters.reference_file);
      if (in) {
        auto check = [&success](const auto cv,  // computed value
                                const auto rv,  // reference value,
                                const auto ev, const auto msg) {
          const auto e = std::abs(cv - rv);
          if (e > ev) {
            mfem_mgis::getErrorStream()
                << "test failed (" << msg << ", " << cv << " vs " << rv
                << ", error " << e << ")\n";
            success = false;
          }
        };
        for (std::vector<mfem_mgis::real>::size_type i = 0; i != r.g0.size();
             ++i) {
          auto g0_ref = mfem_mgis::real{};
          auto g1_ref = mfem_mgis::real{};
          auto tf0_ref = mfem_mgis::real{};
          auto v_ref = mfem_mgis::real{};
          in >> g0_ref >> g1_ref >> tf0_ref >> v_ref;
          check(r.g1[i], g1_ref, eeps, "invalid transverse gradients");
          check(r.tf0[i], tf0_ref, seps, "invalid thermodynamic force value");
          check(r.v[i], v_ref, eeps, "invalid internal state variable");
        }
#ifdef MFEM_USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LAND,
                      MPI_COMM_WORLD);
#endif /* MFEM_USE_MPI */
      }
    }  // end of if (m1.n != 0)
    return success;
  }  // end of checkResults

  [[maybe_unused]] static void saveResults(const std::string& f,
                                           const UniaxialTestResults& r) {
    std::ofstream out(f);
    out.precision(14);
    for (std::vector<mfem_mgis::real>::size_type i = 0; i != r.g0.size(); ++i) {
      out << r.g0[i] << " " << r.g1[i] << " " << r.tf0[i] << " " << r.v[i]
          << '\n';
    }
  }

  [[maybe_unused]] static UniaxialTestResults solve(
      mfem_mgis::NonLinearEvolutionProblem& problem,
      const TestParameters& parameters,
      const mfem_mgis::real t0,
      const mfem_mgis::real t1,
      const mfem_mgis::size_type nsteps) {
    const auto& m1 = problem.getMaterial(1);
    const auto dt = (t1 - t0) / nsteps;
    auto r = mfem_mgis::unit_tests::UniaxialTestResults{};
    extractInitialResults(r, m1, parameters);
    // loop over time step
    auto t = t0;
    for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
      // resolution
      const auto step_timer = mfem_mgis::getTimer("step" + std::to_string(i));
      {
        const auto solve_timer = mfem_mgis::getTimer("solve");
        if (!problem.solve(t, dt)) {
          mfem_mgis::abort("non convergence");
        }
      }
      {
        const auto post_processing_timer =
            mfem_mgis::getTimer("post_processing");
        problem.executePostProcessings(t, dt);
      }
      {
        const auto update_timer = mfem_mgis::getTimer("update");
        problem.update();
      }
      t += dt;
      extractResults(r, m1, parameters);
    }
    return r;
  }

}  // end of namespace mfem_mgis::unit_tests

#endif /* LIB_MFEM_MGIS_UNITTESTINGUTILITIES_HXX */
