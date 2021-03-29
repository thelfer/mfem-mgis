/*!
 * \file   src/LinearSolverFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <utility>
#include "mfem/linalg/solvers.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

namespace mfem_mgis {

  static void setLinearSolverParameters(mfem::IterativeSolver& s,
                                        const Parameters& params) {
    s.iterative_mode = false;
    using Problem = AbstractNonLinearEvolutionProblem;
    if (contains(params, Problem::SolverVerbosityLevel)) {
      s.SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
    }
    if (contains(params, Problem::SolverRelativeTolerance)) {
      s.SetRelTol(get<double>(params, Problem::SolverRelativeTolerance));
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      s.SetAbsTol(get<double>(params, Problem::SolverAbsoluteTolerance));
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      s.SetMaxIter(get<int>(params, Problem::SolverMaximumNumberOfIterations));
    }
  }  // end of setLinearSolverParameters

  template <bool parallel, typename LinearSolverType>
  std::function<std::unique_ptr<LinearSolver>(const Parameters&)>
  buildIterativeSolverGenerator() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return [](const Parameters& p) {
        auto s = std::make_unique<LinearSolverType>(MPI_COMM_WORLD);
        setLinearSolverParameters(*s, p);
        return s;
      };
#else  /* MFEM_USE_MPI */
      return {};
#endif /* MFEM_USE_MPI */
    } else {
      return [](const Parameters& p) {
        auto s = std::make_unique<LinearSolverType>();
        setLinearSolverParameters(*s, p);
        return s;
      };
    }
  }  // end of buildIterativeSolverGenerator

  template <bool parallel>
  static void declareDefaultSolvers(LinearSolverFactory<parallel>& f) {
    f.add("CGSolver",
          buildIterativeSolverGenerator<parallel, mfem::CGSolver>());
    f.add("GMRESSolver",
          buildIterativeSolverGenerator<parallel, mfem::GMRESSolver>());
    if constexpr (!parallel) {
#ifdef MFEM_USE_SUITESPARSE
      f.add("UMFPackSolver", [](const Parameters&) {
        return std::make_unique<mfem::UMFPackSolver>();
      });
#endif
    }
  }  // end of declareDefaultSolvers

#ifdef MFEM_USE_MPI

  LinearSolverFactory<true>& LinearSolverFactory<true>::getFactory() {
    static LinearSolverFactory<true> factory;
    return factory;
  }  // end of getFactory

  void LinearSolverFactory<true>::add(std::string_view n, Generator g) {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      std::string msg("LinearSolverFactory<true>::add: ");
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      mgis::raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<LinearSolver> LinearSolverFactory<true>::generate(
      std::string_view n, const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<true>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared";
      mgis::raise(msg);
    }
    const auto& g = pg->second;
    return g(params);
  }  // end of generate

  LinearSolverFactory<true>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<true>::~LinearSolverFactory() = default;

#endif /* MFEM_USE_MPI */

  LinearSolverFactory<false>& LinearSolverFactory<false>::getFactory() {
    static LinearSolverFactory<false> factory;
    return factory;
  }  // end of getFactory

  void LinearSolverFactory<false>::add(std::string_view n, Generator g) {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      std::string msg("LinearSolverFactory<false>::add: ");
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      mgis::raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<LinearSolver> LinearSolverFactory<false>::generate(
      std::string_view n, const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared";
      mgis::raise(msg);
    }
    const auto& g = pg->second;
    return g(params);
  }  // end of generate

  LinearSolverFactory<false>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<false>::~LinearSolverFactory() = default;

}  // end of namespace mfem_mgis
