/*!
 * \file   src/LinearSolverFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <utility>
#include "mfem/linalg/solvers.hpp"
#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  std::unique_ptr<LinearSolverPreconditioner> setHypreBoomerAMGPreconditioner(
      mfem::IterativeSolver& s,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& opts) {
    auto amg = std::make_unique<mfem::HypreBoomerAMG>();
    checkParameters(opts, {"Strategy"});
    if (contains(opts, "Strategy")) {
      const auto strategy = get<std::string>("Strategy");
      auto& fespace = p.getFiniteElementSpace();
      if (strategy == "Elasticity") {
        amg->SetElasticityOptions(&fespace);
      } else if (strategy == "System") {
        const auto* const m = fespace.GetMesh();
        const auto o = fespace.GetOrdering();
        amg->SetSystemsOptions(m->Dimension(), o == mfem::Ordering::byNODES);
      } else {
        raise(
            "setLinearSolverParameters: "
            "invalid strategy '" +
            strategy + "' for preconditioner HypreBoomerAMG");
      }
    }
    s.SetPreconditioner(*amg);
    return std::move(amg);
  }  // end of setHypreBoomerAMGPreconditioner

#else /* MFEM_USE_MPI */

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreBoomerAMGPreconditioner(
      mfem::IterativeSolver&,
      NonLinearEvolutionProblemImplementation<true>&,
      const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreBoomerAMGPreconditioner

#endif /* MFEM_USE_MPI */

  template <bool parallel>
  std::unique_ptr<LinearSolverPreconditioner> setLinearSolverParameters(
      mfem::IterativeSolver& s,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params) {
    const char* const Preconditioner = "Preconditioner";
    s.iterative_mode = false;
    auto allowed_parameters = getIterativeSolverParametersList();
    allowed_parameters.push_back(Preconditioner);
    checkParameters(params, allowed_parameters);
    setSolverParameters(s, extract(params, getIterativeSolverParametersList()));
    if (contains(params, Preconditioner)) {
      const auto pr = get<Parameters>(params, Preconditioner);
      checkParameters(params, {"Name", "Options"});
      const auto name = get<std::string>(pr, "Name");
      if (name == "HypreBoomerAMG") {
        if constexpr (parallel) {
          return setHypreBoomerAMGPreconditioner(
              s, p, get_if<Parameters>(params, "Options", Parameters{}));
        } else {
          raise(
              "setLinearSolverParameters: "
              "the 'HypreBoomerAMG' is only available in parallel");
        }
      } else {
        raise(
            "setLinearSolverParameters: "
            "unsupported preconditioner '" +
            name + "'");
      }
    }
    return {};
  }  // end of setLinearSolverParameters

  template <bool parallel, typename LinearSolverType>
  std::function<LinearSolverHandler(
      NonLinearEvolutionProblemImplementation<parallel>&, const Parameters&)>
  buildIterativeSolverGenerator() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return [](NonLinearEvolutionProblemImplementation<true>& p,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>(MPI_COMM_WORLD);
        auto prec = setLinearSolverParameters(*s, p, params);
        return LinearSolverHandler{std::move(s), std::move(prec)};
      };
#else  /* MFEM_USE_MPI */
      return {};
#endif /* MFEM_USE_MPI */
    } else {
      return [](NonLinearEvolutionProblemImplementation<false>& p,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>();
        auto prec = setLinearSolverParameters(*s, p, params);
        return LinearSolverHandler{std::move(s), std::move(prec)};
      };
    }
  }  // end of buildIterativeSolverGenerator

#ifdef MFEM_USE_MUMPS

  std::function<LinearSolverHandler(
      NonLinearEvolutionProblemImplementation<true>&, const Parameters&)>
  buildMUMPSSolverGenerator() {
    return [](NonLinearEvolutionProblemImplementation<true>&,
              const Parameters& params) {
      checkParameters(params, {"Symmetric", "PositiveDefinite"});
      auto s = std::make_unique<mfem::MUMPSSolver>();
      const auto symmetric = get_if<bool>(params, "Symmetric", false);
      const auto positive_definite =
          get_if<bool>(params, "PositiveDefinite", false);
      if (symmetric) {
        if (positive_definite) {
          s->SetMatrixSymType(
              mfem::MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
        } else {
          s->SetMatrixSymType(mfem::MUMPSSolver::MatType::SYMMETRIC_INDEFINITE);
        }
      } else {
        s->SetMatrixSymType(mfem::MUMPSSolver::MatType::UNSYMMETRIC);
      }
      return {s, std::unique_ptr<LinearSolverPreconditioner>{}};
    };
  }  // end of builMUMPSGenerator

#endif /* MFEM_USE_MUMPS */

  template <bool parallel>
  static void declareDefaultSolvers(LinearSolverFactory<parallel>& f) {
    f.add("CGSolver",
          buildIterativeSolverGenerator<parallel, mfem::CGSolver>());
    f.add("GMRESSolver",
          buildIterativeSolverGenerator<parallel, mfem::GMRESSolver>());
    if constexpr (parallel) {
#ifdef MFEM_USE_MUMPS
      f.add("MUMPSSolver", buildMUMPSSolverGenerator());
#endif
    }
    if constexpr (!parallel) {
#ifdef MFEM_USE_SUITESPARSE
      f.add("UMFPackSolver", [](NonLinearEvolutionProblemImplementation<false>&,
                                const Parameters&) {
        return LinearSolverHandler{
            std::make_unique<mfem::UMFPackSolver>(),
            std::unique_ptr<LinearSolverPreconditioner>{}};
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
      raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  LinearSolverHandler LinearSolverFactory<true>::generate(
      std::string_view n,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<true>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared";
      raise(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(p, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "error while generating no linear '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      raise(msg);
    }
    return s;
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
      raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  LinearSolverHandler LinearSolverFactory<false>::generate(
      std::string_view n,
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared";
      raise(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(p, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "error while generating no linear '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      raise(msg);
    }
    return s;
  }  // end of generate

  LinearSolverFactory<false>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<false>::~LinearSolverFactory() = default;

}  // end of namespace mfem_mgis
