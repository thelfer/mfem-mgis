/*!
 * \file   src/LinearSolverFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <utility>
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/petsc.hpp"
#include "mfem/config/config.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif
#include "MGIS/Raise.hxx"
#include "MGIS/Contract.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  std::unique_ptr<LinearSolverPreconditioner> setHypreBoomerAMGPreconditioner(
      Context& ctx, FiniteElementSpace<true>& fespace, const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto amg = std::make_unique<mfem::HypreBoomerAMG>();
    checkParameters(opts, {"Strategy", Problem::SolverVerbosityLevel});
    if (contains(opts, "Strategy")) {
      const auto strategy = get<std::string>(opts, "Strategy");
      if (strategy == "Elasticity") {
        amg->SetElasticityOptions(&fespace);
      } else if (strategy == "System") {
        const auto* const m = fespace.GetMesh();
        const auto o = fespace.GetOrdering();
        amg->SetSystemsOptions(m->Dimension(), o == mfem::Ordering::byNODES);
      } else if (strategy != "None") {
        return ctx.registerErrorMessage(
            "setLinearSolverParameters: "
            "invalid strategy '" +
            strategy + "' for preconditioner HypreBoomerAMG");
      }
    }
    if (contains(opts, Problem::SolverVerbosityLevel)) {
      amg->SetPrintLevel(get<int>(opts, Problem::SolverVerbosityLevel));
    }
    return amg;
  }  // end of setHypreBoomerAMGPreconditioner

  std::unique_ptr<LinearSolverPreconditioner> setHypreEuclidPreconditioner(
      Context&, FiniteElementSpace<true>&, const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto euclid = std::make_unique<mfem::HypreEuclid>(MPI_COMM_WORLD);
    checkParameters(opts, {Problem::SolverVerbosityLevel});
    // if (contains(opts, Problem::SolverVerbosityLevel)) {
    // euclid->SetPrintLevel(get<int>(opts, Problem::SolverVerbosityLevel));
    //}
    return euclid;
  }  // end of setHypreEuclidPreconditioner

  std::unique_ptr<LinearSolverPreconditioner> setHypreILUPreconditioner(
      Context&, FiniteElementSpace<true>&, const Parameters& opts) {
#if MFEM_HYPRE_VERSION >= 21900
    using Problem = AbstractNonLinearEvolutionProblem;
    bool verbose = contains(opts, Problem::SolverVerbosityLevel);
    auto ilu = std::make_unique<mfem::HypreILU>();
    bool levelOfFill = contains(opts, "HypreILULevelOfFill");
    checkParameters(opts,
                    {Problem::SolverVerbosityLevel, "HypreILULevelOfFill"});
    if (verbose) {
      ilu->SetPrintLevel(get<int>(opts, Problem::SolverVerbosityLevel));
    }
    if (levelOfFill) {
      auto level = get<int>(opts, "HypreILULevelOfFill");
      HYPRE_ILUSetLevelOfFill(*ilu, level);
    }
    return ilu;
#else  /*  HYPRE_OLD_VERSION */
    MFEM_VERIFY(
        0, "Support for HypreILU is notavailable with this version of MFEM");
    return nullptr;
#endif /* HYPRE_OLD_VERSION */
  }    // end of setHypreILUPreconditioner

  std::unique_ptr<LinearSolverPreconditioner> setHypreParaSailsPreconditioner(
      Context&, FiniteElementSpace<true>&, const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    auto ps = std::make_unique<mfem::HypreParaSails>(MPI_COMM_WORLD);
    checkParameters(opts, {Problem::SolverVerbosityLevel});
    // if (contains(opts, Problem::SolverVerbosityLevel)) {
    //  ps->SetPrintLevel(get<int>(opts, Problem::SolverVerbosityLevel));
    //}
    return ps;
  }  // end of setHypreParaSailsPreconditioner

  std::unique_ptr<LinearSolverPreconditioner> setHypreDiagScalePreconditioner(
      Context&, FiniteElementSpace<true>&, const Parameters& opts) {
    using Problem = AbstractNonLinearEvolutionProblem;
    // auto ps = std::make_unique<mfem::HypreDiagScale>(MPI_COMM_WORLD);
    auto ps = std::make_unique<mfem::HypreDiagScale>();
    checkParameters(opts, {Problem::SolverVerbosityLevel});
    return ps;
  }   // end of setHypreDiagScalePreconditioner
#else /* MFEM_USE_MPI */

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreBoomerAMGPreconditioner(Context&,
                                  FiniteElementSpace<true>&,
                                  const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreBoomerAMGPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreEuclidPreconditioner(Context&,
                               FiniteElementSpace<true>&,
                               const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreEuclidPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreILUPreconditioner(Context&,
                            FiniteElementSpace<true>&,
                            const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreILUPreconditioner

  [[noreturn]] std::unique_ptr<LinearSolverPreconditioner>
  setHypreParaSailsPreconditioner(Context&,
                                  FiniteElementSpace<true>&,
                                  const Parameters&) {
    reportUnsupportedParallelComputations();
  }  // end of setHypreParaSailsPreconditioner

#endif /* MFEM_USE_MPI */

  template <bool parallel>
  std::unique_ptr<LinearSolverPreconditioner> getLinearSolverPreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& pr) {
    checkParameters(pr, {"Name", "Options"});
    const auto name = get<std::string>(pr, "Name");
    if (name == "None") {
      if (contains(pr, "Options")) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "no options expected for preconditioner '" +
            name + "'");
      }
      return {};
    } else if (name == "HypreBoomerAMG") {
      if constexpr (parallel) {
        return setHypreBoomerAMGPreconditioner(
            ctx, fespace, get_if<Parameters>(pr, "Options", Parameters{}));
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreBoomerAMG' is only available in parallel");
      }
    } else if (name == "HypreEuclid") {
      if constexpr (parallel) {
        return setHypreEuclidPreconditioner(
            ctx, fespace, get_if<Parameters>(pr, "Options", Parameters{}));
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreEuclid' is only available in parallel");
      }
    } else if (name == "HypreILU") {
      if constexpr (parallel) {
        return setHypreILUPreconditioner(
            ctx, fespace, get_if<Parameters>(pr, "Options", Parameters{}));
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreILU' is only available in parallel");
      }
    } else if (name == "HypreParaSails") {
      if constexpr (parallel) {
        return setHypreParaSailsPreconditioner(
            ctx, fespace, get_if<Parameters>(pr, "Options", Parameters{}));
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreParaSails' is only available in parallel");
      }
    } else if (name == "HypreDiagScale") {
      if constexpr (parallel) {
        return setHypreDiagScalePreconditioner(
            ctx, fespace, get_if<Parameters>(pr, "Options", Parameters{}));
      } else {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "the 'HypreDiagScale' is only available in parallel");
      }
    } else {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: "
          "unsupported preconditioner '" +
          name + "'");
    }
    return {};
  }  // end of getLinearSolverPreconditioner

  template <bool parallel>
  static std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  setLinearSolverParameters(Context& ctx,
                            mfem::IterativeSolver& s,
                            FiniteElementSpace<parallel>& fespace,
                            const Parameters& params) {
    const char* const Preconditioner = "Preconditioner";
    s.iterative_mode = false;
    auto allowed_parameters = getIterativeSolverParametersList();
    allowed_parameters.push_back(Preconditioner);
    checkParameters(params, allowed_parameters);
    setSolverParameters(s, extract(params, getIterativeSolverParametersList()));
    if (contains(params, Preconditioner)) {
      const auto pr = get<Parameters>(params, Preconditioner);
      auto prec = getLinearSolverPreconditioner<parallel>(ctx, fespace, pr);
      if (isInvalid(prec)) {
        return {};
      }
      s.SetPreconditioner(*prec);
      return prec;
    }
    return std::unique_ptr<LinearSolverPreconditioner>{};
  }  // end of setLinearSolverParameters

#ifdef MFEM_USE_MPI

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHyprePCGSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner};
      auto s = std::make_unique<mfem::HyprePCG>(MPI_COMM_WORLD);
      s->iterative_mode = false;
      checkParameters(params, allowed_parameters);
      if (contains(params, Problem::SolverVerbosityLevel)) {
        s->SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
      }
      if (contains(params, SolverTolerance)) {
        s->SetTol(get<double>(params, SolverTolerance));
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        s->SetMaxIter(
            get<int>(params, Problem::SolverMaximumNumberOfIterations));
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto pr = get<Parameters>(params, Preconditioner);
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, pr);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHyprePCGSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HyprePCG linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHyprePCGSolverGenerator

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreGMRESSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner};
      auto s = std::make_unique<mfem::HypreGMRES>(MPI_COMM_WORLD);
      s->iterative_mode = false;
      checkParameters(params, allowed_parameters);
      if (contains(params, Problem::SolverVerbosityLevel)) {
        s->SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        s->SetAbsTol(get<double>(params, Problem::SolverAbsoluteTolerance));
      }
      if (contains(params, SolverTolerance)) {
        s->SetTol(get<double>(params, SolverTolerance));
      }
      if (contains(params, Krylov_Dimension)) {
        s->SetTol(get<double>(params, Krylov_Dimension));
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        s->SetMaxIter(
            get<int>(params, Problem::SolverMaximumNumberOfIterations));
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto pr = get<Parameters>(params, Preconditioner);
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, pr);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHypreGMRESSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HypreGMRES linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHypreGMRESSolverGenerator

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreFGMRESSolverGenerator() {
    return [](Context& ctx, FiniteElementSpace<true>& fespace,
              const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      const char* const SolverTolerance = "Tolerance";
      const char* const Preconditioner = "Preconditioner";
      const char* const Krylov_Dimension = "KDim";
      const auto allowed_parameters = std::vector<std::string>{
          Problem::SolverVerbosityLevel, SolverTolerance,
          Problem::SolverMaximumNumberOfIterations, Preconditioner};
      auto s = std::make_unique<mfem::HypreFGMRES>(MPI_COMM_WORLD);
      s->iterative_mode = false;
      checkParameters(params, allowed_parameters);
      if (contains(params, Problem::SolverVerbosityLevel)) {
        s->SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        s->SetTol(get<double>(params, Problem::SolverAbsoluteTolerance));
      }
      if (contains(params, SolverTolerance)) {
        s->SetTol(get<double>(params, SolverTolerance));
      }
      if (contains(params, Krylov_Dimension)) {
        s->SetTol(get<double>(params, Krylov_Dimension));
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        s->SetMaxIter(
            get<int>(params, Problem::SolverMaximumNumberOfIterations));
      }
      auto prec = std::unique_ptr<mfem::Solver>{};
      if (contains(params, Preconditioner)) {
        const auto pr = get<Parameters>(params, Preconditioner);
        prec = getLinearSolverPreconditioner<true>(ctx, fespace, pr);
        if (isInvalid(prec)) {
          return {};
        }
        auto* const hypre_prec =
            dynamic_cast<mfem::HypreSolver* const>(prec.get());
        if (hypre_prec == nullptr) {
          return ctx.registerErrorMessage(
              "buildHypreFGMRESSolverGenerator: only Hypre preconditioners "
              "allowed with "
              "the HypreFGMRES linear solver");
        }
        s->SetPreconditioner(*hypre_prec);
      }
      return LinearSolverHandler{std::move(s), std::move(prec)};
    };
  }  // end of buildHypreFGMRESSolverGenerator

#else

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHyprePCGSolverGenerator();

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreGMRESSolverGenerator();

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildHypreFGMRESSolverGenerator();

#endif

  template <bool parallel, typename LinearSolverType>
  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<parallel>&, const Parameters&)>
  buildIterativeSolverGenerator() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return [](Context& ctx, FiniteElementSpace<true>& fespace,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>(MPI_COMM_WORLD);
        auto oprec = setLinearSolverParameters<true>(ctx, *s, fespace, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
#else  /* MFEM_USE_MPI */
      return {};
#endif /* MFEM_USE_MPI */
    } else {
      return [](Context& ctx, FiniteElementSpace<false>& fespace,
                const Parameters& params) {
        auto s = std::make_unique<LinearSolverType>();
        auto oprec = setLinearSolverParameters<false>(ctx, *s, fespace, params);
        if (isInvalid(oprec)) {
          return LinearSolverHandler{};
        }
        return LinearSolverHandler{std::move(s), std::move(*oprec)};
      };
    }
  }  // end of buildIterativeSolverGenerator

#ifdef MFEM_USE_MUMPS

  std::function<LinearSolverHandler(
      Context&, FiniteElementSpace<true>&, const Parameters&)>
  buildMUMPSSolverGenerator() {
    return [](Context&, FiniteElementSpace<true>&, const Parameters& params) {
      checkParameters(params, {"Symmetric", "PositiveDefinite"});
      auto s = std::make_unique<mfem::MUMPSSolver>(MPI_COMM_WORLD);
      const auto symmetric = get_if<bool>(params, "Symmetric", false);
      const auto positive_definite =
          get_if<bool>(params, "PositiveDefinite", false);
      s->SetPrintLevel(1);
      //      if (getMPIrank() == 0) {
      //	mfem_mgis::getOutputStream() << "Global Nbdof " <<
      // p.getFiniteElementSpace().GlobalTrueVSize() << "\n";
      //      }
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
      return LinearSolverHandler{std::move(s),
                                 std::unique_ptr<LinearSolverPreconditioner>{}};
    };
  }  // end of builMUMPSGenerator

#endif /* MFEM_USE_MUMPS */

  template <bool parallel>
  static void declareDefaultSolvers(LinearSolverFactory<parallel>& f) {
    auto ctx = Context{};
    auto check = [&ctx](const bool b) {
      if (!b) {
        auto eh = mgis::ContractViolationHandler{};
        eh.registerErrorMessage(ctx.getErrorMessage().c_str());
      }
    };
    check(f.add(ctx, "CGSolver",
                buildIterativeSolverGenerator<parallel, mfem::CGSolver>()));
    check(f.add(ctx, "GMRESSolver",
                buildIterativeSolverGenerator<parallel, mfem::GMRESSolver>()));
    check(
        f.add(ctx, "BiCGSTABSolver",
              buildIterativeSolverGenerator<parallel, mfem::BiCGSTABSolver>()));
    check(f.add(ctx, "MINRESSolver",
                buildIterativeSolverGenerator<parallel, mfem::MINRESSolver>()));
    check(f.add(ctx, "SLISolver",
                buildIterativeSolverGenerator<parallel, mfem::SLISolver>()));
    if constexpr (parallel) {
#ifdef MFEM_USE_MUMPS
      check(f.add(ctx, "MUMPSSolver", buildMUMPSSolverGenerator()));
#endif
      check(f.add(ctx, "HyprePCG", buildHyprePCGSolverGenerator()));
      check(f.add(ctx, "HypreGMRES", buildHypreGMRESSolverGenerator()));
      check(f.add(ctx, "HypreFGMRES", buildHypreFGMRESSolverGenerator()));
    }
    if constexpr (!parallel) {
#ifdef MFEM_USE_SUITESPARSE
      check(f.add(ctx, "UMFPackSolver",
                  [](Context&, FiniteElementSpace<false>&, const Parameters&) {
                    return LinearSolverHandler{
                        std::make_unique<mfem::UMFPackSolver>(),
                        std::unique_ptr<LinearSolverPreconditioner>{}};
                  }));
#endif
    }
  }  // end of declareDefaultSolvers

#ifdef MFEM_USE_MPI

  LinearSolverFactory<true>& LinearSolverFactory<true>::getFactory() {
    static LinearSolverFactory<true> factory;
    return factory;
  }  // end of getFactory

  bool LinearSolverFactory<true>::add(Context& ctx,
                                      std::string_view n,
                                      Generator g) noexcept {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      auto msg = std::string{"LinearSolverFactory<true>::add: "};
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      return ctx.registerErrorMessage(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
    return true;
  }  // end of add

  LinearSolverHandler LinearSolverFactory<true>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<true>& fespace,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<true>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared. Here is the list of available linear solvers:";
      for (const auto& g : this->generators) {
        msg += " " + g.first;
      }
      return ctx.registerErrorMessage(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(ctx, fespace, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "error while generating no linear '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      return ctx.registerErrorMessage(msg);
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

  bool LinearSolverFactory<false>::add(Context& ctx,
                                       std::string_view n,
                                       Generator g) noexcept {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      auto msg = std::string{"LinearSolverFactory<false>::add: "};
      msg += "a linear solver called '";
      msg += n;
      msg += "' has already been declared";
      return ctx.registerErrorMessage(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
    return true;
  }  // end of add

  LinearSolverHandler LinearSolverFactory<false>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<false>& fespace,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "no linear solver called '";
      msg += n;
      msg += "' declared. Here is the list of available linear solvers:";
      for (const auto& g : this->generators) {
        msg += " " + g.first;
      }
      return ctx.registerErrorMessage(msg);
    }
    const auto& g = pg->second;
    LinearSolverHandler s;
    try {
      s = g(ctx, fespace, params);
    } catch (std::exception& e) {
      std::string msg("LinearSolverFactory<false>::generate: ");
      msg += "error while generating no linear '";
      msg += n;
      msg += "'\n";
      msg += e.what();
      return ctx.registerErrorMessage(msg);
    }
    return s;
  }  // end of generate

  LinearSolverFactory<false>::LinearSolverFactory() {
    declareDefaultSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<false>::~LinearSolverFactory() = default;

}  // end of namespace mfem_mgis
