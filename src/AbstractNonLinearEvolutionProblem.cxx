/*!
 * \file   src/AbstractNonLinearEvolutionProblem.cxx
 * \brief  This file implements the `AbstractNonLinearEvolutionProblem` class
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  const char *const AbstractNonLinearEvolutionProblem::SolverVerbosityLevel =
      "VerbosityLevel";
  const char *const AbstractNonLinearEvolutionProblem::SolverRelativeTolerance =
      "RelativeTolerance";
  const char *const AbstractNonLinearEvolutionProblem::SolverAbsoluteTolerance =
      "AbsoluteTolerance";
  const char *const
      AbstractNonLinearEvolutionProblem::SolverMaximumNumberOfIterations =
          "MaximumNumberOfIterations";

  AbstractNonLinearEvolutionProblem::~AbstractNonLinearEvolutionProblem() =
      default;

  std::optional<size_type> getMaterialIdentifier(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params) noexcept {
    const auto om = get(ctx, params, "Material");
    if (isInvalid(om)) {
      return {};
    }
    return p.getMaterialIdentifier(ctx, *om);
  }

  std::optional<size_type> getBoundaryIdentifier(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params) noexcept {
    const auto ob = get(ctx, params, "Boundary");
    if (isInvalid(ob)) {
      return {};
    }
    return p.getBoundaryIdentifier(ctx, *ob);
  }

  std::optional<std::vector<size_type>> getMaterialsIdentifiers(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) noexcept {
    if (contains(params, "Material") && contains(params, "Materials")) {
      return ctx.registerErrorMessage(
          "getMaterialsIdentifiers: "
          "both `Material` and `Materials` parameters specified");
    }
    if (contains(params, "Material")) {
      if (is<std::vector<Parameter>>(throwing, params, "Material")) {
        return ctx.registerErrorMessage(
            "getMaterialsIdentifiers: invalid `Material` parameter");
      }
      return p.getMaterialsIdentifiers(ctx, get(throwing, params, "Material"));
    } else if (contains(params, "Materials")) {
      return p.getMaterialsIdentifiers(ctx, get(throwing, params, "Materials"));
    }
    if (!b) {
      return ctx.registerErrorMessage(
          "getMaterialsIdentifiers: no parameter named `Material` nor "
          "`Materials` given");
    }
    return p.getAssignedMaterialsIdentifiers();
  }  // end of getMaterialsIdentifiers

  std::optional<std::vector<size_type>> getBoundariesIdentifiers(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) noexcept {
    if (contains(params, "Boundary") && contains(params, "Boundaries")) {
      return ctx.registerErrorMessage(
          "getBoundariesIdentifiers: "
          "both `Boundary` and `Boundaries` parameters specified");
    }
    if (contains(params, "Boundary")) {
      if (is<std::vector<Parameter>>(throwing, params, "Boundary")) {
        return ctx.registerErrorMessage(
            "getBoundariesIdentifiers: invalid `Boundary` parameter");
      }
      return p.getBoundariesIdentifiers(ctx, get(throwing, params, "Boundary"));
    } else if (contains(params, "Boundaries")) {
      return p.getBoundariesIdentifiers(ctx,
                                        get(throwing, params, "Boundaries"));
    }
    if (!b) {
      return ctx.registerErrorMessage(
          "getBoundariesIdentifiers: no parameter named `Boundary` nor "
          "`Boundaries` given");
    }
    return p.getBoundariesIdentifiers(ctx, ".+");
  }  // end of getBoundariesIdentifiers

  std::optional<BehaviourIntegratorsSelection> getBehaviourIntegratorsSelection(
      Context &ctx, const Parameters &params) noexcept {
    if (!contains(params, "BehaviourIntegrator")) {
      return BehaviourIntegratorsSelection{};
    }
    const auto ob = get(ctx, params, "BehaviourIntegrator");
    if (isInvalid(ob)) {
      return {};
    }
    if (is<int>(*ob)) {
      const auto oi = get<int>(ctx, *ob);
      if (isInvalid(oi)) {
        return {};
      }
      if (*oi < 0) {
        return ctx.registerErrorMessage(
            "getBehaviourIntegratorsSelection: negative index of behaviour "
            "integrator");
      }
      return BehaviourIntegratorsSelection{.index = *oi};
    }
    if (is<std::string>(*ob)) {
      const auto os = get<std::string>(ctx, *ob);
      if (isValid(os) && (*os == "All")) {
        return BehaviourIntegratorsSelection{.all = true};
      }
    }
    return ctx.registerErrorMessage(
        "getBehaviourIntegratorsSelection: the `BehaviourIntegrator` parameter "
        "must be an integer or the string `All`");
  }  // end of getBehaviourIntegratorsSelection

  std::optional<std::pair<size_type, size_type>>
  getSelectedBehaviourIntegrators(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const size_type m,
      const BehaviourIntegratorsSelection &s) noexcept {
    const auto on = p.getNumberOfBehaviourIntegrators(ctx, m);
    if (isInvalid(on)) {
      return {};
    }
    if (*on == 0) {
      return ctx.registerErrorMessage(
          "getSelectedBehaviourIntegrators: no behaviour integrator defined "
          "for material '" +
          std::to_string(m) + "'");
    }
    if (s.all) {
      return std::pair<size_type, size_type>{0, *on};
    }
    if (s.index.has_value()) {
      if (*(s.index) >= *on) {
        return ctx.registerErrorMessage(
            "getSelectedBehaviourIntegrators: no behaviour integrator of index "
            "'" +
            std::to_string(*(s.index)) + "' for material '" +
            std::to_string(m) + "'");
      }
      return std::pair<size_type, size_type>{*(s.index), *(s.index) + 1};
    }
    if (*on != 1) {
      return ctx.registerErrorMessage(
          "getSelectedBehaviourIntegrators: material '" + std::to_string(m) +
          "' has several behaviour integrators, the `BehaviourIntegrator` "
          "parameter must be given");
    }
    return std::pair<size_type, size_type>{0, 1};
  }  // end of getSelectedBehaviourIntegrators

  size_type getMaterialIdentifier(attributes::Throwing,
                                  const AbstractNonLinearEvolutionProblem &p,
                                  const Parameters &params) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return getMaterialIdentifier(ctx, p, params) | or_raise;
  }

  size_type getBoundaryIdentifier(attributes::Throwing,
                                  const AbstractNonLinearEvolutionProblem &p,
                                  const Parameters &params) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return getBoundaryIdentifier(ctx, p, params) | or_raise;
  }

  std::vector<size_type> getMaterialsIdentifiers(
      attributes::Throwing,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return getMaterialsIdentifiers(ctx, p, params, b) | or_raise;
  }  // end of getMaterialsIdentifiers

  std::vector<size_type> getBoundariesIdentifiers(
      attributes::Throwing,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return getBoundariesIdentifiers(ctx, p, params, b) | or_raise;
  }  // end of getBoundariesIdentifiers

  BehaviourIntegratorsSelection getBehaviourIntegratorsSelection(
      attributes::Throwing, const Parameters &params) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return getBehaviourIntegratorsSelection(ctx, params) | or_raise;
  }  // end of getBehaviourIntegratorsSelection

#ifdef MFEM_USE_MPI

  MPI_Comm getMPICommunicator(
      const AbstractNonLinearEvolutionProblem &p) noexcept {
    return getMPICommunicator(p.getFiniteElementDiscretization());
  }  // end of getMPICommunicator

  bool isMainProcess(const AbstractNonLinearEvolutionProblem &p) noexcept {
    return isMainProcess(p.getFiniteElementDiscretization());
  }  // end of isMainProcess

#endif /* MFEM_USE_MPI */

}  // end of namespace mfem_mgis
