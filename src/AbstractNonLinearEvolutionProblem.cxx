/*!
 * \file   src/AbstractNonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
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

  size_type getMaterialIdentifier(const AbstractNonLinearEvolutionProblem &p,
                                  const Parameters &params) {
    return p.getMaterialIdentifier(get(params, "Material"));
  }

  size_type getBoundaryIdentifier(const AbstractNonLinearEvolutionProblem &p,
                                  const Parameters &params) {
    return p.getBoundaryIdentifier(get(params, "Boundary"));
  }

  std::vector<size_type> getMaterialsIdentifiers(
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) {
    if (contains(params, "Material") && contains(params, "Materials")) {
      raise(
          "getMaterialsIdentifiers: "
          "both `Material` and `Materials` parameters specified");
    }
    if (contains(params, "Material")) {
      if (is<std::vector<Parameter>>(params, "Material")) {
        raise("getMaterialsIdentifiers: invalid `Material` parameter");
      }
      return p.getMaterialsIdentifiers(get(params, "Material"));
    } else if (contains(params, "Materials")) {
      return p.getMaterialsIdentifiers(get(params, "Materials"));
    }
    if (!b) {
      raise(
          "getMaterialsIdentifiers: no parameter named `Material` nor "
          "`Materials` given");
    }
    return p.getAssignedMaterialsIdentifiers();
  }  // end of getMaterialsIdentifiers

  std::vector<size_type> getBoundariesIdentifiers(
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b) {
    if (contains(params, "Boundary") && contains(params, "Boundaries")) {
      raise(
          "getBoundariesIdentifiers: "
          "both `Boundary` and `Boundaries` parameters specified");
    }
    if (contains(params, "Boundary")) {
      if (is<std::vector<Parameter>>(params, "Boundary")) {
        raise("getBoundariesIdentifiers: invalid `Boundary` parameter");
      }
      return p.getBoundariesIdentifiers(get(params, "Boundary"));
    } else if (contains(params, "Boundaries")) {
      return p.getBoundariesIdentifiers(get(params, "Boundaries"));
    }
    if (!b) {
      raise(
          "getBoundariesIdentifiers: no parameter named `Boundary` nor "
          "`Boundaries` given");
    }
    return p.getBoundariesIdentifiers(".+");
  }  // end of getBoundariesIdentifiers

  AbstractNonLinearEvolutionProblem::~AbstractNonLinearEvolutionProblem() =
      default;

}  // end of namespace mfem_mgis
