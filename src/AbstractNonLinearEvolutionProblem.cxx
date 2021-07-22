/*!
 * \file   src/AbstractNonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  const char *const AbstractNonLinearEvolutionProblem::SolverVerbosityLevel =
      "VerbosityLevel";
  const char *const AbstractNonLinearEvolutionProblem::SolverRelativeTolerance =
      "RelativeTolerance";
  const char *const AbstractNonLinearEvolutionProblem::SolverAbsoluteTolerance =
      "AbsoluteTolerance";
  const char *const AbstractNonLinearEvolutionProblem::SolverMaximumNumberOfIterations =
      "MaximumNumberOfIterations";

  AbstractNonLinearEvolutionProblem::~AbstractNonLinearEvolutionProblem() =
      default;

}  // end of namespace mfem_mgis
