/*!
 * \file   src/AnalyticalTests.cxx
 * \brief
 * \author Thomas Helfer
 * \date   25/03/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static real getL2Error(
      mfem_mgis::NonLinearEvolutionProblem &p,
      std::function<void(const mfem::Vector &, mfem::Vector &)> f) {
    auto &problem = p.template getImplementation<parallel>();
    mfem_mgis::GridFunction<parallel> x(&problem.getFiniteElementSpace());
    const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
    // recover the solution as a grid function
    auto &u1 = problem.getUnknownsAtEndOfTheTimeStep();
    x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
    x.SetFromTrueVector();
    // comparison to analytical solution
    mfem::VectorFunctionCoefficient sol_coef(dim, f);
    return x.ComputeL2Error(sol_coef);
  }

  real computeL2ErrorAgainstAnalyticalSolution(
      NonLinearEvolutionProblem &p,
      std::function<void(mfem::Vector &, const mfem::Vector &)> f) {
    std::function<void(const mfem::Vector &, mfem::Vector &)> mf =
        [f](const mfem::Vector &x, mfem::Vector &v) { f(v, x); };
    const auto &fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return getL2Error<true>(p, mf);
#else
      raise(
          "computeL2ErrorAgainstAnalyticalSolution: "
          "unsupported parallel computations");
#endif
    }
    return getL2Error<false>(p, mf);
  }  //  end of computeL2ErrorAgainstAnalyticalSolution

  bool compareToAnalyticalSolution(
      NonLinearEvolutionProblem &p,
      std::function<void(mfem::Vector &, const mfem::Vector &)> f,
      const Parameters &params) {
    const auto error = computeL2ErrorAgainstAnalyticalSolution(p, f);
    const auto verbosity = get_if<int>(params, "VerbosityLevel", 0);
    if (verbosity > 1) {
      mfem_mgis::getOutputStream() << "L2Error: " << error << "\n";
    }
    return error < get<real>(params, "CriterionThreshold");
  }  // end of compareToAnalyticalSolution

}  // end of namespace mfem_mgis
