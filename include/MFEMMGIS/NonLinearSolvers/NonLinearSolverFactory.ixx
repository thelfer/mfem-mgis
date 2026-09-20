/*!
 * \file   MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.ixx
 * \brief
 * \author th202608
 * \date   20/09/2026
 */

#ifndef LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_IXX
#define LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_IXX

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  template <std::derived_from<AbstractNonLinearSolver> SolverType>
  std::unique_ptr<AbstractNonLinearSolver>
  StandardNonLinearSolverGenerator<SolverType>::operator()(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& parameters) noexcept {
    auto ptr = make_unique<SolverType>(ctx, p);
    if (isInvalid(ptr)) {
      return {};
    }
    if (isInvalid(ptr->setSolverParameters(ctx, parameters))) {
      return {};
    }
    return ptr;
  }  // end of operator()

#endif /* MFEM_USE_MPI */

  template <std::derived_from<AbstractNonLinearSolver> SolverType>
  std::unique_ptr<AbstractNonLinearSolver>
  StandardNonLinearSolverGenerator<SolverType>::operator()(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& parameters) noexcept {
    auto ptr = make_unique<SolverType>(ctx, p);
    if (isInvalid(ptr)) {
      return {};
    }
    if (isInvalid(ptr->setSolverParameters(ctx, parameters))) {
      return {};
    }
    return ptr;
  }  // end of operator()

  template <std::derived_from<AbstractNonLinearSolver> SolverType>
  StandardNonLinearSolverGenerator<
      SolverType>::~StandardNonLinearSolverGenerator() noexcept = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_IXX */
