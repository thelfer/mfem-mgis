/*!
 * \file   MFEMMGIS/AbstractBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for non essential boundary conditions
   */
  struct MFEM_MGIS_EXPORT AbstractBoundaryCondition {
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the nonlinear form integrator describing the
     * the boundary condition
     *
     * \param[in] ctx: execution context
     * \param[in] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     */
    [[nodiscard]] virtual bool addNonlinearFormIntegrator(
        Context&, NonlinearForm<true>&, const mfem::Vector&) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the nonlinear form integrator describing the
     * the boundary condition
     *
     * \param[in] ctx: execution context
     * \param[in] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     */
    [[nodiscard]] virtual bool addNonlinearFormIntegrator(
        Context&, NonlinearForm<false>&, const mfem::Vector&) noexcept = 0;
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the bilinear and linear form integrators
     * describing the boundary condition
     *
     * \param[in] ctx: execution context
     * \param[in] a: form computing the jacobian matrix
     * \param[in] b: form computing the right hand side
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual bool addLinearFormIntegrators(
        Context&,
        BilinearForm<true>&,
        LinearForm<true>&,
        const mfem::Vector&,
        const real,
        const real) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the linear form integrator describing the
     * the boundary condition
     *
     * \param[in] ctx: execution context
     * \param[in] a: form computing the jacobian matrix
     * \param[in] b: form computing the right hand side
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual bool addLinearFormIntegrators(
        Context&,
        BilinearForm<false>&,
        LinearForm<false>&,
        const mfem::Vector&,
        const real,
        const real) noexcept = 0;
    /*!
     * \brief method call at the beginning of each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real) = 0;
    //! \brief destructor
    virtual ~AbstractBoundaryCondition();
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTBOUNDARYCONDITION_HXX */
