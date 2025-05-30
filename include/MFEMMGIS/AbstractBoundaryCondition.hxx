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
     * \param[in] f: form
     */
    virtual void addNonlinearFormIntegrator(NonlinearForm<true>&) = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the nonlinear form integrator describing the
     * the boundary condition
     * \param[in] f: form
     */
    virtual void addNonlinearFormIntegrator(NonlinearForm<false>&) = 0;
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
