/*!
 * \file   include/MFEMMGIS/AbstractDirichletBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTDIRICHLETBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_ABSTRACTDIRICHLETBOUNDARYCONDITION_HXX

#include <memory>
#include <vector>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief abstract class for the definition of Dirichlet boundary conditions.
   */
  struct MFEM_MGIS_EXPORT AbstractDirichletBoundaryCondition {
    /*!
     * \return the list of degrees of freedom treated by this boundary
     * condition
     */
    virtual std::vector<size_type> getHandledDegreesOfFreedom() const = 0;
    /*!
     * \brief update the values of the imposed degrees of freedom
     * \param[in, out] u: unknown vector
     * \param[in] t: time at the end of the time step
     */
    virtual void updateImposedValues(mfem::Vector& u, const real t) const = 0;
    /*!
     * \brief set the increments of the imposed degrees of freedom between
     * the two given times, multiplied by the given factor
     * \param[in, out] du: increment of the unknowns
     * \param[in] ti: time at the beginning of the time step
     * \param[in] te: time at the end of the time step
     * \param[in] f: multiplicative factor
     */
    virtual void setImposedValuesIncrements(mfem::Vector& du,
                                            const real ti,
                                            const real te,
                                            const real f) const = 0;
    //! \brief destructor
    virtual ~AbstractDirichletBoundaryCondition();
  };  // end of struct AbstractDirichletBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTDIRICHLETBOUNDARYCONDITION_HXX */
