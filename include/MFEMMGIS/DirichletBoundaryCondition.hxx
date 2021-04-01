/*!
 * \file   include/MFEMMGIS/DirichletBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITION_HXX

#include <memory>
#include <vector>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief abstract class for the definition of Dirichlet boundary conditions.
   */
  struct MFEM_MGIS_EXPORT DirichletBoundaryCondition {
    /*!
     * \return the list of degrees of freedom treated by this boundary
     * condition
     */
    virtual std::vector<size_type> getHandledDegreesOfFreedom() const = 0;
    /*!
     * \brief update the values of the imposed degrees of freedom
     * \param[in] u: unknown vector
     * \param[in] t: time at the end of the time step
     */
    virtual void updateImposedValues(mfem::Vector&, const real) const = 0;
    /*!
     * \brief update the values of the imposed degrees of freedom
     * \param[in] u: unknown vector
     * \param[in] ti: time at the beginning of the time step
     * \param[in] te: time at the end of the time step
     */
    virtual void setImposedValuesIncrements(mfem::Vector&,
                                            const real,
                                            const real) const = 0;
    //! \brief destructor
    virtual ~DirichletBoundaryCondition();
  };  // end of struct DirichletBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITION_HXX */
