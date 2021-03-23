/*!
 * \file   include/MFEMMGIS/DirichletBoundaryConditionBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITIONBASE_HXX
#define LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITIONBASE_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  //! \brief base class for defining Dirichlet boundary conditions.
  struct MFEM_MGIS_EXPORT DirichletBoundaryConditionBase
      : DirichletBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
     * condition.
     */
    DirichletBoundaryConditionBase(std::shared_ptr<FiniteElementDiscretization>,
                                   const size_type,
                                   const size_type);
    //
    std::vector<size_type> getHandledDegreesOfFreedom() const override;
    //! \brief destructor
    ~DirichletBoundaryConditionBase() override;

   protected:
    //! \brief list of degrees of freedom handled by the boundary conditions
    std::vector<size_type> dofs;
  };  // end of struct DirichletBoundaryConditionBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_DIRICHLETBOUNDARYCONDITIONBASE_HXX */
