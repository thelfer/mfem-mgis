/*!
 * \file   include/MFEMMGIS/UniformDirichletBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/DirichletBoundaryConditionBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief class used to simplify the definition of Dirichlet boundary
   * conditions.
   */
  struct MFEM_MGIS_EXPORT UniformDirichletBoundaryCondition
      : DirichletBoundaryConditionBase {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
     * condition.
     *
     * \note the degree of freedom are set to zero
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization>,
        const size_type,
        const size_type);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
     * condition.
     * \param[in] uvalues: function returning the imposed values
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization>,
        const size_type,
        const size_type,
        std::function<real(const real)>);
    //
    void updateImposedValues(mfem::Vector&, const real) const override;
    //! \brief destructor
    ~UniformDirichletBoundaryCondition() override;

   protected:
    //! \brief function returning
    std::function<real(const real)> ufct;
  };  // end of struct UniformDirichletBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX */
