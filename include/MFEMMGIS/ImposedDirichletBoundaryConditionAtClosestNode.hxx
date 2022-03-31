/*!
 * \file   include/MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX
#define LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX

#include <array>
#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"

namespace mfem_mgis {

  // forwar declaration
  struct AbstractNonLinearEvolutionProblem;
  struct Parameters;

  /*!
   * \brief an helper structure to block the closest point to the given poistion
   * along the a specified component.
   */
  struct MFEM_MGIS_EXPORT ImposedDirichletBoundaryConditionAtClosestNode
      : public DirichletBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization>,
        const std::array<real, 2u>,
        const size_type);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     * \param[in] uvalues: function returning the imposed values
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization>,
        const std::array<real, 2u>,
        const size_type,
        std::function<real(const real)>);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization>,
        const std::array<real, 3u>,
        const size_type);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     * \param[in] uvalues: function returning the imposed values
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization>,
        const std::array<real, 3u>,
        const size_type,
        std::function<real(const real)>);
    //
    std::vector<size_type> getHandledDegreesOfFreedom() const override;
    void updateImposedValues(mfem::Vector&, const real) const override;
    void setImposedValuesIncrements(mfem::Vector&,
                                    const real,
                                    const real) const override;
    //! \brief destructor
    ~ImposedDirichletBoundaryConditionAtClosestNode() override;

   protected:
    //! \brief function returning the value of the imposed displacement
    std::function<real(const real)> ufct;
    //! \brief degree of freedomon blocked
    const std::optional<size_type> dof;
  };  // end of struct ImposedDirichletBoundaryConditionAtClosestNode

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX */
