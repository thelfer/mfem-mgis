/*!
 * \file   include/MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx
 * \brief  This file declares the
 * `ImposedDirichletBoundaryConditionAtClosestNode` class
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX
#define LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX

#include <array>
#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractDirichletBoundaryCondition.hxx"

namespace mfem_mgis {

  // forwar declaration
  struct AbstractNonLinearEvolutionProblem;
  struct Parameters;

  /*!
   * \brief a helper structure to impose the value of a specified component
   * of the unknowns at the node closest to the given position.
   */
  struct MFEM_MGIS_EXPORT ImposedDirichletBoundaryConditionAtClosestNode
      : public AbstractDirichletBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::array<real, 2u> pt,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     * \param[in] uvalues: function returning the imposed values
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::array<real, 2u> pt,
        const size_type c,
        std::function<real(const real)> uvalues);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::array<real, 3u> pt,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] pt: position of the point
     * \param[in] c: component blocked
     * \param[in] uvalues: function returning the imposed values
     */
    ImposedDirichletBoundaryConditionAtClosestNode(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::array<real, 3u> pt,
        const size_type c,
        std::function<real(const real)> uvalues);
    /*!
     * \return the list of degrees of freedom treated by this boundary
     * condition, empty if the closest node is not handled by the current
     * process
     */
    std::vector<size_type> getHandledDegreesOfFreedom() const override;
    /*!
     * \brief update the value of the imposed degree of freedom
     * \param[in, out] u: unknown vector
     * \param[in] t: time at the end of the time step
     */
    void updateImposedValues(mfem::Vector& u, const real t) const override;
    /*!
     * \brief set the increment of the imposed degree of freedom between the
     * two given times, multiplied by the given factor
     * \param[in, out] du: increment of the unknowns
     * \param[in] ti: time at the beginning of the time step
     * \param[in] te: time at the end of the time step
     * \param[in] f: multiplicative factor
     */
    void setImposedValuesIncrements(mfem::Vector& du,
                                    const real ti,
                                    const real te,
                                    const real f) const override;
    //! \brief destructor
    ~ImposedDirichletBoundaryConditionAtClosestNode() override;

   protected:
    //! \brief function returning the value of the imposed displacement
    std::function<real(const real)> ufct;
    //! \brief blocked degree of freedom, if handled by the current process
    const std::optional<size_type> dof;
  };  // end of struct ImposedDirichletBoundaryConditionAtClosestNode

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_IMPOSEDDIRICHLETBOUNDARYCONDITIONATCLOSESTNODE_HXX */
