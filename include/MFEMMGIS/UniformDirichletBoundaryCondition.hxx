/*!
 * \file   include/MFEMMGIS/UniformDirichletBoundaryCondition.hxx
 * \brief  This file declares the `UniformDirichletBoundaryCondition` class
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#ifndef LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/DirichletBoundaryConditionBase.hxx"

namespace mfem_mgis {

  // forwar declaration
  struct AbstractNonLinearEvolutionProblem;
  struct Parameters;

  /*!
   * \brief class used to simplify the definition of Dirichlet boundary
   * conditions.
   */
  struct MFEM_MGIS_EXPORT UniformDirichletBoundaryCondition
      : DirichletBoundaryConditionBase {
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     * \param[in] params: parameters defining the boundary condition:
     * `Boundary` or `Boundaries`, `Component` and, optionally,
     * `LoadingEvolution`
     */
    UniformDirichletBoundaryCondition(AbstractNonLinearEvolutionProblem& p,
                                      const Parameters& params);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknowns treated by this boundary
     * condition.
     *
     * \note the degrees of freedom are set to zero
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const size_type bid,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknowns treated by this boundary
     * condition.
     * \param[in] uvalues: function returning the imposed values
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const size_type bid,
        const size_type c,
        std::function<real(const real)> uvalues);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: regular expression selecting the boundaries by name
     * \param[in] c: component of the unknowns treated by this boundary
     * condition.
     *
     * \note the degrees of freedom are set to zero
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view bid,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: regular expression selecting the boundaries by name
     * \param[in] c: component of the unknowns treated by this boundary
     * condition.
     * \param[in] uvalues: function returning the imposed values
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view bid,
        const size_type c,
        std::function<real(const real)> uvalues);
    /*!
     * \brief update the values of the imposed degrees of freedom
     * \param[in, out] u: unknown vector
     * \param[in] t: time at the end of the time step
     */
    void updateImposedValues(mfem::Vector& u, const real t) const override;
    /*!
     * \brief set the increments of the imposed degrees of freedom between
     * the two given times, multiplied by the given factor
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
    ~UniformDirichletBoundaryCondition() override;

   protected:
    //! \brief function returning the value of the imposed displacement
    std::function<real(const real)> ufct;
  };  // end of struct UniformDirichletBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UNIFORMDIRICHLETBOUNDARYCONDITION_HXX */
