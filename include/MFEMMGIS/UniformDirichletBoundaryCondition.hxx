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
     * \param[in] params: parameters defining the boundary condition
     */
    UniformDirichletBoundaryCondition(AbstractNonLinearEvolutionProblem& p,
                                      const Parameters& params);
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
        std::shared_ptr<FiniteElementDiscretization> fed,
        const size_type bid,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
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
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
     * condition.
     *
     * \note the degree of freedom are set to zero
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view bid,
        const size_type c);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] c: component of the unknows treated by this boundary
     * condition.
     * \param[in] uvalues: function returning the imposed values
     */
    UniformDirichletBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view bid,
        const size_type c,
        std::function<real(const real)> uvalues);
    //
    void updateImposedValues(mfem::Vector& u, const real t) const override;
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
