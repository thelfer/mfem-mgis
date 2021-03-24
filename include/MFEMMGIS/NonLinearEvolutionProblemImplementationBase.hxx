/*!
 * \file   NonLinearEvolutionProblemImplementationBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX

#include <memory>
#include <vector>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // struct forward declaration
  struct DirichletBoundaryCondition;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementationBase
      : AbstractNonLinearEvolutionProblem {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    NonLinearEvolutionProblemImplementationBase(
        std::shared_ptr<FiniteElementDiscretization>);
    //! \return the underlying finite element discretization
    FiniteElementDiscretization& getFiniteElementDiscretization() override;
    //! \return the underlying finite element discretization
    std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() override;
    //! \return the unknowns at the beginning of the time step
    mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() override;
    //! \return the unknowns at the beginning of the time step
    const mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() const override;
    //! \return the unknowns at the end of the time step
    mfem::Vector& getUnknownsAtEndOfTheTimeStep() override;
    //! \return the unknowns at the end of the time step
    const mfem::Vector& getUnknownsAtEndOfTheTimeStep() const override;
    /*!
     * \brief add a Dirichlet boundary condition
     * \param[in] bc: boundary condition
     */
    void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    //! \brief revert the state to the beginning of the time step.
    void revert() override;
    //! \brief update the state to the end of the time step.
    void update() override;

    //! \brief destructor
    virtual ~NonLinearEvolutionProblemImplementationBase();

   protected:
    /*!
     * \brief set the time time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real);
    /*!
     * \brief method called before each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real);
    /*!
     * \brief declare the degrees of freedom handled by Dirichlet boundary
     * conditions.
     * \param[in] dofs: list of degrees of freedom
     */
    virtual void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type>) = 0;
    //! \brief underlying finite element discretization
    const std::shared_ptr<FiniteElementDiscretization> fe_discretization;
    //! \brief list of boundary conditions
    std::vector<std::unique_ptr<DirichletBoundaryCondition>>
        dirichlet_boundary_conditions;
    /*!
     * \brief a boolean value to specifiy if the initialization phase is still
     * open.
     *
     * This initialization phase ends at the first call to the `setup` method,
     * i.e. at the first call of the `solve` method.
     */
    bool initialization_phase = true;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
  };  // end of struct NonLinearEvolutionProblemImplementationBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX */
