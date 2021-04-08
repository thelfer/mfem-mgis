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
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;
  // forward declaration
  struct DirichletBoundaryCondition;
  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct Material;
  // forward declaration
  struct BehaviourIntegrator;
  // forward declaration
  struct MultiMaterialNonLinearIntegrator;
  // forward declaration
  struct NewtonSolver;
  // forward declaration
  enum struct IntegrationType;

  /*!
   * \brief class for solving non linear evolution problems.
   *
   * By default, we use of the `MultiMaterialNonLinearIntegrator` class to
   * compute the inner forces contributions to the residual.
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementationBase
      : AbstractNonLinearEvolutionProblem {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief name of the parameter used to activate/desactivate
     * this use of the `MultiMaterialNonLinearIntegrator` class
     */
    static const char* const UseMultiMaterialNonLinearIntegrator;
    //! \return the list of valid parameters
    static std::vector<std::string> getParametersList();
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblemImplementationBase(
        std::shared_ptr<FiniteElementDiscretization>,
        const Hypothesis,
        const Parameters&);
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(const std::vector<real>&);
    /*!
     * \brief set the linear solver
     * \param[in] s: linear solver
     */
    virtual void updateLinearSolver(std::unique_ptr<LinearSolver>);
    //
    FiniteElementDiscretization& getFiniteElementDiscretization() override;
    std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() override;
    mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() override;
    const mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() const override;
    mfem::Vector& getUnknownsAtEndOfTheTimeStep() override;
    const mfem::Vector& getUnknownsAtEndOfTheTimeStep() const override;
    void setSolverParameters(const Parameters&) override;
    std::vector<size_type> getMaterialIdentifiers() const override;
    const Material& getMaterial(const size_type) const override;
    Material& getMaterial(const size_type) override;
    const BehaviourIntegrator& getBehaviourIntegrator(
        const size_type) const override;
    BehaviourIntegrator& getBehaviourIntegrator(const size_type) override;
    void addBehaviourIntegrator(const std::string&,
                                const size_type,
                                const std::string&,
                                const std::string&) override;
    void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    void revert() override;
    void update() override;
    bool solve(const real, const real) override;
    //! \brief destructor
    virtual ~NonLinearEvolutionProblemImplementationBase();

   protected:
    /*!
     * \brief declare the degrees of freedom handled by Dirichlet boundary
     * conditions.
     * \param[in] dofs: list of degrees of freedom
     */
    virtual void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type>) = 0;
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
     * \brief compute prediction
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void computePrediction(const real, const real);
    /*!
     * \brief integrate the behaviour for given estimate of the unknowns at
     * the end of the time step.
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     */
    virtual bool integrate(const mfem::Vector&, const IntegrationType) = 0;
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
    //! \brief newton solver
    std::unique_ptr<NewtonSolver> solver;
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    /*!
     * \brief pointer to the underlying domain integrator
     * The memory associated with this pointer must be released in derived class
     */
    MultiMaterialNonLinearIntegrator* const mgis_integrator;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;

  };  // end of struct NonLinearEvolutionProblemImplementationBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX */
