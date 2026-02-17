
/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX

#include <memory>
#include <vector>
#include "mfem/linalg/vector.hpp"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

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
  struct AbstractBehaviourIntegrator;
  // forward declaration
  struct MultiMaterialNonLinearIntegrator;
  // forward declaration
  struct NewtonSolver;
  // forward declaration
  enum struct IntegrationType;
  // forward declaration
  struct LinearSolverHandler;

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
    /*!
     * \brief set the linear solver
     * \param[in] s: linear solver
     * \param[in] p: linear solver preconditioner
     */
    virtual void updateLinearSolver(
        std::unique_ptr<LinearSolver>,
        std::unique_ptr<LinearSolverPreconditioner>);
    /*!
     * \brief set the linear solver
     * \param[in] s: linear solver handler
     */
    virtual void updateLinearSolver(LinearSolverHandler);
    //
    FiniteElementDiscretization& getFiniteElementDiscretization() override;
    const FiniteElementDiscretization& getFiniteElementDiscretization()
        const override;
    std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() override;
    mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() override;
    const mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() const override;
    mfem::Vector& getUnknownsAtEndOfTheTimeStep() override;
    const mfem::Vector& getUnknownsAtEndOfTheTimeStep() const override;
    [[nodiscard]] mfem::Vector& getUnknowns(
        const TimeStepStage) noexcept override;
    [[nodiscard]] const mfem::Vector& getUnknowns(
        const TimeStepStage) const noexcept override;
    void setSolverParameters(const Parameters&) override;
    void setMaterialsNames(const std::map<size_type, std::string>&) override;
    void setBoundariesNames(const std::map<size_type, std::string>&) override;
    std::vector<size_type> getAssignedMaterialsIdentifiers() const override;
    size_type getMaterialIdentifier(const Parameter&) const override;
    size_type getBoundaryIdentifier(const Parameter&) const override;
    std::vector<size_type> getMaterialsIdentifiers(
        const Parameter&) const override;
    std::vector<size_type> getBoundariesIdentifiers(
        const Parameter&) const override;
    const Material& getMaterial(const Parameter&) const override;
    Material& getMaterial(const Parameter&) override;
    const AbstractBehaviourIntegrator& getBehaviourIntegrator(
        const size_type) const override;
    AbstractBehaviourIntegrator& getBehaviourIntegrator(
        const size_type) override;
    void addBehaviourIntegrator(const std::string&,
                                const Parameter&,
                                const std::string&,
                                const std::string&) override;
    [[nodiscard]] std::vector<size_type> getEssentialDegreesOfFreedom()
        const override;
    [[nodiscard]] std::optional<LinearizedOperators> getLinearizedOperators(
        Context&, const mfem::Vector&) noexcept override;
    [[nodiscard]] const std::vector<
        std::unique_ptr<DirichletBoundaryCondition>>&
    getDirichletBoundaryConditions() noexcept override;
    void setup(const real, const real) override;
    void setPredictionPolicy(const PredictionPolicy &) noexcept override;
    [[nodiscard]] NonLinearResolutionOutput solve(const real,
                                                  const real) override;
    void revert() override;
    void update() override;
    //! \brief destructor
    virtual ~NonLinearEvolutionProblemImplementationBase();

   protected:
    /*!
     * \brief declare the degrees of freedom handled by Dirichlet boundary
     * conditions.
     * \param[in] dofs: list of degrees of freedom
     *
     * \note copy is required to create a mutable mfem::Array
     */
    virtual void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type>) = 0;
    /*!
     * \brief set the time time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real);
    /*!
     * \brief compute prediction
     * \param[in] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<real> computePrediction(
        Context&, const real, const real) noexcept = 0;
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
#ifdef MFEM_USE_PETSC
    //! \brief newton solver
    std::unique_ptr<mfem::PetscNonlinearSolver> petsc_solver;
#endif /* MFEM_USE_PETSC */
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    //! \brief linear solver preconditioner
    std::unique_ptr<LinearSolverPreconditioner> linear_solver_preconditioner;
    /*!
     * \brief pointer to the underlying domain integrator
     * The memory associated with this pointer must be released in derived class
     */
    MultiMaterialNonLinearIntegrator* const mgis_integrator = nullptr;
    //! \brief registred boundary conditions
    std::vector<std::unique_ptr<AbstractBoundaryCondition>> boundary_conditions;
    //! \brief prediction policy
    PredictionPolicy prediction_policy;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
  };  // end of struct NonLinearEvolutionProblemImplementationBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX */
