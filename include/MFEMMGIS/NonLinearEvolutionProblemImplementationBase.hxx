/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx
 * \brief  This file declares the `NonLinearEvolutionProblemImplementationBase`
 * class
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

  // forward declarations
  struct Parameters;
  struct DirichletBoundaryCondition;
  struct FiniteElementDiscretization;
  struct Material;
  struct AbstractBehaviourIntegrator;
  struct MultiMaterialNonLinearIntegrator;
  struct NewtonSolver;
  enum struct IntegrationType;
  struct LinearSolverHandler;
  struct AbstractNonLinearSolver;

  /*!
   * \brief class for solving non linear evolution problems.
   *
   * By default, we use the `MultiMaterialNonLinearIntegrator` class to
   * compute the inner forces contributions to the residual.
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementationBase
      : AbstractNonLinearEvolutionProblem {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    //! \brief name of the parameter used to select a nonlinear solver
    static const char* const NonLinearSolver;
    /*!
     * \brief name of the parameter used to activate or deactivate
     * the use of the `MultiMaterialNonLinearIntegrator` class
     */
    static const char* const UseMultiMaterialNonLinearIntegrator;
    //! \return the list of valid parameters
    static std::vector<std::string> getParametersList();
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblemImplementationBase(
        Context& ctx,
        std::shared_ptr<FiniteElementDiscretization> fed,
        const Hypothesis h,
        const Parameters& p);
    /*!
     * \brief set the macroscopic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(const std::vector<real>& g);
    /*!
     * \brief set the linear solver
     * \param[in, out] ctx: execution context
     * \param[in] s: linear solver
     * \return true on success
     */
    [[nodiscard]] virtual bool updateLinearSolver(
        Context& ctx, std::unique_ptr<LinearSolver> s) noexcept;
    /*!
     * \brief set the linear solver
     * \param[in, out] ctx: execution context
     * \param[in] s: linear solver
     * \param[in] p: linear solver preconditioner
     * \return true on success
     */
    [[nodiscard]] virtual bool updateLinearSolver(
        Context& ctx,
        std::unique_ptr<LinearSolver> s,
        std::unique_ptr<LinearSolverPreconditioner> p) noexcept;
    /*!
     * \brief set the linear solver
     * \param[in, out] ctx: execution context
     * \param[in] s: linear solver handler
     * \return true on success
     */
    [[nodiscard]] virtual bool updateLinearSolver(
        Context& ctx, LinearSolverHandler s) noexcept;
    /*!
     * \brief method called before each resolution
     *
     * This method must be called before `solve`: this is not done automatically
     * as this is the case for `NonLinearEvolutionProblem::solve`. This is
     * mostly motivated by unit testing, see `PredictionTest` for an example.
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool setup(Context& ctx,
                                     const real t,
                                     const real dt) noexcept;
    //
    [[nodiscard]] FiniteElementDiscretization&
    getFiniteElementDiscretization() noexcept override;
    [[nodiscard]] const FiniteElementDiscretization&
    getFiniteElementDiscretization() const noexcept override;
    [[nodiscard]] std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() noexcept override;
    [[nodiscard]] bool setMaterialsNames(
        Context& ctx,
        const std::map<size_type, std::string>& ids) noexcept override;
    [[nodiscard]] bool setBoundariesNames(
        Context& ctx,
        const std::map<size_type, std::string>& ids) noexcept override;
    [[nodiscard]] mfem::Vector& getUnknowns(
        const TimeStepStage ts) noexcept override;
    [[nodiscard]] const mfem::Vector& getUnknowns(
        const TimeStepStage ts) const noexcept override;
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined, empty if the multi material support is
     * disabled.
     */
    [[nodiscard]] std::vector<size_type> getAssignedMaterialsIdentifiers()
        const noexcept override;
    [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
        Context& ctx, const Parameter& m) const noexcept override;
    [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
        Context& ctx, const Parameter& m) const noexcept override;
    [[nodiscard]] std::optional<std::vector<size_type>> getMaterialsIdentifiers(
        Context& ctx, const Parameter& m) const noexcept override;
    [[nodiscard]] std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context& ctx,
                             const Parameter& m) const noexcept override;
    OptionalReference<const Material> getMaterial(
        Context& ctx,
        const Parameter& m,
        const size_type b) const noexcept override;
    OptionalReference<Material> getMaterial(
        Context& ctx, const Parameter& m, const size_type b) noexcept override;
    [[nodiscard]] std::optional<size_type> getNumberOfBehaviourIntegrators(
        Context& ctx, const Parameter& m) const noexcept override;
    OptionalReference<const AbstractBehaviourIntegrator> getBehaviourIntegrator(
        Context& ctx,
        const Parameter& m,
        const size_type b) const noexcept override;
    OptionalReference<AbstractBehaviourIntegrator> getBehaviourIntegrator(
        Context& ctx, const Parameter& m, const size_type b) noexcept override;
    std::optional<std::map<size_type, size_type>> addBehaviourIntegrator(
        Context& ctx,
        const std::string& n,
        const Parameter& m,
        const std::string& l,
        const std::string& b) noexcept override;
    std::optional<std::map<size_type, size_type>> addBehaviourIntegrator(
        Context& ctx,
        const std::string& n,
        const Parameter& m,
        const std::string& l,
        const std::string& b,
        const Parameters& params) noexcept override;
    [[nodiscard]] std::vector<size_type> getEssentialDegreesOfFreedom()
        const override;
    [[nodiscard]] bool areStiffnessOperatorsFromLastIterationAvailable()
        const noexcept override;
    [[nodiscard]] std::optional<LinearizedOperators> getLinearizedOperators(
        Context& ctx, const mfem::Vector& U) noexcept override;
    [[nodiscard]] const std::vector<
        std::unique_ptr<AbstractDirichletBoundaryCondition>>&
    getDirichletBoundaryConditions() const noexcept override;
    [[nodiscard]] const std::vector<std::unique_ptr<AbstractBoundaryCondition>>&
    getBoundaryConditions() const noexcept override;
    void setPredictionPolicy(const PredictionPolicy& p) noexcept override;
    [[nodiscard]] bool setSolverParameters(
        Context& ctx, const Parameters& params) noexcept override;
    [[nodiscard]] PredictionPolicy getPredictionPolicy()
        const noexcept override;
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return the output of the non linear resolution
     * \note `setup` must be called before
     */
    [[nodiscard]] NonLinearResolutionOutput solve(
        Context& ctx, const real t, const real dt) noexcept override;
    [[nodiscard]] bool revert(Context& ctx) noexcept override;
    [[nodiscard]] bool update(Context& ctx) noexcept override;
    //! \brief destructor
    ~NonLinearEvolutionProblemImplementationBase() override;

   protected:
    /*!
     * \brief declare the degrees of freedom handled by Dirichlet boundary
     * conditions.
     * \param[in] dofs: list of degrees of freedom
     *
     * \note copy is required to create a mutable mfem::Array
     */
    virtual void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type> dofs) = 0;
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real dt);
    /*!
     * \brief compute a prediction of the unknowns at the end of the time step
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return the norm of the initial residual
     */
    [[nodiscard]] virtual std::optional<real> computePrediction(
        Context& ctx, const real t, const real dt) noexcept = 0;
    //! \brief underlying finite element discretization
    const std::shared_ptr<FiniteElementDiscretization> fe_discretization;
    //! \brief list of Dirichlet boundary conditions
    std::vector<std::unique_ptr<AbstractDirichletBoundaryCondition>>
        dirichlet_boundary_conditions;
    /*!
     * \brief a boolean value to specify if the initialization phase is still
     * open.
     *
     * This initialization phase ends at the first call to the `setup` method,
     * i.e. at the first call of the `NonLinearEvolutionProblem::solve` method.
     */
    bool initialization_phase = true;
    //! \brief if the stiffness operators of the last iteration are available
    bool hasStiffnessOperatorsBeenComputed = false;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
    //! \brief nonlinear solver
    std::unique_ptr<AbstractNonLinearSolver> solver;
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
    //! \brief registered boundary conditions
    std::vector<std::unique_ptr<AbstractBoundaryCondition>> boundary_conditions;
    //! \brief prediction policy
    PredictionPolicy prediction_policy;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
  };  // end of struct NonLinearEvolutionProblemImplementationBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATIONBASE_HXX */
