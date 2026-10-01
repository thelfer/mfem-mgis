/*!
 * \file   include/MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx
 * \brief  This file declares the `AbstractNonLinearEvolutionProblem` class
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX

#include <map>
#include <string>
#include <memory>
#include <vector>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/LinearSolverHandler.hxx"
#include "MFEMMGIS/NonLinearResolutionOutput.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameter;
  struct Parameters;
  struct FiniteElementDiscretization;
  struct AbstractDirichletBoundaryCondition;
  struct AbstractBehaviourIntegrator;
  struct Material;
  struct AbstractBoundaryCondition;
  enum struct IntegrationType;

  /*!
   * \brief strategy used to make a prediction of the solution at the end of the
   * time step
   */
  enum struct PredictionStrategy {
    /*!
     * \brief By default, a nonlinear evolution problem uses the solution at the
     * beginning of the time step, modified by applying Dirichlet boundary
     * conditions, as the initial guess of the solution at the end of the time
     * step.
     *
     * \warning In mechanics, this may lead to very high increments of the
     * deformation gradients or the strain in the neighboring elements of
     * boundaries where evolving unknowns are imposed.
     */
    DEFAULT_PREDICTION,
    /*!
     * The `BEGINNING_OF_TIME_STEP_PREDICTION` strategy determines the
     * increment of the unknown \f$\Delta\,{u}\f$ by solving the following
     * linear system:
     *
     * \f[
     *   \mathbb{K}\,\cdot\,\Delta\,\mathbb{u} =
     *   \ets{\mathbb{F}_{e}}-\bts{\mathbb{F}_{i}}
     * \f]
     *
     * where:
     *
     * - \f$\mathbb{K}\f$ denotes the prediction operator (See
     * `PredictionPolicy::prediction_operator`)
     * - \f$\ets{\mathbb{F}_{e}}\f$ denotes the external forces at the
     *   end of the time step.
     * - \f$\bts{\mathbb{F}_{i}}\f$ denotes the inner forces at the beginning
     *   of the time step.
     * - \f$\Delta\,\mathbb{u}\f$ is subject to the increment of the
     *   imposed Dirichlet boundary conditions.
     *
     * \note Although the wording explicitly refers to mechanics, this equation
     * applies to all physics.
     */
    BEGINNING_OF_TIME_STEP_PREDICTION,
    /*!
     * The `CONSTANT_GRADIENTS_INTEGRATION_PREDICTION` strategy
     * determines the increment of the unknown \f$\Delta\,{u}\f$ by solving
     * the following linear system:
     *
     * \f[
     *   \tilde{\mathbb{K}}\,\cdot\,\Delta\,\mathbb{u} =
     *   \ets{\mathbb{F}_{e}}-\ets{\tilde{\mathbb{F}}_{i}}
     * \f]
     *
     * where:
     *
     * - \f$\tilde{\mathbb{K}}\f$ denotes the operator selected by
     *   `PredictionPolicy::integration_operator`.
     * - \f$\ets{\mathbb{F}_{e}}\f$ denotes the external forces at the
     *   end of the time step.
     * - \f$\ets{\tilde{\mathbb{F}}_{i}}\f$ denotes an approximation of the
     *   inner forces at the end of the time step computed by assuming that the
     *   gradients are constant over the time (and thus equal to their values
     *   at the beginning of the time step).
     * - \f$\Delta\,\mathbb{u}\f$ is subject to the increment of the
     *   imposed Dirichlet boundary conditions.
     *
     * \note the behaviour integration allows taking into account:
     * - the evolution of stress-free strain (thermal expansion, swelling,
     * etc.),
     * - the viscoplastic relaxation of the stress.
     *
     * \note Although the wording explicitly refers to mechanics, this equation
     * applies to all physics.
     */
    CONSTANT_GRADIENTS_INTEGRATION_PREDICTION
  };

  /*!
   * \brief list of prediction operators available
   */
  enum struct PredictionOperator {
    //! \brief The elastic operator
    ELASTIC,
    /*!
     * \brief The secant operator is typically defined by the elastic operator
     * affected by damage.
     *
     * \note Most behaviours do not compute the secant operator or simply return
     * the elastic operator.
     */
    SECANT,
    /*!
     * \brief The tangent operator, defined by the time-continuous derivative of
     * the thermodynamic force with respect to the gradients.
     *
     * \note Most behaviours do not compute the tangent operator.
     */
    TANGENT,
    /*!
     * The `LAST_ITERATE_OPERATOR` relies on the operator
     * computed at the last iteration of the previous time step.
     *
     * \note At the first time step, the elastic operator is used.
     * \warning not implemented yet: the elastic operator is always used.
     */
    LAST_ITERATE_OPERATOR
  };

  /*!
   * \brief list of operators available after the behaviour integration
   */
  enum struct IntegrationOperator {
    //! \brief The elastic operator
    ELASTIC,
    /*!
     * \brief The secant operator is typically defined by the elastic operator
     * affected by damage.
     *
     * \note Most behaviours do not compute the secant operator or simply return
     * the elastic operator.
     */
    SECANT,
    /*!
     * \brief The tangent operator, defined by the time-continuous derivative of
     * the thermodynamic force with respect to the gradients.
     *
     * \note Most behaviours do not compute the tangent operator.
     */
    TANGENT,
    /*!
     * \brief The consistent tangent operator, defined by the derivative of
     * the thermodynamic force with respect to the gradients at the end of the
     * time step. See Simo and Taylor, Consistent tangent operators for
     * rate-independent elastoplasticity, 1985.
     *
     * \note the consistent tangent operator takes into account the details
     * related to the algorithm used to integrate the constitutive equations.
     */
    CONSISTENT_TANGENT
  };

  /*!
   * \brief Prediction policy for estimating the increment of the unknowns over
   * a time step.
   */
  struct [[nodiscard]] PredictionPolicy {
    //! \brief selected prediction strategy
    PredictionStrategy strategy = PredictionStrategy::DEFAULT_PREDICTION;
    /*!
     * \brief selected prediction operator
     *
     * \note this member is only used if `strategy` is set to
     * `PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION`
     */
    PredictionOperator prediction_operator = PredictionOperator::ELASTIC;
    /*!
     * \brief selected operator
     *
     * \note this member is only used if `strategy` is set to
     * `PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`
     */
    IntegrationOperator integration_operator =
        IntegrationOperator::CONSISTENT_TANGENT;
    /*!
     * \brief boolean stating if the behaviour shall be integrated using a null
     * time increment.
     *
     * \note a null time increment allows neglecting viscoplastic relaxation.
     */
    bool null_time_increment = false;
  };  // end of PredictionPolicy

  /*!
   * \brief operators generated by the linearization of the behaviour
   * integrators
   */
  struct [[nodiscard]] LinearizedOperators {
    //! \brief stiffness matrix
    std::unique_ptr<BilinearFormIntegrator> K;
    //! \brief opposite of the inner forces
    std::unique_ptr<LinearFormIntegrator> mFi;
  };

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT AbstractNonLinearEvolutionProblem {
    //! \brief string associated to the `VerbosityLevel` parameter
    static const char *const SolverVerbosityLevel;
    //! \brief string associated to the `RelativeTolerance` parameter
    static const char *const SolverRelativeTolerance;
    //! \brief string associated to the `AbsoluteTolerance` parameter
    static const char *const SolverAbsoluteTolerance;
    //! \brief string associated to the `MaximumNumberOfIterations` parameter
    static const char *const SolverMaximumNumberOfIterations;
    /*!
     * \brief set the names of the materials
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     * \return true on success
     */
    [[nodiscard]] virtual bool setMaterialsNames(
        Context &ctx, const std::map<size_type, std::string> &ids) noexcept = 0;
    /*!
     * \brief set the names of the boundaries
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     * \return true on success
     */
    [[nodiscard]] virtual bool setBoundariesNames(
        Context &ctx, const std::map<size_type, std::string> &ids) noexcept = 0;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual FiniteElementDiscretization &
    getFiniteElementDiscretization() noexcept = 0;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual const FiniteElementDiscretization &
    getFiniteElementDiscretization() const noexcept = 0;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() noexcept = 0;
    /*!
     * \brief return the unknowns at the given time step stage
     * \param[in] ts: time step stage
     * \return the unknowns at the given time step stage
     */
    [[nodiscard]] virtual mfem::Vector &getUnknowns(
        const TimeStepStage ts) noexcept = 0;
    /*!
     * \brief return the unknowns at the given time step stage
     * \param[in] ts: time step stage
     * \return the unknowns at the given time step stage
     */
    [[nodiscard]] virtual const mfem::Vector &getUnknowns(
        const TimeStepStage ts) const noexcept = 0;
    /*!
     * \brief set the solver parameters
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters
     * \return true on success
     */
    [[nodiscard]] virtual bool setSolverParameters(
        Context &ctx, const Parameters &params) noexcept = 0;
    /*!
     * \brief set the linear solver
     * \param[in, out] ctx: execution context
     * \param[in] s: linear solver
     * \return true on success
     */
    [[nodiscard]] virtual bool setLinearSolver(
        Context &ctx, LinearSolverHandler s) noexcept = 0;
    /*!
     * \brief set the linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] params: parameters
     * \return true on success
     */
    [[nodiscard]] virtual bool setLinearSolver(
        Context &ctx,
        std::string_view n,
        const Parameters &params) noexcept = 0;
    //! \return if the stiffness operators from the last iteration are available
    [[nodiscard]] virtual bool areStiffnessOperatorsFromLastIterationAvailable()
        const noexcept = 0;
    /*!
     * \brief set the prediction policy
     * \param[in] p: prediction policy
     */
    virtual void setPredictionPolicy(const PredictionPolicy &p) noexcept = 0;
    //! \return the prediction policy
    [[nodiscard]] virtual PredictionPolicy getPredictionPolicy()
        const noexcept = 0;
    /*!
     * \return the list of the degrees of freedom handled by Dirichlet boundary
     * conditions.
     */
    [[nodiscard]] virtual std::vector<size_type> getEssentialDegreesOfFreedom()
        const = 0;
    /*!
     * \brief integrate all the behaviours for a given estimate of the unknowns
     * at the end of the time step.
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \param[in] odt: optional value for the time step. If invalid, the real
     * time step is used.
     * \return true on success
     */
    virtual bool integrate(const mfem::Vector &u,
                           const IntegrationType it,
                           const std::optional<real> odt) = 0;
    /*!
     * \brief return the linearised operators
     * \param[in, out] ctx: execution context
     * \param[in] U: current estimate of the unknowns
     * \return the linearised operators
     * \note the integrate method shall be called appropriately before calling
     * those operators
     */
    [[nodiscard]] virtual std::optional<LinearizedOperators>
    getLinearizedOperators(Context &ctx, const mfem::Vector &U) noexcept = 0;
    //! \return the Dirichlet boundary conditions
    [[nodiscard]] virtual const std::vector<
        std::unique_ptr<AbstractDirichletBoundaryCondition>>
        &getDirichletBoundaryConditions() const noexcept = 0;
    //! \return the standard boundary conditions
    [[nodiscard]] virtual const std::vector<
        std::unique_ptr<AbstractBoundaryCondition>>
        &getBoundaryConditions() const noexcept = 0;
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return the output of the non linear resolution
     */
    [[nodiscard]] virtual NonLinearResolutionOutput solve(
        Context &ctx, const real t, const real dt) noexcept = 0;
    /*!
     * \brief add a new behaviour integrator
     * \return a mapping between the material id and the identifier of the
     * behaviour integrator
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material ids
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual std::optional<std::map<size_type, size_type>>
    addBehaviourIntegrator(Context &ctx,
                           const std::string &n,
                           const Parameter &m,
                           const std::string &l,
                           const std::string &b) noexcept = 0;
    /*!
     * \brief add a new behaviour integrator
     * \return a mapping between the material id and the identifier of the
     * behaviour integrator
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material ids
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] params: additional parameters
     */
    virtual std::optional<std::map<size_type, size_type>>
    addBehaviourIntegrator(Context &ctx,
                           const std::string &n,
                           const Parameter &m,
                           const std::string &l,
                           const std::string &b,
                           const Parameters &params) noexcept = 0;
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    virtual std::vector<size_type> getAssignedMaterialsIdentifiers()
        const noexcept = 0;
    /*!
     * \brief return the material identifier described by the given parameter
     * \return the material identifier described by the given parameter.
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     *
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] virtual std::optional<size_type> getMaterialIdentifier(
        Context &ctx, const Parameter &m) const noexcept = 0;
    /*!
     * \brief return the boundary identifier described by the given parameter
     * \return the boundary identifier described by the given parameter.
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: boundary name or identifier
     *
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] virtual std::optional<size_type> getBoundaryIdentifier(
        Context &ctx, const Parameter &m) const noexcept = 0;
    /*!
     * \brief return the materials identifiers described by the given parameter
     * \return the list of materials identifiers described by the given
     * parameter.
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings or integers.
     *
     * Integers are directly interpreted as materials identifiers.
     *
     * Strings are interpreted as regular expressions which allows the
     * selection of materials by names.
     */
    [[nodiscard]] virtual std::optional<std::vector<size_type>>
    getMaterialsIdentifiers(Context &ctx,
                            const Parameter &m) const noexcept = 0;
    /*!
     * \brief return the boundaries identifiers described by the given
     * parameter
     * \return the list of boundaries identifiers described by the given
     * parameter.
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: boundary name or identifier
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings or integers.
     *
     * Integers are directly interpreted as boundaries identifiers.
     *
     * Strings are interpreted as regular expressions which allows the
     * selection of boundaries by names.
     */
    [[nodiscard]] virtual std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context &ctx,
                             const Parameter &m) const noexcept = 0;
    /*!
     * \brief get the material
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     * \param[in] b: behaviour integrator id
     */
    virtual OptionalReference<const Material> getMaterial(
        Context &ctx, const Parameter &m, const size_type b) const noexcept = 0;
    /*!
     * \brief get the material
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     * \param[in] b: behaviour integrator id
     */
    virtual OptionalReference<Material> getMaterial(
        Context &ctx, const Parameter &m, const size_type b) noexcept = 0;
    /*!
     * \brief get the number of behaviour integrators of a material
     * \return the number of behaviour integrators for the given material
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     */
    [[nodiscard]] virtual std::optional<size_type>
    getNumberOfBehaviourIntegrators(Context &ctx,
                                    const Parameter &m) const noexcept = 0;
    /*!
     * \brief get a behaviour integrator
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     * \param[in] b: behaviour integrator id
     */
    virtual OptionalReference<const AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context &ctx,
                           const Parameter &m,
                           const size_type b) const noexcept = 0;
    /*!
     * \brief get a behaviour integrator
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material name or identifier
     * \param[in] b: behaviour integrator id
     */
    virtual OptionalReference<AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context &ctx,
                           const Parameter &m,
                           const size_type b) noexcept = 0;
    /*!
     * \brief add a boundary condition
     * \param[in, out] ctx: execution context
     * \param[in] f: boundary condition
     * \return true on success
     */
    [[nodiscard]] virtual bool addBoundaryCondition(
        Context &ctx,
        std::unique_ptr<AbstractBoundaryCondition> f) noexcept = 0;
    /*!
     * \brief add a Dirichlet boundary condition
     * \param[in, out] ctx: execution context
     * \param[in] bc: boundary condition
     * \return true on success
     */
    [[nodiscard]] virtual bool addBoundaryCondition(
        Context &ctx,
        std::unique_ptr<AbstractDirichletBoundaryCondition> bc) = 0;
    /*!
     * \brief add a new post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: post-processing
     * \return true on success
     */
    [[nodiscard]] virtual bool addPostProcessing(
        Context &ctx,
        const std::function<void(const real, const real)> &p) noexcept = 0;
    /*!
     * \brief add a new post-processing
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the post-processing
     * \param[in] p: parameters
     * \return true on success
     */
    [[nodiscard]] virtual bool addPostProcessing(
        Context &ctx, std::string_view n, const Parameters &p) noexcept = 0;
    /*!
     * \brief execute the registered postprocessings at the initial time of the
     * simulation
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] virtual bool executeInitialPostProcessings(
        Context &ctx, const real t) noexcept = 0;
    /*!
     * \brief execute the registered postprocessings
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool executePostProcessings(
        Context &ctx, const real t, const real dt) noexcept = 0;
    /*!
     * \brief revert the state to the beginning of the time step.
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] virtual bool revert(Context &ctx) noexcept = 0;
    /*!
     * \brief update the state to the end of the time step.
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] virtual bool update(Context &ctx) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearEvolutionProblem();
  };  // end of struct AbstractNonLinearEvolutionProblem

  /*!
   * \brief get the material identifier from the `Material` parameter
   * \return the material identifier from the parameters from the `Material`
   * parameter.
   *
   * The `Material` parameter must be either a string or an integer.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params) noexcept;

  /*!
   * \brief get the boundary identifier from the `Boundary` parameter
   * \return the boundary identifier from the parameters from the `Boundary`
   * parameter.
   *
   * The `Boundary` parameter must be either a string or an integer.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
      Context &ctx,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params) noexcept;

  /*!
   * \brief get the materials identifiers from the parameters
   * \return the materials identifiers from the parameters if one of the
   * `Material` or `Materials` parameters exists. If no such parameter exists,
   * the identifiers of the materials with a behaviour integrator are returned
   * if `b` is true, or an error is raised.
   *
   * The `Material` parameter must be either a string or an integer.
   * The `Materials` parameter must be either a string, an integer or a vector
   * of parameters which must be either strings or integers.
   *
   * \note `Material` and `Materials` can't be specified at the same time.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Material` or `Materials` parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<size_type>>
  getMaterialsIdentifiers(Context &ctx,
                          const AbstractNonLinearEvolutionProblem &p,
                          const Parameters &params,
                          const bool b = true) noexcept;

  /*!
   * \brief get the boundaries identifiers from the parameters
   * \return the boundaries identifiers from the parameters if one of the
   * `Boundary` or `Boundaries` parameters exists. If no such parameter
   * exists, the identifiers of all the named boundaries are returned if `b`
   * is true, or an error is raised.
   *
   * The `Boundary` parameter must be either a string or an integer.
   * The `Boundaries` parameter must be either a string, an integer or a
   * vector of parameters which must be either strings or integers.
   *
   * \note `Boundary` and `Boundaries` can't be specified at the same time.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Boundary` or `Boundaries` parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<size_type>>
  getBoundariesIdentifiers(Context &ctx,
                           const AbstractNonLinearEvolutionProblem &p,
                           const Parameters &params,
                           const bool b = true) noexcept;

  /*!
   * \brief get the material identifier from the `Material` parameter
   * \return the material identifier from the parameters from the `Material`
   * parameter.
   *
   * The `Material` parameter must be either a string or an integer.
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getMaterialIdentifier(
      attributes::Throwing throwing,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params);

  /*!
   * \brief get the boundary identifier from the `Boundary` parameter
   * \return the boundary identifier from the parameters from the `Boundary`
   * parameter.
   *
   * The `Boundary` parameter must be either a string or an integer.
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getBoundaryIdentifier(
      attributes::Throwing throwing,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params);

  /*!
   * \brief get the materials identifiers from the parameters
   * \return the materials identifiers from the parameters if one of the
   * `Material` or `Materials` parameters exists. If no such parameter exists,
   * the identifiers of the materials with a behaviour integrator are returned
   * if `b` is true, or an error is raised.
   *
   * The `Material` parameter must be either a string or an integer.
   * The `Materials` parameter must be either a string, an integer or a vector
   * of parameters which must be either strings or integers.
   *
   * \note `Material` and `Materials` can't be specified at the same time.
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Material` or `Materials` parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type> getMaterialsIdentifiers(
      attributes::Throwing throwing,
      const AbstractNonLinearEvolutionProblem &p,
      const Parameters &params,
      const bool b = true);

  /*!
   * \brief get the boundaries identifiers from the parameters
   * \return the boundaries identifiers from the parameters if one of the
   * `Boundary` or `Boundaries` parameters exists. If no such parameter
   * exists, the identifiers of all the named boundaries are returned if `b`
   * is true, or an error is raised.
   *
   * The `Boundary` parameter must be either a string or an integer.
   * The `Boundaries` parameter must be either a string, an integer or a
   * vector of parameters which must be either strings or integers.
   *
   * \note `Boundary` and `Boundaries` can't be specified at the same time.
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Boundary` or `Boundaries` parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type>
  getBoundariesIdentifiers(attributes::Throwing throwing,
                           const AbstractNonLinearEvolutionProblem &p,
                           const Parameters &params,
                           const bool b = true);

#ifdef MFEM_USE_MPI

  /*!
   * \brief get the MPI communicator of the nonlinear evolution problem
   * \return the MPI communicator associated with the nonlinear evolution
   * problem
   *
   * \param[in] p: nonlinear evolution problem
   *
   * \note If a sequential computation is described, `MPI_COMM_WORLD` is
   * returned.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] MPI_Comm getMPICommunicator(
      const AbstractNonLinearEvolutionProblem &p) noexcept;
  /*!
   * \brief check if the current process is the main one
   * \return if the current process is the main one (the process of rank 0)
   * \param[in] p: nonlinear evolution problem
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isMainProcess(
      const AbstractNonLinearEvolutionProblem &p) noexcept;

#endif /* MFEM_USE_MPI */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX */
