/*!
 * \file   include/MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx
 * \brief
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
  struct DirichletBoundaryCondition;
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
     * The `PREDICTION` strategy determines the increment of the
     * unknown \f$\Delta\,{u}\f$ by solving the following linear
     * system:
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
     * - \f$\bts{\mathbb{F}_{e}}\f$ denotes the external forces at the
     *   beginning of the time step.
     * - \f$\bts{\mathbb{F}_{i}}\f$ denotes the inner forces at the beginning
     *   of the time step.
     * - \f$\Delta\,\mathbb{u}\f$ is submitted to the the increment of the
     *   imposed Dirichlet boundary conditions.
     *
     * \note Although the wording explicitly refers to mechanics, this equation
     * applies to all physics.
     */
    BEGINNING_OF_TIME_STEP_PREDICTION,
    /*!
     * The `CONSTANT_GRADIENTS_INTEGRATION` strategy
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
     * - \f$\tilde{\mathbb{K}}\f$ denotes the elastic operator.
     * - \f$\bts{\mathbb{F}_{e}}\f$ denotes the external forces at the
     *   beginning of the time step.
     * - \f$\ets{\tilde{\mathbb{F}}_{i}}\f$ denotes an approximation inner
     *   forces at the end of the time step computed by assuming that the
     *   gradients are constant over the time (and thus egal to their values at
     *   the beginning of the time step).
     * - \f$\Delta\,\mathbb{u}\f$ is submitted to the the increment of the
     *   imposed Dirichlet boundary conditions.
     *
     * \note the behaviour integration allows taking into account:
     * - the evolution of stress-free strain (thermal expansion, swelling,
     * etc..),
     * - the viscoplastic relaxation of the stress (
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
     * time step. See \cite simo_consistent_1985 for details.
     *
     * \note the consistent tangent operator takes into account the details
     * related to the algorithm used to integration the constitutive equations.
     */
    CONSISTENT_TANGENT,
  };

  /*!
   * \brief Prediction policy for estimating the increment of the unknowns over
   * a time step.
   */
  struct [[nodiscard]] PredictionPolicy {
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
     * `PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION`
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
   * \brief operators generated by the linearization of this integrator
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
    static const char *const SolverVerbosityLevel;
    static const char *const SolverRelativeTolerance;
    static const char *const SolverAbsoluteTolerance;
    static const char *const SolverMaximumNumberOfIterations;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual FiniteElementDiscretization &
    getFiniteElementDiscretization() = 0;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual const FiniteElementDiscretization &
    getFiniteElementDiscretization() const = 0;
    //! \return the underlying finite element discretization
    [[nodiscard]] virtual std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() = 0;
    //! \return the unknowns at the beginning of the time step
    [[nodiscard, deprecated("use getUnknowns instead")]]  //
    virtual mfem::Vector &
    getUnknownsAtBeginningOfTheTimeStep() = 0;
    //! \return the unknowns at the beginning of the time step
    [[nodiscard, deprecated("use getUnknowns instead")]]  //
    virtual const mfem::Vector &
    getUnknownsAtBeginningOfTheTimeStep() const = 0;
    //! \return the unknowns at the end of the time step
    [[nodiscard, deprecated("use getUnknowns instead")]]  //
    virtual mfem::Vector &
    getUnknownsAtEndOfTheTimeStep() = 0;
    //! \return the unknowns at the end of the time step
    [[nodiscard, deprecated("use getUnknowns instead")]]  //
    virtual const mfem::Vector &
    getUnknownsAtEndOfTheTimeStep() const = 0;
    //! \return the unknowns at the end of the time step
    [[nodiscard]] virtual mfem::Vector &getUnknowns(
        const TimeStepStage) noexcept = 0;
    //! \return the unknowns at the end of the time step
    [[nodiscard]] virtual const mfem::Vector &getUnknowns(
        const TimeStepStage) const noexcept = 0;
    /*!
     * \brief set the solver parameters
     * \param[in] params: parameters
     *
     * The following parameters are allowed:
     */
    virtual void setSolverParameters(const Parameters &) = 0;
    /*!
     * \brief set the linear solver
     * \param[in] ctx: execution context
     * \param[in] s: linear solver
     */
    [[nodiscard]] virtual bool setLinearSolver(
        Context &, LinearSolverHandler) noexcept = 0;
    /*!
     * \brief set the linear solver
     * \param[in] n: name of the linear solver
     * \param[in] params: parameters
     */
    virtual void setLinearSolver(std::string_view, const Parameters &) = 0;
    //! \return if the stiffness operators from the last iteration are available
    [[nodiscard]] virtual bool areStiffnessOperatorsFromLastIterationAvailable()
        const noexcept = 0;
    /*!
     * \param[in] p: prediction policy
     */
    virtual void setPredictionPolicy(const PredictionPolicy &) noexcept = 0;
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
     * \brief integrate all the behaviours for given estimate of the unknowns at
     * the end of the time step.
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \param[in] odt: optional value for the time step. If invalid, the real
     * time step is used.
     */
    virtual bool integrate(const mfem::Vector &,
                           const IntegrationType,
                           const std::optional<real>) = 0;
    /*!
     * \brief return the linearised operators
     * \note the integrate method shall be call appropriately before calling
     * those operators
     */
    [[nodiscard]] virtual std::optional<LinearizedOperators>
    getLinearizedOperators(Context &, const mfem::Vector &) noexcept = 0;
    //! \brief return the Dirichlet boundary conditions
    [[nodiscard]] virtual const std::vector<
        std::unique_ptr<DirichletBoundaryCondition>>
        &getDirichletBoundaryConditions() const noexcept = 0;
    //! \brief return the standard boundary conditions
    [[nodiscard]] virtual const std::vector<
        std::unique_ptr<AbstractBoundaryCondition>>
        &getBoundaryConditions() const noexcept = 0;
    /*!
     * \brief method called before each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real) = 0;
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual NonLinearResolutionOutput solve(const real, const real) = 0;
    /*!
     * \brief add a new behaviour integrator
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material ids
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string &,
                                        const Parameter &,
                                        const std::string &,
                                        const std::string &) = 0;
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    virtual void setMaterialsNames(
        const std::map<size_type, std::string> &) = 0;
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    virtual void setBoundariesNames(
        const std::map<size_type, std::string> &) = 0;
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    virtual std::vector<size_type> getAssignedMaterialsIdentifiers() const = 0;
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    virtual size_type getMaterialIdentifier(const Parameter &) const = 0;
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    virtual size_type getBoundaryIdentifier(const Parameter &) const = 0;
    /*!
     * \return the list of materials identifiers described by the given
     * parameter.
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     *
     * Integers are directly intepreted as materials identifiers.
     *
     * Strings are intepreted as regular expressions which allows the selection
     * of materials by names.
     */
    virtual std::vector<size_type> getMaterialsIdentifiers(
        const Parameter &) const = 0;
    /*!
     * \return the list of boundaries identifiers described by the given
     * parameter.
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     *
     * Integers are directly intepreted as boundaries identifiers.
     *
     * Strings are intepreted as regular expressions which allows the selection
     * of boundaries by names.
     */
    virtual std::vector<size_type> getBoundariesIdentifiers(
        const Parameter &) const = 0;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual const Material &getMaterial(const Parameter &) const = 0;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual Material &getMaterial(const Parameter &) = 0;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual const AbstractBehaviourIntegrator &getBehaviourIntegrator(
        const size_type) const = 0;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual AbstractBehaviourIntegrator &getBehaviourIntegrator(
        const size_type) = 0;
    /*!
     * \brief add a boundary condition
     * \param[in] f: boundary condition
     */
    [[deprecated]] virtual void addBoundaryCondition(
        std::unique_ptr<AbstractBoundaryCondition>) = 0;
    /*!
     * \brief add a boundary condition
     * \param[in] ctx: execution context
     * \param[in] f: boundary condition
     */
    [[nodiscard]] virtual bool addBoundaryCondition(
        Context &, std::unique_ptr<AbstractBoundaryCondition>) noexcept = 0;
    /*!
     * \brief add a Dirichlet boundary condition
     * \param[in] bc: boundary condition
     */
    virtual void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) = 0;
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(
        const std::function<void(const real, const real)> &) = 0;
    /*!
     * \brief add a new post-processing
     * \param[in] n: name of the post-processing
     * \param[in] p: parameters
     */
    virtual void addPostProcessing(std::string_view, const Parameters &) = 0;
    /*!
     * \brief execute the registred postprocessings
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void executePostProcessings(const real, const real) = 0;
    //! \brief revert the state to the beginning of the time step.
    virtual void revert() = 0;
    //! \brief update the state to the end of the time step.
    virtual void update() = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearEvolutionProblem();
  };  // end of struct AbstractNonLinearEvolutionProblem

  /*!
   * \return the material identifier from the parameters from the `Material`
   * parameter.
   *
   * The `Material` parameter must be either a string and an integer.
   *
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  size_type getMaterialIdentifier(const AbstractNonLinearEvolutionProblem &,
                                  const Parameters &);

  /*!
   * \return the boundary identifier from the parameters from the `Boundary`
   * parameter.
   *
   * The `Boundary` parameter must be either a string and an integer.
   *
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   */
  size_type getBoundaryIdentifier(const AbstractNonLinearEvolutionProblem &,
                                  const Parameters &);

  /*!
   * \return the materials identifiers from the parameters if the one of
   * `Material` or `Materials` parameters exist.  If not such parameter exist,
   * all the materials identifiers are returned if `b` is true, or an error is
   * raised.
   *
   * The `Material` parameter must be either a string and an integer.
   * The `Materials` parameter must be either a string, an integer or a vector
   * of parameters which must be either strings or integers.
   *
   * \note `Material` and `Materials` can't be specificed at the same time.
   *
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Material` or `Materials` parameters
   */
  std::vector<size_type> getMaterialsIdentifiers(
      const AbstractNonLinearEvolutionProblem &,
      const Parameters &,
      const bool = true);

  /*!
   * \return the boundaries identifiers from the parameters if the one of
   * `Boundary` or `Boundaries` parameters exist. If not such parameter exist,
   * all the boundaries identifiers are returned if `b` is true, or an error is
   * raised.
   *
   * The `Boundary` parameter must be either a string and an integer.
   * The `Boundaries` parameter must be either a string, an integer or a
   * vector of parameters which must be either strings or integers.
   *
   * \note `Boundary` and `Boundaries` can't be specificed at the same time.
   *
   * \param[in] p: non linear problem
   * \param[in] params: parameters
   * \param[in] b: allowing missing `Boundary` or `Boundaries` parameters
   */
  std::vector<size_type> getBoundariesIdentifiers(
      const AbstractNonLinearEvolutionProblem &,
      const Parameters &,
      const bool = true);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX */
