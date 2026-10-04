/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief  This file declares the `NonLinearEvolutionProblem` class
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX

#include <vector>
#include <utility>

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem
      : AbstractNonLinearEvolutionProblem {
    //! \brief name of the `Hypothesis` parameter
    static const char *const HypothesisParameter;
    //! \return the list of valid parameters
    static std::vector<std::string> getParametersList();
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: parameters.
     *
     * The following parameters are the most common (see also the
     * FiniteElementDiscretization class for details):
     *
     * - `Parallel` (boolean): if true, a parallel computation is to be
     *    performed. This value is assumed to be false by default.
     * - `MeshFileName` (string): mesh file.
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. The default value is `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `Hypothesis` (string): modelling hypothesis
     * - `UseMultiMaterialNonLinearIntegrator` (boolean): if false, do not
     *   add the `MultiMaterialNonLinearIntegrator`. True by default.
     */
    NonLinearEvolutionProblem(Context &ctx, const Parameters &p);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] p: parameters used to initialize the problem and the
     * underlying finite element discretization
     *
     * The following parameters are the most common (see also the
     * FiniteElementDiscretization class for details):
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. The default value is `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `Hypothesis` (string): modelling hypothesis
     * - `UseMultiMaterialNonLinearIntegrator` (boolean): if false, do not
     *   add the `MultiMaterialNonLinearIntegrator`. True by default.
     */
    NonLinearEvolutionProblem(Context &ctx,
                              MeshDiscretization &m,
                              const Parameters &p);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters used to initialize the problem and the
     * underlying finite element discretization
     *
     * The following parameters are the most common (see also the
     * FiniteElementDiscretization class for details):
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. The default value is `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UseMultiMaterialNonLinearIntegrator` (boolean): if false, do not
     *   add the `MultiMaterialNonLinearIntegrator`. True by default.
     */
    NonLinearEvolutionProblem(Context &ctx,
                              MeshDiscretization &m,
                              const Hypothesis h,
                              const Parameters &p = Parameters());
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblem(Context &ctx,
                              std::shared_ptr<FiniteElementDiscretization> fed,
                              const Parameters &p);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblem(Context &ctx,
                              std::shared_ptr<FiniteElementDiscretization> fed,
                              const Hypothesis h,
                              const Parameters &p = Parameters());
    //! \brief move constructor
    NonLinearEvolutionProblem(NonLinearEvolutionProblem &&) noexcept = default;
    //! \brief deleted copy constructor
    NonLinearEvolutionProblem(const NonLinearEvolutionProblem &) noexcept =
        delete;
    /*!
     * \return the internal implementation
     * \warning this method shall be used with care: calling it if the
     * implementation does not match the parallel template argument aborts the
     * code
     */
    template <bool parallel>
    [[nodiscard]] NonLinearEvolutionProblemImplementation<parallel>
        &getImplementation() noexcept;
    /*!
     * \return the internal implementation
     * \warning this method shall be used with care: calling it if the
     * implementation does not match the parallel template argument aborts the
     * code
     */
    template <bool parallel>
    [[nodiscard]] const NonLinearEvolutionProblemImplementation<parallel>
        &getImplementation() const noexcept;
    //
    [[nodiscard]] FiniteElementDiscretization &
    getFiniteElementDiscretization() noexcept override;
    [[nodiscard]] const FiniteElementDiscretization &
    getFiniteElementDiscretization() const noexcept override;
    [[nodiscard]] std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() noexcept override;
    [[nodiscard]] bool setMaterialsNames(
        Context &ctx,
        const std::map<size_type, std::string> &ids) noexcept override;
    [[nodiscard]] bool setBoundariesNames(
        Context &ctx,
        const std::map<size_type, std::string> &ids) noexcept override;
    [[nodiscard]] mfem::Vector &getUnknowns(
        const TimeStepStage ts) noexcept override;
    [[nodiscard]] const mfem::Vector &getUnknowns(
        const TimeStepStage ts) const noexcept override;
    [[nodiscard]] bool setSolverParameters(
        Context &ctx, const Parameters &params) noexcept override;
    [[nodiscard]] bool setLinearSolver(Context &ctx,
                                       LinearSolverHandler s) noexcept override;
    [[nodiscard]] bool setLinearSolver(
        Context &ctx,
        std::string_view n,
        const Parameters &params) noexcept override;
    [[nodiscard]] bool areStiffnessOperatorsFromLastIterationAvailable()
        const noexcept override;
    void setPredictionPolicy(const PredictionPolicy &p) noexcept override;
    [[nodiscard]] PredictionPolicy getPredictionPolicy()
        const noexcept override;
    [[nodiscard]] bool addBoundaryCondition(
        Context &ctx,
        std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept
        override;
    [[nodiscard]] bool addBoundaryCondition(
        Context &ctx,
        std::unique_ptr<AbstractBoundaryCondition> f) noexcept override;
    /*!
     * \brief add a uniform Dirichlet boundary condition
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters defining the boundary condition
     * \return true on success
     */
    [[nodiscard]] bool addUniformDirichletBoundaryCondition(
        Context &ctx, const Parameters &params) noexcept;
    [[nodiscard]] bool addPostProcessing(
        Context &ctx,
        const std::function<void(const real, const real)> &p) noexcept override;
    [[nodiscard]] bool addPostProcessing(Context &ctx,
                                         std::string_view n,
                                         const Parameters &p) noexcept override;
    [[nodiscard]] bool executeInitialPostProcessings(
        Context &ctx, const real t) noexcept override;
    [[nodiscard]] bool executePostProcessings(Context &ctx,
                                              const real t,
                                              const real dt) noexcept override;
    std::optional<std::map<size_type, size_type>> addBehaviourIntegrator(
        Context &ctx,
        const std::string &n,
        const Parameter &m,
        const std::string &l,
        const std::string &b) noexcept override;
    std::optional<std::map<size_type, size_type>> addBehaviourIntegrator(
        Context &ctx,
        const std::string &n,
        const Parameter &m,
        const std::string &l,
        const std::string &b,
        const Parameters &params) noexcept override;
    [[nodiscard]] std::vector<size_type> getAssignedMaterialsIdentifiers()
        const noexcept override;
    [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
        Context &ctx, const Parameter &m) const noexcept override;
    [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
        Context &ctx, const Parameter &m) const noexcept override;
    [[nodiscard]] std::optional<std::vector<size_type>> getMaterialsIdentifiers(
        Context &ctx, const Parameter &m) const noexcept override;
    [[nodiscard]] std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context &ctx,
                             const Parameter &m) const noexcept override;
    OptionalReference<const Material> getMaterial(
        Context &ctx,
        const Parameter &m,
        const size_type b) const noexcept override;
    OptionalReference<Material> getMaterial(
        Context &ctx, const Parameter &m, const size_type b) noexcept override;
    [[nodiscard]] std::optional<size_type> getNumberOfBehaviourIntegrators(
        Context &ctx, const Parameter &m) const noexcept override;
    OptionalReference<const AbstractBehaviourIntegrator> getBehaviourIntegrator(
        Context &ctx,
        const Parameter &m,
        const size_type b) const noexcept override;
    OptionalReference<AbstractBehaviourIntegrator> getBehaviourIntegrator(
        Context &ctx, const Parameter &m, const size_type b) noexcept override;
    std::vector<size_type> getEssentialDegreesOfFreedom() const override;
    [[nodiscard]] bool integrate(const mfem::Vector &u,
                                 const IntegrationType it,
                                 const std::optional<real> odt) override;
    [[nodiscard]] std::optional<LinearizedOperators> getLinearizedOperators(
        Context &ctx, const mfem::Vector &U) noexcept override;
    [[nodiscard]] const std::vector<
        std::unique_ptr<AbstractDirichletBoundaryCondition>>
        &getDirichletBoundaryConditions() const noexcept override;
    [[nodiscard]] virtual const std::vector<
        std::unique_ptr<AbstractBoundaryCondition>>
        &getBoundaryConditions() const noexcept override;
    /*!
     * \brief call `setup`, then solve the non linear problem over the given
     * time step
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return the output of the non linear resolution
     */
    [[nodiscard]] NonLinearResolutionOutput solve(
        Context &ctx, const real t, const real dt) noexcept override;
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   protected:
    /*!
     * \brief method called before each resolution
     * This method is automatically called by `solve`
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool setup(Context &ctx,
                                     const real t,
                                     const real dt) noexcept;
    //! \brief implementation of the non linear problem
    std::unique_ptr<AbstractNonLinearEvolutionProblem> pimpl;
  };  // end of struct NonLinearEvolutionProblem

  /*!
   * \brief resolve the dependencies (material properties and external state
   * variables) of the first nonlinear evolution problem given using the
   * gradients, thermodynamic forces and internal state variables from the
   * second.
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] p: problem to be treated
   * \param[in] provider: problem used to resolve the dependencies
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool resolveBehaviourIntegratorsDependencies(
      Context &ctx,
      NonLinearEvolutionProblem &p,
      const NonLinearEvolutionProblem &provider) noexcept;

#ifdef MFEM_USE_MPI

  //! \brief parallel specialisation
  template <>
  NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() noexcept;

  //! \brief parallel specialisation
  template <>
  const NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() const noexcept;

#else /* MFEM_USE_MPI */

  //! \brief parallel specialisation
  template <>
  [[noreturn]] NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() noexcept;

  //! \brief parallel specialisation
  template <>
  [[noreturn]] const NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() const noexcept;

#endif /* MFEM_USE_MPI */

  //! \brief sequential specialisation
  template <>
  NonLinearEvolutionProblemImplementation<false>
      &NonLinearEvolutionProblem::getImplementation() noexcept;

  //! \brief sequential specialisation
  template <>
  const NonLinearEvolutionProblemImplementation<false>
      &NonLinearEvolutionProblem::getImplementation() const noexcept;

  /*!
   * \brief describe a boundary by its faces
   * \return a description of the boundary by a vector of pairs
   * associating for each face its identifier and the identifier of the
   * adjacent element.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<std::pair<size_type, size_type>>
  buildFacesDescription(NonLinearEvolutionProblem &p, const size_type bid);

  /*!
   * \brief list the elements having degrees of freedom on a boundary
   * \return a vector of pairs associating the index of each element having
   * at least one degree of freedom on the boundary with the local indexes of
   * these degrees of freedom, grouped by component.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   * \note in parallel, all processes must call this function: the degrees of
   * freedom of the boundary are synchronized between them.
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::vector<
      std::pair<size_type, std::vector<std::vector<size_type>>>>
  getElementsDegreesOfFreedomOnBoundary(NonLinearEvolutionProblem &p,
                                        const size_type bid);

  /*!
   * \brief compute the resultant of the inner forces on the given boundary
   * \param[in, out] ctx: execution context
   * \param[out] F: resultant
   * \param[in] p: non linear evolution problem
   * \param[in] elts_dofs: a structure which gives for each element having at
   * least one node on the boundary the list of the nodes of this element on the
   * boundary.
   * \param[in] selection: selection of the behaviour integrators of each
   * material, whose inner forces are summed
   * \return true on success
   *
   * \note in parallel, the resultant is only the contribution of the given
   * process
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeResultantForceOnBoundary(
      Context &ctx,
      mfem::Vector &F,
      NonLinearEvolutionProblem &p,
      const std::vector<
          std::pair<size_type, std::vector<std::vector<size_type>>>>
          &elts_dofs,
      const BehaviourIntegratorsSelection &selection = {}) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX */
