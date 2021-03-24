/*!
 * \file   src/NonLinearEvolutionProblemImplementation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/PostProcessingFactory.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  static void setSolverParametersImplementation(NewtonSolver& solver,
                                                const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    if (contains(params, Problem::SolverVerbosityLevel)) {
      solver.SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
    }
    if (contains(params, Problem::SolverRelativeTolerance)) {
      solver.SetRelTol(get<double>(params, Problem::SolverRelativeTolerance));
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      solver.SetAbsTol(get<double>(params, Problem::SolverAbsoluteTolerance));
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      solver.SetMaxIter(
          get<int>(params, Problem::SolverMaximumNumberOfIterations));
    }
  }  // end of setSolverParametersImplementation

  /*!
   * \brief post-processing defined by an std::function
   * \tparam parallel: boolean stating if the computations are performed in
   * parallel.
   */
  template <bool parallel>
  struct StdFunctionPostProcessing final : PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] fct: function executing the postprocessing
     */
    StdFunctionPostProcessing(
        const std::function<void(const real, const real)>& fct)
        : f(fct) {}  // end of StdFunctionPostProcessing
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real t,
                 const real dt) override {
      this->f(t, dt);
    }  // end of execute
    //! \brief destructor
    ~StdFunctionPostProcessing() override = default;

   private:
    //! \brief function executing the post-processing
    std::function<void(const real, const real)> f;
  };  // end of struct StdFunctionPostProcessing

#ifdef MFEM_USE_MPI

  NonLinearEvolutionProblemImplementation<true>::
      NonLinearEvolutionProblemImplementation(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : NonLinearEvolutionProblemImplementationBase(fed, h, p),
        mfem::ParNonlinearForm(&(fed->getFiniteElementSpace<true>())),
        solver(fed->getFiniteElementSpace<true>().GetComm()) {
    if (this->fe_discretization->getMesh<true>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      std::unique_ptr<PostProcessing<true>> p) {
    this->postprocessings.push_back(std::move(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      std::string_view n, const Parameters& p) {
    const auto& f = PostProcessingFactory<true>::getFactory();
    this->postprocessings.push_back(f.generate(n, *this, p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      const std::function<void(const real, const real)>& p) {
    this->addPostProcessing(
        std::make_unique<StdFunctionPostProcessing<true>>(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<true>::executePostProcessings(
      const real t, const real dt) {
    for (auto& p : this->postprocessings) {
      p->execute(*this, t, dt);
    }
  }  // end of executePostProcessings

  const FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  void NonLinearEvolutionProblemImplementation<true>::setSolverParameters(
      const Parameters& params) {
    setSolverParametersImplementation(this->solver, params);
  }  // end of setSolverParameters

  NewtonSolver& NonLinearEvolutionProblemImplementation<true>::getSolver() {
    return this->solver;
  }  // end of getSolver

  void NonLinearEvolutionProblemImplementation<true>::solve(const real t,
                                                            const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->setup(t, dt);
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  void NonLinearEvolutionProblemImplementation<true>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  NonLinearEvolutionProblemImplementation<
      true>::~NonLinearEvolutionProblemImplementation() = default;
#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblemImplementation<false>::
      NonLinearEvolutionProblemImplementation(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : NonLinearEvolutionProblemImplementationBase(fed, h, p),
        mfem::NonlinearForm(&(fed->getFiniteElementSpace<false>())) {
    if (this->fe_discretization->getMesh<false>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  void NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      std::unique_ptr<PostProcessing<false>> p) {
    this->postprocessings.push_back(std::move(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      const std::function<void(const real, const real)>& p) {
    this->addPostProcessing(
        std::make_unique<StdFunctionPostProcessing<false>>(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<false>::addPostProcessing(
      std::string_view n, const Parameters& p) {
    const auto& f = PostProcessingFactory<false>::getFactory();
    this->postprocessings.push_back(f.generate(n, *this, p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<false>::executePostProcessings(
      const real t, const real dt) {
    for (auto& p : this->postprocessings) {
      p->execute(*this, t, dt);
    }
  }  // end of executePostProcessings

  const FiniteElementSpace<false>& NonLinearEvolutionProblemImplementation<
      false>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<false>&
  NonLinearEvolutionProblemImplementation<false>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  NewtonSolver& NonLinearEvolutionProblemImplementation<false>::getSolver() {
    return this->solver;
  }  // end of getSolver

  void NonLinearEvolutionProblemImplementation<false>::setSolverParameters(
      const Parameters& params) {
    setSolverParametersImplementation(this->solver, params);
  }  // end of setSolverParameters

  void NonLinearEvolutionProblemImplementation<false>::solve(const real t,
                                                             const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->setup(t, dt);
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  void NonLinearEvolutionProblemImplementation<false>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  NonLinearEvolutionProblemImplementation<
      false>::~NonLinearEvolutionProblemImplementation() = default;

}  // end of namespace mfem_mgis
