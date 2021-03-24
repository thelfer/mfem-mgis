/*!
 * \file   src/NonLinearEvolutionProblemImplementation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

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
          const Parameters&)
      : NonLinearEvolutionProblemImplementationBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h),
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
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblemImplementation

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      std::unique_ptr<PostProcessing<true>> p) {
    this->postprocessings.push_back(std::move(p));
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

  void NonLinearEvolutionProblemImplementation<true>::addBehaviourIntegrator(
      const std::string& n,
      const size_type l,
      const std::string& m,
      const std::string& b) {
    MultiMaterialEvolutionProblemBase::addBehaviourIntegrator(n, l, m, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblemImplementation<true>::getMaterial(
      const size_type m) const {
    return MultiMaterialEvolutionProblemBase::getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblemImplementation<true>::getMaterial(
      const size_type m) {
    return MultiMaterialEvolutionProblemBase::getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator& NonLinearEvolutionProblemImplementation<
      true>::getBehaviourIntegrator(const size_type m) const {
    return MultiMaterialEvolutionProblemBase::getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator& NonLinearEvolutionProblemImplementation<
      true>::getBehaviourIntegrator(const size_type m) {
    return MultiMaterialEvolutionProblemBase::getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblemImplementation<true>::setup(const real t,
                                                            const real dt) {
    NonLinearEvolutionProblemImplementationBase::setup(t, dt);
    MultiMaterialEvolutionProblemBase::setup(t, dt);
  }  // end of setup

  void NonLinearEvolutionProblemImplementation<true>::revert() {
    NonLinearEvolutionProblemImplementationBase::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblemImplementation<true>::update() {
    NonLinearEvolutionProblemImplementationBase::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblemImplementation<true>::setTimeIncrement(
      const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

  const FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

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
          const Parameters&)
      : NonLinearEvolutionProblemImplementationBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h),
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
    this->AddDomainIntegrator(this->mgis_integrator);
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

  void NonLinearEvolutionProblemImplementation<false>::executePostProcessings(
      const real t, const real dt) {
    for (auto& p : this->postprocessings) {
      p->execute(*this, t, dt);
    }
  }  // end of executePostProcessings

  void NonLinearEvolutionProblemImplementation<false>::addBehaviourIntegrator(
      const std::string& n,
      const size_type l,
      const std::string& m,
      const std::string& b) {
    MultiMaterialEvolutionProblemBase::addBehaviourIntegrator(n, l, m, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblemImplementation<false>::getMaterial(
      const size_type m) const {
    return MultiMaterialEvolutionProblemBase::getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblemImplementation<false>::getMaterial(
      const size_type m) {
    return MultiMaterialEvolutionProblemBase::getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator& NonLinearEvolutionProblemImplementation<
      false>::getBehaviourIntegrator(const size_type m) const {
    return MultiMaterialEvolutionProblemBase::getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator& NonLinearEvolutionProblemImplementation<
      false>::getBehaviourIntegrator(const size_type m) {
    return MultiMaterialEvolutionProblemBase::getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblemImplementation<false>::setup(const real t,
                                                             const real dt) {
    NonLinearEvolutionProblemImplementationBase::setup(t, dt);
    MultiMaterialEvolutionProblemBase::setup(t, dt);
  }  // end of setup

  void NonLinearEvolutionProblemImplementation<false>::revert() {
    NonLinearEvolutionProblemImplementationBase::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblemImplementation<false>::update() {
    NonLinearEvolutionProblemImplementationBase::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblemImplementation<false>::setTimeIncrement(
      const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

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
