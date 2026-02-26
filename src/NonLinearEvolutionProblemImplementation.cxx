/*!
 * \file   src/NonLinearEvolutionProblemImplementation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "mfem/fem/datacollection.hpp"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif MFEM_USE_PETSC
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pbilinearform.hpp"
#endif /* MFEM_USE_MPI */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/PostProcessingFactory.hxx"
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  struct PredictionResult {
    //! \brief prediction of the increment of the unknowns
    std::unique_ptr<mfem_mgis::GridFunction<parallel>> du;
    //! \brief initial residual, if available
    const real initial_residual_norm;
  };

  template <bool parallel>
  [[nodiscard]] std::optional<PredictionResult<parallel>> computePrediction(
      mfem_mgis::Context& ctx,
      mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& p,
      LinearSolver& ls,
      const real t,
      const real dt) noexcept {
    auto& fed = p.getFiniteElementDiscretization();
    auto& fespace = fed.template getFiniteElementSpace<parallel>();
    const auto& u0 = p.getUnknowns(mfem_mgis::bts);
    auto success = p.integrate(
        u0, mfem_mgis::IntegrationType::PREDICTION_ELASTIC_OPERATOR);
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LAND,
                    fespace.GetComm());
    }
#endif /* MFEM_USE_MPI */
    if (!success) {
      return ctx.registerErrorMessage("integration failure");
    }
    auto operators = p.getLinearizedOperators(ctx, u0);
    if (isInvalid(operators)) {
      return {};
    }
    //
    auto du = std::make_unique<mfem_mgis::GridFunction<parallel>>(&fespace);
    *du = 0.0;
    auto a = mfem_mgis::BilinearForm<parallel>(&fespace);
    a.AddDomainIntegrator(operators->K.release());
    auto b = mfem_mgis::LinearForm<parallel>(&fespace);
    if (operators->mFi.get() != nullptr) {
      b.AddDomainIntegrator(operators->mFi.release());
    }
    //
    for (const auto& bc : p.getBoundaryConditions()) {
      if (!bc->addLinearFormIntegrators(ctx, a, b, u0, t, dt)) {
        return {};
      }
    }
    //
    a.Assemble();
    b.Assemble();
    //
    auto A = [] {
#ifdef MFEM_USE_MPI
      if constexpr (parallel) {
        return mfem::HypreParMatrix{};
      } else {
        return mfem::SparseMatrix{};
      }
#else
      return mfem::SparseMatrix{};
#endif /* MFEM_USE_MPI */
    }();
    mfem::Vector B, X;
    //
    auto essential_dofs = p.getEssentialDegreesOfFreedom();
    auto edofs_list = mfem::Array<mfem_mgis::size_type>(essential_dofs.data(),
                                                        essential_dofs.size());
    if constexpr (parallel) {
      auto du_values = mfem::Vector(fespace.GetTrueVSize());
      du_values = 0.0;
      //
      for (const auto& bc : p.getDirichletBoundaryConditions()) {
        bc->setImposedValuesIncrements(du_values, t, t + dt);
      }
      du->Distribute(du_values);
    } else {
      for (const auto& bc : p.getDirichletBoundaryConditions()) {
        bc->setImposedValuesIncrements(*du, t, t + dt);
      }
    }
    //
    a.FormLinearSystem(edofs_list, *du, b, A, X, B);
    const auto norm = [&B, &fespace] {
      if constexpr (parallel) {
        return std::sqrt(mfem::InnerProduct(fespace.GetComm(), B, B));
      } else {
        return std::sqrt(mfem::InnerProduct(B, B));
      }
    }();
    ls.SetOperator(A);
    ls.Mult(B, X);
    // check for convergence
    const auto usesIterativeLinearSolver =
        dynamic_cast<const IterativeSolver*>(&ls) != nullptr;
    if (usesIterativeLinearSolver) {
      const auto& isolver = static_cast<const mfem::IterativeSolver&>(ls);
      if (!isolver.GetConverged()) {
        return ctx.registerErrorMessage("linear solver did not converge");
      }
    }
    //
    a.RecoverFEMSolution(X, b, *du);
    //
    return PredictionResult<parallel>{.du = std::move(du),
                                      .initial_residual_norm = norm};
  }  // end of computePrediction

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
        unknowns0(&(fed->getFiniteElementSpace<true>()), this->u0),
        unknowns1(&(fed->getFiniteElementSpace<true>()), this->u1) {
    if (this->fe_discretization->getMesh<true>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
#ifdef MFEM_USE_PETSC
    if (usePETSc()) {
      this->petsc_solver = std::make_unique<mfem::PetscNonlinearSolver>(
          this->getFiniteElementSpace().GetComm(), *this);
      this->petsc_solver->iterative_mode = true;
    } else {
      this->solver = std::make_unique<NewtonSolver>(*this);
    }
#else  /* MFEM_USE_PETSC */
    this->solver = std::make_unique<NewtonSolver>(*this);
#endif /* MFEM_USE_PETSC */
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  template <bool parallel>
  void export_prediction(
      mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& p,
      mfem_mgis::GridFunction<parallel>& du) {
    auto& fed = p.getFiniteElementDiscretization();
    auto exporter = mfem::ParaViewDataCollection{
        parallel ? "result-prediction-test-parallel"
                 : "result-prediction-test"};
    exporter.SetMesh(&(fed.template getMesh<parallel>()));
    exporter.SetDataFormat(mfem::VTKFormat::BINARY);
    exporter.RegisterField("DisplacementPrediction", &du);
    exporter.SetCycle(1);
    exporter.SetTime(1);
    exporter.Save();
  }

  std::optional<real>
  NonLinearEvolutionProblemImplementation<true>::computePrediction(
      Context& ctx, const real t, const real dt) noexcept {
    if (this->linear_solver.get() == nullptr) {
      return ctx.registerErrorMessage("linear solver not initialized");
    }
    const auto oresult = ::mfem_mgis::computePrediction<true>(
        ctx, *this, *(this->linear_solver), t, dt);
#pragma message("required ?")
    this->solver->setLinearSolver(*(this->linear_solver));
    if (isInvalid(oresult)) {
      return {};
    }
    this->u1 = this->u0;
    this->u1 += oresult->du->GetTrueVector();
    return oresult->initial_residual_norm;
  }  // end of computePrediction

  void NonLinearEvolutionProblemImplementation<true>::Mult(
      const mfem::Vector& u, mfem::Vector& r) const {
    return mfem_mgis::NonlinearForm<true>::Mult(u, r);
  }  // end of Mult

  void NonLinearEvolutionProblemImplementation<true>::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of addBoundaryCondition

  void NonLinearEvolutionProblemImplementation<true>::addBoundaryCondition(
      std::unique_ptr<AbstractBoundaryCondition> f) {
    auto ctx = Context{};
    if (!this->addBoundaryCondition(ctx, std::move(f))) {
      raise(ctx.getErrorMessage());
    }
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<true>::addBoundaryCondition(
      Context& ctx, std::unique_ptr<AbstractBoundaryCondition> f) noexcept {
    if (!f->addNonlinearFormIntegrator(ctx, *this, this->u1)) {
      return false;
    }
    this->boundary_conditions.push_back(std::move(f));
    return true;
  }  // end of addBoundaryCondition

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      std::unique_ptr<PostProcessing<true>> p) {
    this->postprocessings.push_back(std::move(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblemImplementation<true>::addPostProcessing(
      std::string_view n, const Parameters& p) {
    const auto& f = PostProcessingFactory<true>::getFactory();
    this->postprocessings.push_back(f.generate(n, *this, p));
  }  // end of addPostProcessing

  bool NonLinearEvolutionProblemImplementation<true>::setLinearSolver(
      Context& ctx, LinearSolverHandler s) noexcept {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid linear solver");
    }
    this->updateLinearSolver(std::move(s));
    return true;
  }

  void NonLinearEvolutionProblemImplementation<true>::setLinearSolver(
      std::string_view n, const Parameters& p) {
    auto ctx = Context{};
    const auto& f = LinearSolverFactory<true>::getFactory();
    auto& fespace = this->getFiniteElementSpace();
    auto s = f.generate(ctx, n, fespace, p);
    if (isInvalid(s)) {
      raise(ctx.getErrorMessage());
    }
    this->updateLinearSolver(std::move(s));
  }  // end of setLinearSolver

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

  Mesh<true>& NonLinearEvolutionProblemImplementation<true>::getMesh() {
    return this->fe_discretization->getMesh<true>();
  }  // end of getMesh

  const Mesh<true>& NonLinearEvolutionProblemImplementation<true>::getMesh()
      const {
    return this->fe_discretization->getMesh<true>();
  }  // end of getMesh

  FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  const FiniteElementSpace<true>&
  NonLinearEvolutionProblemImplementation<true>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  bool NonLinearEvolutionProblemImplementation<true>::integrate(
      const mfem::Vector& u, const IntegrationType it) {
    if (this->mgis_integrator == nullptr) {
      return true;
    }
    const auto& pu = this->Prolongate(u);
    const auto& fespace = this->getFiniteElementSpace();
    mfem::Array<int> vdofs;
    mfem::Vector ue;
    bool noerror = true;
    for (size_type i = 0; noerror && (i != fespace.GetNE()); ++i) {
      const auto& e = *(fespace.GetFE(i));
      auto& tr = *(fespace.GetElementTransformation(i));
      fespace.GetElementVDofs(i, vdofs);
      pu.GetSubVector(vdofs, ue);
      noerror = this->mgis_integrator->integrate(e, tr, ue, it);
    }
    MPI_Allreduce(MPI_IN_PLACE, &noerror, 1, MPI_C_BOOL, MPI_LAND,
                  MPI_COMM_WORLD);
    return noerror;
  }  // end of integrate

  void NonLinearEvolutionProblemImplementation<true>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  GridFunction<true>& NonLinearEvolutionProblemImplementation<
      true>::getUnknownsAsGridFunction(const TimeStepStage ts) noexcept {
    if (ts == bts) {
      return this->unknowns0;
    }
    return this->unknowns1;
  }  // end of getUnknownsAsGridFunction

  const GridFunction<true>& NonLinearEvolutionProblemImplementation<
      true>::getUnknownsAsGridFunction(const TimeStepStage ts) const noexcept {
    if (ts == bts) {
      return this->unknowns0;
    }
    return this->unknowns1;
  }  // end of getUnknownsAsGridFunction

  NonLinearEvolutionProblemImplementation<
      true>::~NonLinearEvolutionProblemImplementation() = default;

#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblemImplementation<false>::
      NonLinearEvolutionProblemImplementation(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : NonLinearEvolutionProblemImplementationBase(fed, h, p),
        mfem::NonlinearForm(&(fed->getFiniteElementSpace<false>())),
        unknowns0(&(fed->getFiniteElementSpace<false>()), this->u0),
        unknowns1(&(fed->getFiniteElementSpace<false>()), this->u1) {
    this->solver = std::make_unique<NewtonSolver>(*this);
    if (this->fe_discretization->getMesh<false>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      raise(
          "NonLinearEvolutionProblemImplementationBase::"
          "NonLinearEvolutionProblemImplementationBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    if (this->mgis_integrator != nullptr) {
      this->AddDomainIntegrator(this->mgis_integrator);
    }
  }  // end of NonLinearEvolutionProblemImplementation

  std::optional<real>
  NonLinearEvolutionProblemImplementation<false>::computePrediction(
      Context& ctx, const real t, const real dt) noexcept {
    if (this->linear_solver.get() == nullptr) {
      return ctx.registerErrorMessage("linear solver not initialized");
    }
    const auto oresult = ::mfem_mgis::computePrediction<false>(
        ctx, *this, *(this->linear_solver), t, dt);
    this->solver->setLinearSolver(*(this->linear_solver));
    if (isInvalid(oresult)) {
      return {};
    }
    this->u1 = this->u0;
    this->u1 += *(oresult->du);
    return oresult->initial_residual_norm;
  }  // end of computePrediction

  void NonLinearEvolutionProblemImplementation<false>::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of addBoundaryCondition

  void NonLinearEvolutionProblemImplementation<false>::addBoundaryCondition(
      std::unique_ptr<AbstractBoundaryCondition> f) {
    auto ctx = Context{};
    if (!this->addBoundaryCondition(ctx, std::move(f))) {
      raise(ctx.getErrorMessage());
    }
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementation<false>::addBoundaryCondition(
      Context& ctx, std::unique_ptr<AbstractBoundaryCondition> f) noexcept {
    if (!f->addNonlinearFormIntegrator(ctx, *this, this->u1)) {
      return false;
    }
    this->boundary_conditions.push_back(std::move(f));
    return true;
  }  // end of addBoundaryCondition

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

  Mesh<false>& NonLinearEvolutionProblemImplementation<false>::getMesh() {
    return this->fe_discretization->getMesh<false>();
  }  // end of getMesh

  const Mesh<false>& NonLinearEvolutionProblemImplementation<false>::getMesh()
      const {
    return this->fe_discretization->getMesh<false>();
  }  // end of getMesh

  const FiniteElementSpace<false>& NonLinearEvolutionProblemImplementation<
      false>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<false>&
  NonLinearEvolutionProblemImplementation<false>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  bool NonLinearEvolutionProblemImplementation<false>::setLinearSolver(
      Context& ctx, LinearSolverHandler s) noexcept {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid linear solver");
    }
    this->updateLinearSolver(std::move(s));
    return true;
  }

  void NonLinearEvolutionProblemImplementation<false>::setLinearSolver(
      std::string_view n, const Parameters& p) {
    auto ctx = Context{};
    const auto& f = LinearSolverFactory<false>::getFactory();
    auto& fespace = this->getFiniteElementSpace();
    auto s = f.generate(ctx, n, fespace, p);
    if (isInvalid(s)) {
      raise(ctx.getErrorMessage());
    }
    this->updateLinearSolver(std::move(s));
  }  // end of setLinearSolver

  bool NonLinearEvolutionProblemImplementation<false>::integrate(
      const mfem::Vector& u, const IntegrationType it) {
    if (this->mgis_integrator == nullptr) {
      return true;
    }
    const auto& pu = this->Prolongate(u);
    const auto& fespace = this->getFiniteElementSpace();
    mfem::Array<int> vdofs;
    mfem::Vector ue;
    for (size_type i = 0; i != fespace.GetNE(); ++i) {
      const auto& e = *(fespace.GetFE(i));
      auto& tr = *(fespace.GetElementTransformation(i));
      fespace.GetElementVDofs(i, vdofs);
      pu.GetSubVector(vdofs, ue);
      if (!this->mgis_integrator->integrate(e, tr, ue, it)) {
        return false;
      }
    }
    return true;
  }  // end of integrate

  void NonLinearEvolutionProblemImplementation<false>::
      markDegreesOfFreedomHandledByDirichletBoundaryConditions(
          std::vector<size_type> dofs) {
    auto tmp = mfem::Array<size_type>(dofs.data(), dofs.size());
    this->SetEssentialTrueDofs(tmp);
  }  // end of markDegreesOfFreedomHandledByDirichletBoundaryConditions

  GridFunction<false>& NonLinearEvolutionProblemImplementation<
      false>::getUnknownsAsGridFunction(const TimeStepStage ts) noexcept {
    if (ts == bts) {
      return this->unknowns0;
    }
    return this->unknowns1;
  }  // end of getUnknownsAsGridFunction

  const GridFunction<false>& NonLinearEvolutionProblemImplementation<
      false>::getUnknownsAsGridFunction(const TimeStepStage ts) const noexcept {
    if (ts == bts) {
      return this->unknowns0;
    }
    return this->unknowns1;
  }  // end of getUnknownsAsGridFunction

  NonLinearEvolutionProblemImplementation<
      false>::~NonLinearEvolutionProblemImplementation() = default;

}  // end of namespace mfem_mgis
