/*!
 * \file   L2Projection.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/01/2026
 */

#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/mesh/submesh/submesh.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/mesh/submesh/psubmesh.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pbilinearform.hpp"
#endif

#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/L2Projection.hxx"

namespace mfem_mgis {

  template <bool parallel>
  [[nodiscard]] static bool checkConsistency_impl(
      Context& ctx,
      const FiniteElementDiscretization& fed0,
      const FiniteElementDiscretization& fed1) noexcept {
    //
    const auto& mesh0 = fed0.getMesh<parallel>();
    const auto& mesh1 = fed1.getMesh<parallel>();
    if (&mesh0 != &mesh1) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: "
          "underling meshes are distinct");
    }
    //
    const auto& fes0 = fed0.getFiniteElementSpace<parallel>();
    const auto& fes1 = fed1.getFiniteElementSpace<parallel>();
    if ((fes0.FEColl() != fes1.FEColl()) ||
        (fes0.GetOrdering() != fes1.GetOrdering())) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: "
          "underling finite element spaces are distinct");
    }
    return true;
  }  // end of checkConsistency_impl

  [[nodiscard]] static bool checkConsistency(
      Context& ctx,
      const FiniteElementDiscretization& fed0,
      const FiniteElementDiscretization& fed1) {
    const auto parallel = fed0.describesAParallelComputation();
    if (parallel) {
      if (!fed1.describesAParallelComputation()) {
        return ctx.registerErrorMessage(
            "inconsistent finite element discretization space: one is "
            "parallel, while the other is not");
      }
#ifdef MFEM_USE_MPI
      return checkConsistency_impl<true>(ctx, fed0, fed1);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    if (fed1.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "inconsistent finite element discretization space: one is "
          "parallel, while the other is not");
    }
    return checkConsistency_impl<false>(ctx, fed0, fed1);
  }  // end of checkConsistency

  template <bool parallel>
  [[nodiscard]] static bool checkConsistency(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return true;
    }
    const auto& fed0 =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed0.getMesh<parallel>();
    const auto& elts_attributes = mesh.GetElementAttributes();
    auto ids = std::vector<size_type>();
    ids.reserve(fcts.size());
    for (const auto& f : fcts) {
      const auto& qspace = f.getPartialQuadratureSpace();
      const auto& fed = qspace.getFiniteElementDiscretization();
      if (!checkConsistency(ctx, fed0, fed)) {
        return false;
      }
      const auto id = qspace.getId();
      if (elts_attributes.Find(id) == -1) {
        return ctx.registerErrorMessage("material id '" + std::to_string(id) +
                                        "' is not a element attribute");
      }
      if (std::find(ids.begin(), ids.end(), id) != ids.end()) {
        return ctx.registerErrorMessage(
            "several quadrature functions associated with material id '" +
            std::to_string(id) + "' given");
      }
      ids.push_back(id);
    }
    return true;
  }  // end of checkConsistency

  template <bool parallel>
  [[nodiscard]] static std::optional<bool> areDefinedOnWholeMesh_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    //
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    // thanks to checkConsistency we know that the element attributes contains
    // all the material ids, so we just need to check that the number of
    // attributes equals to number of functions
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    const auto& elts_attributes = mesh.GetElementAttributes();
    return elts_attributes.Size() == static_cast<size_type>(fcts.size());
  }  // end of areDefinedOnWholeMesh_impl

  [[nodiscard]] static std::optional<bool> areDefinedOnWholeMesh(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return areDefinedOnWholeMesh_impl<true>(ctx, fcts);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return areDefinedOnWholeMesh_impl<false>(ctx, fcts);
  }

  template <bool parallel>
  [[nodiscard]] static std::optional<L2ProjectionResult<parallel>>
  createL2ProjectionResult_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    const auto onAll = areDefinedOnWholeMesh(ctx, fcts);
    if (isInvalid(onAll)) {
      return {};
    }
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    auto r = L2ProjectionResult<parallel>{};
    if (*onAll) {
      auto oresult = makeGridFunction<parallel>(ctx, fcts, mesh);
      if (isInvalid(oresult)) {
        return {};
      }
      r.fe_space.swap(oresult->first);
      r.result.swap(oresult->second);
      return r;
    }
    // create sub mesh
    auto ids = mfem::Array<int>{};
    for (const auto& f : fcts) {
      const auto& qspace = f.getPartialQuadratureSpace();
      const auto id = qspace.getId();
      ids.Append(id);
    }
    r.submesh = std::make_unique<SubMesh<parallel>>(
        SubMesh<parallel>::CreateFromDomain(mesh, ids));
    //
    auto oresult = makeGridFunction<parallel>(ctx, fcts, *(r.submesh));
    if (isInvalid(oresult)) {
      return {};
    }
    r.fe_space.swap(oresult->first);
    r.result.swap(oresult->second);
    return r;
  }  // end of createL2ProjectionResult_impl

  template <>
  std::optional<L2ProjectionResult<true>> createL2ProjectionResult<true>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return createL2ProjectionResult_impl<true>(ctx, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of createL2ProjectionResult

  template <>
  std::optional<L2ProjectionResult<false>> createL2ProjectionResult<false>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return createL2ProjectionResult_impl<false>(ctx, fcts);
  }  // end of createL2ProjectionResult

  struct L2ProjectionFormIntegratorBase {
    L2ProjectionFormIntegratorBase(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : functions(makeFunctionsMapping(throwing, fcts)) {}

   private:
    [[nodiscard]] static std::
        unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        makeFunctionsMapping(
            attributes::Throwing,
            const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
      auto functions =
          std::unordered_map<size_type,
                             ImmutablePartialQuadratureFunctionView>{};
      for (const auto& f : fcts) {
        const auto mid = f.getPartialQuadratureSpace().getId();
        if (!functions.insert({mid, f}).second) {
          raise("multiple functions defined for material '" +
                std::to_string(mid) + "'");
        }
      }
      return functions;
    }  // end of makeFunctionsMapping

   protected:
    // functions
    std::unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        functions;
    mfem::Vector shape;
  };

  struct ScalarL2ProjectionRHSFormIntegrator final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    using L2ProjectionFormIntegratorBase::L2ProjectionFormIntegratorBase;
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto v =
            tr.Weight() * f.getIntegrationPointValue(tr.ElementNo, i);
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect
  };

  struct ScalarL2ProjectionRHSFormIntegratorII final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ScalarL2ProjectionRHSFormIntegratorII(
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : L2ProjectionFormIntegratorBase(fcts), elts_mapping(m) {}

    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto n = elts_mapping[tr.ElementNo];
        const auto v = tr.Weight() * f.getIntegrationPointValue(n, i);
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect

   private:
    const mfem::Array<int> elts_mapping;
  };

  struct ComponentL2ProjectionRHSFormIntegrator final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ComponentL2ProjectionRHSFormIntegrator(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
        const size_type c)
        : L2ProjectionFormIntegratorBase(fcts), component(c) {}
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto fv =
            f.getIntegrationPointValues(tr.ElementNo, i)[this->component];
        const auto v = tr.Weight() * fv;
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect
   private:
    const size_type component;
  };

  struct ComponentL2ProjectionRHSFormIntegratorII final
      : public L2ProjectionFormIntegratorBase,
        public mfem::LinearFormIntegrator {
    //
    ComponentL2ProjectionRHSFormIntegratorII(
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
        const size_type c)
        : L2ProjectionFormIntegratorBase(fcts), elts_mapping(m), component(c) {}

    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& elvect) override {
      const auto m = tr.Attribute;
      const auto& f = this->functions.at(m);
      const int dof = e.GetDof();
      //
      shape.SetSize(dof);  // vector of size dof
      elvect.SetSize(dof);
      elvect = 0.0;
      //
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(e, tr);
      //
      for (int i = 0; i < ir.GetNPoints(); i++) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        const auto n = elts_mapping[tr.ElementNo];
        const auto fv = f.getIntegrationPointValues(n, i)[this->component];
        const auto v = tr.Weight() * fv;
        e.CalcPhysShape(tr, shape);
        add(elvect, ip.weight * v, shape, elvect);
      }
    }  // end of AssembleRHSElementVect

   private:
    const mfem::Array<int> elts_mapping;
    const size_type component;
  };

  struct VectorL2ProjectionBilinearFormIntegrator final
      : public mfem::BilinearFormIntegrator {
    //
    VectorL2ProjectionBilinearFormIntegrator(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : nc(fcts.at(0).getNumberOfComponents()) {}
    //
    void AssembleElementMatrix(const mfem::FiniteElement& el,
                               mfem::ElementTransformation& Trans,
                               mfem::DenseMatrix& elmat) override {
      const int nd = el.GetDof();
      elmat.SetSize(nd * (this->nc));
      shape.SetSize(nd);
      partelmat.SetSize(nd);
      //
      const auto* ir = GetIntegrationRule(el, Trans);
      if (ir == NULL) {
        const auto order = 2 * el.GetOrder() + Trans.OrderW();
        ir = &mfem::IntRules.Get(el.GetGeomType(), order);
      }
      //
      elmat = 0.0;
      for (size_type s = 0; s < ir->GetNPoints(); s++) {
        const auto& ip = ir->IntPoint(s);
        Trans.SetIntPoint(&ip);
        el.CalcPhysShape(Trans, shape);
        const auto norm = ip.weight * Trans.Weight();
        MultVVt(shape, partelmat);
        partelmat *= norm;
        for (size_type k = 0; k < this->nc; k++) {
          elmat.AddMatrix(partelmat, nd * k, nd * k);
        }
      }
    }  // end of AssembleElementMatrix
   private:
    mfem::Vector shape;
    mfem::DenseMatrix partelmat;
    //! \brief number of components
    const int nc;
  };

  template <bool parallel>
  static bool updateL2Projection_impl(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    // checks
    if ((r.fe_space == nullptr) || (r.result == nullptr)) {
      return ctx.registerErrorMessage("uninitialized result");
    }
    if (r.result->FESpace() != r.fe_space.get()) {
      return ctx.registerErrorMessage("inconsistent finite element space");
    }
    if (l.linear_solver == nullptr) {
      return ctx.registerErrorMessage("uninitialized linear solver");
    }
    // check consistency
    if (fcts.empty()) {
      return ctx.registerErrorMessage(
          "empty list of partial quadrature functions given");
    }
    if (!checkConsistency<parallel>(ctx, fcts)) {
      return {};
    }
    if (r.fe_space->GetVDim() != getNumberOfComponents(fcts.front())) {
      return ctx.registerErrorMessage("inconsistent number of components");
    }
    //
    const auto& mesh = *(r.fe_space->GetMesh());
    const auto& elts_attributes = mesh.GetElementAttributes();
    for (int i = 0; i != elts_attributes.Size(); ++i) {
      const auto id = elts_attributes[i];
      const auto found = [&id, &fcts] {
        for (const auto& f : fcts) {
          const auto& qspace = f.getPartialQuadratureSpace();
          if (id == qspace.getId()) {
            return true;
          }
        }
        return false;
      }();
      if (!found) {
        return ctx.registerErrorMessage(
            "no partial quadrature function defined on material id '" +
            std::to_string(id) + "'");
      }
    }
    //
    *(r.result) = real{};
    //
    // Here begin the tricky part. In order to reduce the computational cost,
    // we will make the projection component by component. But of course, if the
    // functions are scalar, we do want to avoid creating a temporary finite
    // element space or a copy of the projection in *(r.result).
    //
    // Thus we introduce:
    // -  local_fespace and local_gridfunction which are only allocated
    // -  fespace and x which points with the finite element space and grid
    //    function used for the resolution
    //
    auto local_fespace = std::unique_ptr<FiniteElementSpace<parallel>>{};
    auto* const fespace = [&r, &local_fespace, &mesh] {
      if (r.fe_space->GetVDim() == 1) {
        return r.fe_space.get();
      }
      if constexpr (parallel) {
        local_fespace = std::make_unique<FiniteElementSpace<parallel>>(
            r.fe_space->GetParMesh(), r.fe_space->FEColl(), 1,
            r.fe_space->GetOrdering());
      } else {
        local_fespace = std::make_unique<FiniteElementSpace<parallel>>(
            r.fe_space->GetMesh(), r.fe_space->FEColl(), 1,
            r.fe_space->GetOrdering());
      }
      return local_fespace.get();
    }();
    //
    auto local_gridfunction = std::unique_ptr<GridFunction<parallel>>{};
    auto* const x = [&r, &local_gridfunction, &fespace] {
      if (r.fe_space->GetVDim() == 1) {
        return r.result.get();
      }
      local_gridfunction = std::make_unique<GridFunction<parallel>>(fespace);
      return local_gridfunction.get();
    }();
    // resolutions
    // Mass matrix
    BilinearForm<parallel> a(fespace);
    a.AddDomainIntegrator(new mfem::MassIntegrator);
    a.Assemble();
    for (size_type c = 0; c != r.fe_space->GetVDim(); ++c) {
      // initialize the solution, if need
      if (r.fe_space->GetVDim() != 1) {
        *x = real{};
      }
      // right-hand side
      LinearForm<parallel> b(fespace);
      if (r.fe_space->GetVDim() == 1) {
        if (r.submesh != nullptr) {
          b.AddDomainIntegrator(new ScalarL2ProjectionRHSFormIntegratorII(
              r.submesh->GetParentElementIDMap(), fcts));
        } else {
          b.AddDomainIntegrator(new ScalarL2ProjectionRHSFormIntegrator(fcts));
        }
      } else {
        if (r.submesh != nullptr) {
          b.AddDomainIntegrator(new ComponentL2ProjectionRHSFormIntegratorII(
              r.submesh->GetParentElementIDMap(), fcts, c));
        } else {
          b.AddDomainIntegrator(
              new ComponentL2ProjectionRHSFormIntegrator(fcts, c));
        }
      }
      b.Assemble();
      // resolution
      mfem::Array<int> boundary_dofs;
      std::conditional_t<parallel, mfem::HypreParMatrix, mfem::SparseMatrix> A;
      mfem::Vector B, X;
      a.FormLinearSystem(boundary_dofs, *x, b, A, X, B);
      l.linear_solver->SetOperator(A);
      l.linear_solver->Mult(B, X);
      auto* const isolver =
          dynamic_cast<IterativeSolver*>(l.linear_solver.get());
      if (isolver != nullptr) {
        if (!isolver->GetConverged()) {
          return ctx.registerErrorMessage("linear solver failed");
        }
      }
      a.RecoverFEMSolution(X, b, *x);
      // copy solution to r.result, if needed
      if (r.fe_space->GetVDim() != 1) {
        auto& dest = *(r.result);
        const auto& src = *(x);
        const auto bynodes =
            (r.fe_space->GetOrdering() == mfem::Ordering::byNODES);
        if (bynodes) {
          const auto n = x->Size();
          const auto offset = x->Size() * c;
          for (size_type idx = 0; idx != n; ++idx) {
            dest[offset + idx] = src[idx];
          }
        } else {
          const auto n = r.fe_space->GetVDim();
          for (size_type idx = 0; idx != x->Size(); ++idx) {
            dest[idx * n + c] = src[idx];
          }
        }
      }
    }
    return true;
  }  // end of updateL2Projection_impl

  template <>
  bool updateL2Projection<true>(
      Context& ctx,
      L2ProjectionResult<true>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return updateL2Projection_impl<true>(ctx, r, l, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of updateL2Projection

  template <>
  bool updateL2Projection<false>(
      Context& ctx,
      L2ProjectionResult<false>& r,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return updateL2Projection_impl<false>(ctx, r, l, fcts);
  }  // end of updateL2Projection

  template <bool parallel>
  static std::optional<L2ProjectionResult<parallel>> computeL2Projection_impl(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    auto ores = createL2ProjectionResult<parallel>(ctx, fcts);
    if (isInvalid(ores)) {
      return {};
    }
    const auto ok = updateL2Projection<parallel>(ctx, *ores, l, fcts);
    if (!ok) {
      return {};
    }
    return ores;
  }  // end of computeL2Projection_impl

  template <>
  std::optional<L2ProjectionResult<true>> computeL2Projection<true>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
#ifdef MFEM_USE_MPI
    return computeL2Projection_impl<true>(ctx, l, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of computeL2Projection

  template <>
  std::optional<L2ProjectionResult<false>> computeL2Projection<false>(
      Context& ctx,
      LinearSolverHandler& l,
      const std::vector<ImmutablePartialQuadratureFunctionView>&
          fcts) noexcept {
    return computeL2Projection_impl<false>(ctx, l, fcts);
  }  // end of computeL2Projection

}  // end of namespace mfem_mgis