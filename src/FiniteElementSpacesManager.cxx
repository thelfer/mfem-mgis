/*!
 * \file   src/FiniteElementSpacesManager.cxx
 * \brief  this file implements the `FiniteElementSpacesManager` class
 * \author Thomas Helfer
 * \date   09/09/2026
 */

#include <mfem/mesh/mesh.hpp>
#include <mfem/fem/fespace.hpp>
#include <mfem/fem/gridfunc.hpp>
#include "mfem/mesh/submesh/submesh.hpp"
#ifdef MFEM_USE_MPI
#include <mfem/mesh/pmesh.hpp>
#include <mfem/fem/pfespace.hpp>
#include "mfem/mesh/submesh/psubmesh.hpp"
#endif
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"

namespace mfem_mgis {

  const char* const FiniteElementSpacesManager::FiniteElementFamily =
      "FiniteElementFamily";
  const char* const FiniteElementSpacesManager::FiniteElementOrder =
      "FiniteElementOrder";

  std::vector<std::string>
  FiniteElementSpacesManager::getFiniteElementCollectionParametersList() {
    return {FiniteElementSpacesManager::FiniteElementFamily,
            FiniteElementSpacesManager::FiniteElementOrder};
  }  // end of getFiniteElementCollectionParametersList

  std::vector<std::string> FiniteElementSpacesManager::getParametersList() {
    auto d = MeshDiscretization::getParametersList();
    const auto names =
        FiniteElementSpacesManager::getFiniteElementCollectionParametersList();
    d.insert(d.end(), names.begin(), names.end());
    return d;
  }  // end of getParametersList

  [[nodiscard]] static std::shared_ptr<const FiniteElementCollection>
  buildFiniteElementCollection(attributes::Throwing,
                               MeshDiscretization& m,
                               const Parameters& params) {
    checkParameters(
        throwing, params,
        FiniteElementSpacesManager::getFiniteElementCollectionParametersList());
    const auto& fe_family = get_if<std::string>(
        throwing, params, FiniteElementSpacesManager::FiniteElementFamily,
        "H1");
    const auto fe_order = get_if<int>(
        throwing, params, FiniteElementSpacesManager::FiniteElementOrder, 1);
    // building the finite element collection
    if (fe_family != "H1") {
      raise(
          "FiniteElementSpacesManager::FiniteElementSpacesManager: "
          "unsupported finite element family '" +
          fe_family + "'");
    }
    return std::make_shared<mfem::H1_FECollection>(fe_order,
                                                   getSpaceDimension(m));
  }  // end of buildFiniteElementCollectionAndSpace

  struct FiniteElementSpacesManager::Implementation {
    /*!
     * \brief constructor from parameters
     * \param[in] ctx: execution context
     * \param[in] parameters: parameters
     */
    Implementation(Context& ctx, const Parameters& parameters)
        : mesh(ctx,
               extract(throwing,
                       parameters,
                       MeshDiscretization::getParametersList())) {
      CatchTimeSection(ctx, "FiniteElementSpacesManager::Constructor");
      checkParameters(throwing, parameters,
                      FiniteElementSpacesManager::getParametersList());
      this->fec = buildFiniteElementCollection(
          throwing, this->mesh,
          remove(parameters, MeshDiscretization::getParametersList()));
    }  // end of Implementation
    /*!
     * \brief constructor from a mesh discretization
     * \param[in] m: mesh discretization
     * \param[in] parameters: parameters
     */
    Implementation(const MeshDiscretization& m, const Parameters& parameters)
        : mesh(m) {
      checkParameters(throwing, parameters,
                      FiniteElementSpacesManager::
                          getFiniteElementCollectionParametersList());
      this->fec =
          buildFiniteElementCollection(throwing, this->mesh, parameters);
    }  // end of Implementation
    /*!
     * \brief constructor from a parallel mesh
     * \param[in] m: mesh
     * \param[in] c: finite element collection
     */
    Implementation(const MeshDiscretization& m,
                   std::shared_ptr<const FiniteElementCollection> c)
        : mesh(m), fec(c) {
      if (c.get() == nullptr) {
        raise("invalid finite element collection pointer");
      }
    }  // end of Implementation
    //! \return the mesh discretization
    MeshDiscretization getMeshDiscretization() const noexcept {
      return this->mesh;
    }
    /*!
     * \brief create a new finite element space or reuse an existing one
     * \param[in] ctx: execution context
     * \param[in] nc: vectorial dimension
     *
     * \note if a finite element space is created, it is stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context& ctx, const size_type nc) noexcept {
      if (nc < 1) {
        return ctx.registerErrorMessage(
            "invalid number of components (vectorial dimension), (" +
            std::to_string(nc) + ") given");
      }
      if constexpr (parallel) {
#ifdef MFEM_USE_MPI
        if (!this->mesh.describesAParallelComputation()) {
          return ctx.registerErrorMessage(
              "can't create a parallel finite element space on top of a "
              "sequential mesh");
        }
        auto p = this->parallel_fespaces.find(nc);
        if (p != this->parallel_fespaces.end()) {
          return p->second;
        }
        auto ptr = make_shared<FiniteElementSpace<true>>(
            ctx, this->mesh.getMutableMeshPointer<true>().get(),
            this->fec.get(), nc);
        this->parallel_fespaces.insert({nc, ptr});
        return ptr;
#else
        reportUnsupportedParallelComputations();
#endif /* */
      } else {
        if (this->mesh.describesAParallelComputation()) {
          return ctx.registerErrorMessage(
              "can't create a sequential finite element space on top of a "
              "parallel mesh");
        }
        auto p = this->sequential_fespaces.find(nc);
        if (p != this->sequential_fespaces.end()) {
          return p->second;
        }
        auto ptr = make_shared<FiniteElementSpace<false>>(
            ctx, this->mesh.getMutableMeshPointer<false>().get(),
            this->fec.get(), nc);
        this->sequential_fespaces.insert({nc, ptr});
        return ptr;
      }
    }  // end of getFiniteElementSpace
    /*!
     * \brief create a new parallel finite element space or reuse an existing
     * one
     * \param[in] ctx: execution context
     * \param[in] args: arguments defining the finite element space
     *
     * \note if a the list of materials identifiers contains the whole set of
     * material identifiers, the finite element space will be created on the
     * whole mesh and no submesh is created.
     *
     * \note if a sub mesh is created, it is stored internally by the underlying
     * mesh description.
     * \note if a finite element space is created, it is
     * stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(
        Context& ctx,
        const FiniteElementSpacesManager::
            GetFiniteElementSpaceOnSubMeshArguments& args) noexcept {
      const auto nc = args.number_of_components;
      if (nc < 1) {
        return ctx.registerErrorMessage("invalid number of components ('" +
                                        std::to_string(nc) + "')");
      }
      if (args.location == MeshDiscretization::Location::ON_MATERIALS) {
        const auto oids =
            this->mesh.getMaterialsIdentifiers(ctx, args.identifiers);
        if (isInvalid(oids)) {
          return {};
        }
        if (oids->empty()) {
          return ctx.registerErrorMessage("empty list of material identifiers");
        }
        const auto n = static_cast<size_type>(oids->size());
        if (n == getMaterialsAttributes(this->mesh).Size()) {
          return this->template getFiniteElementSpace<parallel>(ctx, nc);
        }
      }
      auto os = this->mesh.template getMutableSubMeshReference<parallel>(
          ctx, args.identifiers, args.location);
      if (isInvalid(os)) {
        return {};
      }
      return this->template getFiniteElementSpace<parallel>(ctx, *os, nc);
    }  // end of getFiniteElementSpace
    /*!
     * \brief create a new parallel finite element space or reuse an existing
     * one
     * \param[in] ctx: execution context
     * \param[in] m: mesh
     * \param[in] nc: number of components
     *
     * \note if a the list of materials identifiers contains the whole set of
     * material identifiers, the finite element space will be created on the
     * whole mesh and no submesh is created.
     *
     * \note if a sub mesh is created, it is stored internally by the underlying
     * mesh description.
     * \note if a finite element space is created, it is
     * stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context& ctx,
                          const Mesh<parallel>& m,
                          const size_type nc) noexcept {
      auto mptr = this->mesh.template getMutableMeshPointer<parallel>(ctx, m);
      if (isInvalid(mptr)) {
        return {};
      }
      auto omanager =
          this->template getFiniteElementSpacesManager<parallel>(ctx, m);
      if (isInvalid(omanager)) {
        return {};
      }
      auto p = omanager->find(nc);
      if (p == omanager->end()) {
        auto ptr = make_shared<FiniteElementSpace<parallel>>(
            ctx, mptr.get(), this->fec.get(), nc);
        if (isInvalid(ptr)) {
          return {};
        }
        omanager->insert({nc, ptr});
        return ptr;
      }
      return p->second;
    }  // end of getFiniteElementSpace

    /*!
     * \brief set of the nodal finite element space to the underlying mesh
     * \param[in] ctx: execution context
     *
     * \note if a scalar finite element space has already been declared, it is
     * reused.
     */
    [[nodiscard]] bool setNodalFiniteElementSpace(Context& ctx) noexcept {
      if (this->mesh.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        return this->template setNodalFiniteElementSpace<true>(
            ctx, *(this->mesh.getMutableMeshPointer<true>()));
#else
        reportUnsupportedParallelComputations();
#endif
      } else {
        return this->template setNodalFiniteElementSpace<false>(
            ctx, *(this->mesh.getMutableMeshPointer<false>()));
      }
    }  // end of setNodalFiniteElementSpace

    /*!
     * \brief set of the nodal finite element space to the underlying mesh
     * \param[in] ctx: execution context
     * \param[in] m: mesh
     *
     * \note the given mesh must be handled by the mesh discretization
     * \note if a scalar finite element space has already been declared, it is
     * reused.
     */
    template <bool parallel>
    [[nodiscard]] bool setNodalFiniteElementSpace(
        Context& ctx, const Mesh<parallel>& m) noexcept {
      const auto d = getSpaceDimension(this->mesh);
      auto mptr = this->mesh.template getMutableMeshPointer<parallel>(ctx, m);
      if (isInvalid(mptr)) {
        return false;
      }
      if constexpr (parallel) {
#ifdef MFEM_USE_MPI
        auto ptr = this->getFiniteElementSpace<true>(ctx, m, d);
        if (isInvalid(ptr)) {
          return false;
        }
        const auto* const nodes = m.GetNodes();
        if (nodes == nullptr) {
          mptr->SetNodalFESpace(ptr.get());
        } else {
          // nodes is a pointer to a grid function, even in parallel
          if (nodes->FESpace() != ptr.get()) {
            mptr->SetNodalFESpace(ptr.get());
          }
        }
#else
        reportUnsupportedParallelComputations();
#endif
      } else {
        auto ptr = this->getFiniteElementSpace<false>(ctx, m, d);
        if (isInvalid(ptr)) {
          return false;
        }
        const auto* const nodes = m.GetNodes();
        if (nodes == nullptr) {
          mptr->SetNodalFESpace(ptr.get());
        } else {
          if (nodes->FESpace() != ptr.get()) {
            mptr->SetNodalFESpace(ptr.get());
          }
        }
      }
      return true;
    }
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept {
      return *(this->fec);
    }  // end of getFiniteElementCollection
    //! \return the finite element collection
    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept {
      return this->fec;
    }  // end of getFiniteElementCollectionPointer
    /*!
     * \return if the given element space is also managed by this finite element
     * space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<true>& s) const noexcept {
      if (!this->mesh.describesAParallelComputation()) {
        return false;
      }
#ifdef MFEM_USE_MPI
      const auto nc = s.GetVDim();
      const auto* const m = s.GetParMesh();
      if (m == &(this->mesh.getMesh<true>())) {
        const auto p = this->parallel_fespaces.find(nc);
        if (p != this->parallel_fespaces.end()) {
          return p->second.get() == &s;
        }
        return false;
      }
      if (!this->mesh.manages(*m)) {
        return false;
      }
      const auto* const sm = dynamic_cast<const SubMesh<true>*>(m);
      if (sm == nullptr) {
        return false;
      }
      const auto pm = this->parallel_fespaces_on_submeshes.find(sm);
      if (pm == this->parallel_fespaces_on_submeshes.end()) {
        return false;
      }
      const auto pfes = pm->second.find(nc);
      if (pfes == pm->second.end()) {
        return false;
      }
      return pfes->second.get() == &s;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }  // end of manages
    /*!
     * \return if the given element space is also managed by this finite element
     * space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<false>& s) const noexcept {
      if (this->mesh.describesAParallelComputation()) {
        return false;
      }
      const auto nc = s.GetVDim();
      const auto* const m = s.GetMesh();
      if (m == &(this->mesh.getMesh<false>())) {
        const auto p = this->sequential_fespaces.find(nc);
        if (p == this->sequential_fespaces.end()) {
          return false;
        }
        return p->second.get() == &s;
      }
      if (!this->mesh.manages(*m)) {
        return false;
      }
      const auto* const sm = dynamic_cast<const SubMesh<false>*>(m);
      if (sm == nullptr) {
        return false;
      }
      const auto pm = this->sequential_fespaces_on_submeshes.find(sm);
      if (pm == this->sequential_fespaces_on_submeshes.end()) {
        return false;
      }
      const auto pfes = pm->second.find(nc);
      if (pfes == pm->second.end()) {
        return false;
      }
      return pfes->second.get() == &s;
    }  // end of manages

   private:
    /*!
     * \return the manager associated with a given mesh
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     *
     * \note The mesh must be managed by the mesh description
     */
    template <bool parallel>
    [[nodiscard]] OptionalReference<
        std::map<size_type, std::shared_ptr<FiniteElementSpace<parallel>>>>
    getFiniteElementSpacesManager(Context& ctx, const Mesh<parallel>& m) {
      if (!this->mesh.manages(m)) {
        return ctx.registerErrorMessage(
            "mesh is not managed by the underlying mesh description");
      }
      if constexpr (parallel) {
        if (!this->mesh.describesAParallelComputation()) {
          return ctx.registerErrorMessage(
              "can't create a manager of parallel finite element spaces on a "
              "sequential mesh");
        }
#ifdef MFEM_USE_MPI
        if (&(this->mesh.getMesh<true>()) == &m) {
          return {&(this->parallel_fespaces)};
        }
        const auto* const sm = dynamic_cast<const SubMesh<true>*>(&m);
        ctx.assertOrTerminate(
            sm != nullptr,
            "can't downcast the given mesh pointer to a sub mesh pointer");
        return {&(this->parallel_fespaces_on_submeshes[sm])};
#else
        reportUnsupportedParallelComputations();
#endif
      } else {
        if (this->mesh.describesAParallelComputation()) {
          return ctx.registerErrorMessage(
              "can't create a manager of sequential finite element spaces on a "
              "parallel mesh");
        }
        if (&(this->mesh.getMesh<false>()) == &m) {
          return {&(this->sequential_fespaces)};
        }
        const auto* const sm = dynamic_cast<const SubMesh<false>*>(&m);
        ctx.assertOrTerminate(
            sm != nullptr,
            "can't downcast the given mesh pointer to a sub mesh pointer");
        return {&(this->sequential_fespaces_on_submeshes[sm])};
      }
    }  // end of getFiniteElementSpacesManager
    //! \brief mesh
    MeshDiscretization mesh;
    //! \brief finite element collection
    std::shared_ptr<const FiniteElementCollection> fec;
    //! \brief parallel finite element spaces
#ifdef MFEM_USE_MPI
    std::map<size_type, std::shared_ptr<FiniteElementSpace<true>>>
        parallel_fespaces;
    std::map<const SubMesh<true>*,
             std::map<size_type, std::shared_ptr<FiniteElementSpace<true>>>>
        parallel_fespaces_on_submeshes;
#endif /* MFEM_USE_MPI */
    std::map<size_type, std::shared_ptr<FiniteElementSpace<false>>>
        sequential_fespaces;
    std::map<const SubMesh<false>*,
             std::map<size_type, std::shared_ptr<FiniteElementSpace<false>>>>
        sequential_fespaces_on_submeshes;
  };

  FiniteElementSpacesManager::FiniteElementSpacesManager(
      Context& ctx, const Parameters& parameters)
      : pimpl(std::make_shared<Implementation>(ctx, parameters)) {
  }  // end of FiniteElementSpacesManager

  FiniteElementSpacesManager::FiniteElementSpacesManager(
      Context&, const MeshDiscretization& m, const Parameters& parameters)
      : pimpl(std::make_shared<Implementation>(m, parameters)) {
  }  // end of FiniteElementSpacesManager

  FiniteElementSpacesManager::FiniteElementSpacesManager(
      Context&,
      const MeshDiscretization& m,
      std::shared_ptr<const FiniteElementCollection> c)
      : pimpl(std::make_shared<Implementation>(m, c)) {
  }  // end of FiniteElementSpacesManager

  FiniteElementSpacesManager::FiniteElementSpacesManager(
      FiniteElementSpacesManager&&) noexcept = default;

  FiniteElementSpacesManager::FiniteElementSpacesManager(
      const FiniteElementSpacesManager&) noexcept = default;

  MeshDiscretization FiniteElementSpacesManager::getMeshDiscretization()
      const noexcept {
    return this->pimpl->getMeshDiscretization();
  }  // end of getMeshDiscretization

  bool FiniteElementSpacesManager::setNodalFiniteElementSpace(
      Context& ctx) const noexcept {
    return this->pimpl->setNodalFiniteElementSpace(ctx);
  }  // end of setNodalFiniteElementSpace

  bool FiniteElementSpacesManager::setNodalFiniteElementSpace(
      Context& ctx, const Mesh<true>& m) const noexcept {
    return this->pimpl->setNodalFiniteElementSpace<true>(ctx, m);
  }  // end of setNodalFiniteElementSpace

  bool FiniteElementSpacesManager::setNodalFiniteElementSpace(
      Context& ctx, const Mesh<false>& m) const noexcept {
    return this->pimpl->setNodalFiniteElementSpace<false>(ctx, m);
  }  // end of setNodalFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<true>>
  FiniteElementSpacesManager::getParallelFiniteElementSpace(
      Context& ctx,
      const GetFiniteElementSpaceOnSubMeshArguments& args) const noexcept {
#ifdef MFEM_USE_MPI
    return this->pimpl->getFiniteElementSpace<true>(ctx, args);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of getParallelFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<false>>
  FiniteElementSpacesManager::getSequentialFiniteElementSpace(
      Context& ctx,
      const GetFiniteElementSpaceOnSubMeshArguments& args) const noexcept {
    return this->pimpl->getFiniteElementSpace<false>(ctx, args);
  }  // end of getSequentialFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<true>>
  FiniteElementSpacesManager::getParallelFiniteElementSpace(
      Context& ctx, const size_type nc) const noexcept {
#ifdef MFEM_USE_MPI
    return this->pimpl->getFiniteElementSpace<true>(ctx, nc);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of getParallelFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<true>>
  FiniteElementSpacesManager::getParallelFiniteElementSpace(
      Context& ctx, const Mesh<true>& m, const size_type nc) const noexcept {
#ifdef MFEM_USE_MPI
    return this->pimpl->getFiniteElementSpace<true>(ctx, m, nc);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of getParallelFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<false>>
  FiniteElementSpacesManager::getSequentialFiniteElementSpace(
      Context& ctx, const size_type nc) const noexcept {
    return this->pimpl->getFiniteElementSpace<false>(ctx, nc);
  }  // end of getSequentialFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<false>>
  FiniteElementSpacesManager::getSequentialFiniteElementSpace(
      Context& ctx, const Mesh<false>& m, const size_type nc) const noexcept {
    return this->pimpl->getFiniteElementSpace<false>(ctx, m, nc);
  }  // end of getSequentialFiniteElementSpace

  const FiniteElementCollection&
  FiniteElementSpacesManager::getFiniteElementCollection() const noexcept {
    return this->pimpl->getFiniteElementCollection();
  }  // end of getFiniteElementCollection

  std::shared_ptr<const FiniteElementCollection>
  FiniteElementSpacesManager::getFiniteElementCollectionPointer()
      const noexcept {
    return this->pimpl->getFiniteElementCollectionPointer();
  }  // end of getFiniteElementCollectionPointer

  bool FiniteElementSpacesManager::manages(
      const FiniteElementSpace<true>& s) const noexcept {
    return this->pimpl->manages(s);
  }  // end of manages

  bool FiniteElementSpacesManager::manages(
      const FiniteElementSpace<false>& s) const noexcept {
    return this->pimpl->manages(s);
  }  // end of manages

}  // end of namespace mfem_mgis