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
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return if two partial quadrature spaces are defined on the same mesh, or
   * on two meshes for which the elements of the location of the first space
   * can be mapped, see `MeshDiscretization::getElementsMapping`
   * \param[in] s1: first quadrature space
   * \param[in] s2: second quadrature space
   *
   * \note this call is collective when the meshes are distinct.
   */
  [[nodiscard]] static bool areQuadratureSpacesOnMappableMeshes(
      const PartialQuadratureSpaceImplementation& s1,
      const PartialQuadratureSpaceImplementation& s2) noexcept {
    auto ctx = Context{};
    const auto fed = s1.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto om1 = s1.getMesh<true>(ctx);
      const auto om2 = s2.getMesh<true>(ctx);
      if (!areValid(om1, om2)) {
        return false;
      }
      if (&(*om1) == &(*om2)) {
        return true;
      }
      // collective call
      const auto omapping =
          fed.getElementsMapping(ctx, *om1, *om2, s1.getLocation());
      return isValid(omapping);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto om1 = s1.getMesh<false>(ctx);
    const auto om2 = s2.getMesh<false>(ctx);
    if (!areValid(om1, om2)) {
      return false;
    }
    if (&(*om1) == &(*om2)) {
      return true;
    }
    // collective call
    const auto omapping =
        fed.getElementsMapping(ctx, *om1, *om2, s1.getLocation());
    return isValid(omapping);
  }  // end of areQuadratureSpacesOnMappableMeshes

  static bool areMappable(
      const PartialQuadratureSpaceImplementation& s1,
      const PartialQuadratureSpaceImplementation& s2) noexcept {
    if (&s1 == &s2) {
      return true;
    }
    const auto& fed1 = s1.getFiniteElementDiscretization();
    const auto& fed2 = s2.getFiniteElementDiscretization();
    if (fed1.getFiniteElementSpacesManager() !=
        fed2.getFiniteElementSpacesManager()) {
      return false;
    }
    if (s1.getLocation() != s2.getLocation()) {
      return false;
    }
    const auto success = [&s1, &s2] {
      if (s1.getNumberOfIntegrationPoints() !=
          s2.getNumberOfIntegrationPoints()) {
        return false;
      }
      const auto& offsets1 = s1.getOffsets();
      const auto& offsets2 = s2.getOffsets();
      if (offsets1.size() != offsets2.size()) {
        return false;
      }
      auto p1 = offsets1.begin();
      auto p2 = offsets2.begin();
      for (; p1 != offsets1.end(); ++p1, ++p2) {
        if ((p1->first != p2->first) || (p1->second != p2->second)) {
          return false;
        }
      }
      return true;
    }();
    if (!isTrueOnAllProcesses(s1.getMeshDiscretization(), success)) {
      return false;
    }
    // checked last since it may build the elements mapping
    return areQuadratureSpacesOnMappableMeshes(s1, s2);
  }  // end of areMappable

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  const char* const FiniteElementSpacesManager::FiniteElementFamily =
      "FiniteElementFamily";
  const char* const FiniteElementSpacesManager::FiniteElementOrder =
      "FiniteElementOrder";
  const char* const FiniteElementSpacesManager::FiniteElementSpaceOrdering =
      "FiniteElementSpaceOrdering";

  std::vector<std::string>
  FiniteElementSpacesManager::getFiniteElementCollectionParametersList() {
    return {FiniteElementSpacesManager::FiniteElementFamily,
            FiniteElementSpacesManager::FiniteElementOrder,
            FiniteElementSpacesManager::FiniteElementSpaceOrdering};
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

  /*!
   * \return the ordering of the degrees of freedom of the finite element spaces
   * \param[in] params: parameters
   */
  [[nodiscard]] static mfem::Ordering::Type getFiniteElementSpaceOrdering(
      attributes::Throwing, const Parameters& params) {
    const auto ordering = get_if<std::string>(
        throwing, params,
        FiniteElementSpacesManager::FiniteElementSpaceOrdering, "byNODES");
    if (ordering == "byNODES") {
      return mfem::Ordering::byNODES;
    }
    if (ordering != "byVDIM") {
      raise(
          "FiniteElementSpacesManager::FiniteElementSpacesManager: "
          "invalid finite element space ordering '" +
          ordering + "' (expected 'byNODES' or 'byVDIM')");
    }
    return mfem::Ordering::byVDIM;
  }  // end of getFiniteElementSpaceOrdering

  struct FiniteElementSpacesManager::Implementation {
    /*!
     * \brief constructor from parameters
     * \param[in, out] ctx: execution context
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
      this->ordering = getFiniteElementSpaceOrdering(throwing, parameters);
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
      this->ordering = getFiniteElementSpaceOrdering(throwing, parameters);
    }  // end of Implementation
    /*!
     * \brief constructor from a mesh discretization and a finite element
     * collection
     * \param[in] m: mesh discretization
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
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
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
            this->fec.get(), nc, this->ordering);
        if (isInvalid(ptr)) {
          return {};
        }
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
            this->fec.get(), nc, this->ordering);
        if (isInvalid(ptr)) {
          return {};
        }
        this->sequential_fespaces.insert({nc, ptr});
        return ptr;
      }
    }  // end of getFiniteElementSpace
    /*!
     * \brief create a new finite element space or reuse an existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] args: arguments defining the finite element space
     *
     * \note if the list of materials identifiers contains the whole set of
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
     * \brief create a new finite element space on the given mesh or reuse an
     * existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] nc: number of components
     *
     * \note the given mesh must be handled by the mesh discretization
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
            ctx, mptr.get(), this->fec.get(), nc, this->ordering);
        if (isInvalid(ptr)) {
          return {};
        }
        omanager->insert({nc, ptr});
        return ptr;
      }
      return p->second;
    }  // end of getFiniteElementSpace

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
      if (!this->mesh.manages(*m)) {
        return false;
      }
      if (m == &(this->mesh.getMesh<true>())) {
        const auto p = this->parallel_fespaces.find(nc);
        if (p != this->parallel_fespaces.end()) {
          return p->second.get() == &s;
        }
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
      if (!this->mesh.manages(*m)) {
        return false;
      }
      if (m == &(this->mesh.getMesh<false>())) {
        const auto p = this->sequential_fespaces.find(nc);
        if (p == this->sequential_fespaces.end()) {
          return false;
        }
        return p->second.get() == &s;
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

    /*!
     * \return the identifier of a quadrature space in the list of
     * quadrature spaces of a given location
     *
     * \param[in] s: quadrature space
     */
    [[nodiscard]] std::optional<size_type> getIdentifier(
        Context& ctx,
        const std::shared_ptr<
            const ::mfem_mgis::internals::PartialQuadratureSpaceImplementation>&
            s) noexcept {
      using PerLocationIdentifiersList = std::vector<std::vector<
          std::shared_ptr<const ::mfem_mgis::internals::
                              PartialQuadratureSpaceImplementation>>>;
      if (isInvalid(s)) {
        return ctx.registerErrorMessage("invalid partial quadrature space");
      }
      if (s->getMeshDiscretization() != this->mesh) {
        return ctx.registerErrorMessage(
            "partial quadrature space is not defined on the proper mesh");
      }
      if (this->mesh.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        const auto& ofes = s->getFiniteElementSpace<true>(ctx);
        if (isInvalid(ofes)) {
          return {};
        }
        if (!(this->manages(*ofes))) {
          return ctx.registerErrorMessage(
              "partial quadrature space is not defined on a finite element "
              "space managed by this finite element spaces manager");
        }
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      } else {
        const auto& ofes = s->getFiniteElementSpace<false>(ctx);
        if (isInvalid(ofes)) {
          return {};
        }
        if (!(this->manages(*ofes))) {
          return ctx.registerErrorMessage(
              "partial quadrature space is not defined on a finite element "
              "space managed by this finite element spaces manager");
        }
      }
      const auto l = s->getLocation();
      auto p = this->quadraturesSpaceIdentifiersStorage.find(l);
      if (p == this->quadraturesSpaceIdentifiersStorage.end()) {
        p = this->quadraturesSpaceIdentifiersStorage
                .insert({l, PerLocationIdentifiersList{}})
                .first;
      }
      return this->getIdentifier(p->second, s);
    }  // end of getIdendifier

    /*!
     * \return if two quadrature space are equivalent
     *
     * \param[in, out] ctx: execution context
     * \param[in] s1: first quadrature space
     * \param[in] s2: second quadrature space
     */
    [[nodiscard]] std::optional<bool> areMappable(
        Context& ctx,
        const std::shared_ptr<
            const ::mfem_mgis::internals::PartialQuadratureSpaceImplementation>&
            qspace1,
        const std::shared_ptr<
            const ::mfem_mgis::internals::PartialQuadratureSpaceImplementation>&
            qspace2) noexcept {
      if (!areValid(qspace1, qspace2)) {
        return ctx.registerErrorMessage("invalid quadrature spaces");
      }
      if (qspace1->getLocation() != qspace2->getLocation()) {
        return false;
      }
      const auto oid1 = this->getIdentifier(ctx, qspace1);
      const auto oid2 = this->getIdentifier(ctx, qspace2);
      if (!areValid(oid1, oid2)) {
        return {};
      }
      return *oid1 == *oid2;
    }  // end of areMappable

   private:
    /*!
     * \return the identifier of a quadrature space in the list of
     * quadrature spaces of a given location
     *
     * \param[in, out] identifiers: list of quadrature spaces
     * \param[in] s: quadrature space
     */
    [[nodiscard]] size_type getIdentifier(
        std::vector<std::vector<std::shared_ptr<
            const ::mfem_mgis::internals::
                PartialQuadratureSpaceImplementation>>>& identifiers,
        const std::shared_ptr<
            const ::mfem_mgis::internals::PartialQuadratureSpaceImplementation>&
            s) noexcept {
      // check if the quadrature space has already been declared
      for (size_type id = 0; const auto& spaces : identifiers) {
        for (const auto& s2 : spaces) {
          if (s.get() == s2.get()) {
            return id;
          }
        }
        ++id;
      }
      // add the quadrature space
      for (size_type id = 0; auto& spaces : identifiers) {
        assert(!spaces.empty());  // by design, spaces can't be empty
        if (::mfem_mgis::internals::areMappable(*s, *(spaces.at(0)))) {
          spaces.push_back(s);
          return id;
        }
        ++id;
      }
      // the quadrature space is not equivalent to any of the previously
      // registred space
      const auto id = identifiers.size();
      identifiers.push_back(
          std::vector<
              std::shared_ptr<const ::mfem_mgis::internals::
                                  PartialQuadratureSpaceImplementation>>(1, s));
      return id;
    }  // end of getIdentifier
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
    //! \brief ordering of the degrees of freedom of the finite element spaces
    mfem::Ordering::Type ordering = mfem::Ordering::byNODES;
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
    //! \brief a simple alias
    using IdentifiersMap =
        std::map<LocationIdentifier,
                 std::vector<std::vector<std::shared_ptr<
                     const ::mfem_mgis::internals::
                         PartialQuadratureSpaceImplementation>>>>;
    //! \brief list of equivalent quadrature spaces, sorted by locations
    IdentifiersMap quadraturesSpaceIdentifiersStorage;
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

  std::optional<bool> FiniteElementSpacesManager::areMappable(
      Context& ctx,
      const PartialQuadratureSpace& s1,
      const PartialQuadratureSpace& s2) const noexcept {
    const auto& fed1 = s1.getFiniteElementDiscretization();
    const auto& fed2 = s2.getFiniteElementDiscretization();
    const auto& m1 = fed1.getFiniteElementSpacesManager();
    const auto& m2 = fed2.getFiniteElementSpacesManager();
    if ((*this != m1) || (*this != m2)) {
      return false;
    }
    return this->pimpl->areMappable(ctx, s1.getImplementation(),
                                    s2.getImplementation());
  }  // end of areMappable

  std::optional<size_type> FiniteElementSpacesManager::getIdentifier(
      Context& ctx, const PartialQuadratureSpace& s) const noexcept {
    return this->pimpl->getIdentifier(ctx, s.getImplementation());
  }  // end of getIdentifier

  bool FiniteElementSpacesManager::operator==(
      const FiniteElementSpacesManager& rhs) const noexcept {
    return this->pimpl.get() == rhs.pimpl.get();
  }  // end of operator==

}  // end of namespace mfem_mgis