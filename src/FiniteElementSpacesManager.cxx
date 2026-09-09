/*!
 * \file   src/FiniteElementSpacesManager.cxx
 * \brief  this file implements the `FiniteElementSpacesManager` class
 * \author Thomas Helfer
 * \date   09/09/2026
 */

#include <mfem/mesh/mesh.hpp>
#include <mfem/fem/fespace.hpp>
#ifdef MFEM_USE_MPI
#include <mfem/mesh/pmesh.hpp>
#include <mfem/fem/pfespace.hpp>
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
     * \brief set of the nodal finite element space to the underlying mesh
     * \param[in] ctx: execution context
     *
     * \note if a scalar finite element space has already been declared, it is
     * reused.
     */
    [[nodiscard]] bool setNodalFiniteElementSpace(Context& ctx) noexcept {
      const auto d = getSpaceDimension(this->mesh);
      if (this->mesh.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        auto& m = *(this->mesh.getMutableMeshPointer<true>());
        auto ptr = this->getFiniteElementSpace<true>(ctx, d);
        if (isInvalid(ptr)) {
          return false;
        }
        m.SetNodalFESpace(ptr.get());
#else
        reportUnsupportedParallelComputations();
#endif /* */
      } else {
        auto& m = *(this->mesh.getMutableMeshPointer<false>());
        auto ptr = this->getFiniteElementSpace<false>(ctx, d);
        if (isInvalid(ptr)) {
          return false;
        }
        m.SetNodalFESpace(ptr.get());
      }
      return true;
    }
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept {
      return *(this->fec);
    }  // end of getFiniteElementCollection

    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept {
      return this->fec;
    }  // end of getFiniteElementCollectionPointer

   private:
    //! \brief mesh
    MeshDiscretization mesh;
    //! \brief finite element collection
    std::shared_ptr<const FiniteElementCollection> fec;
    //! \brief parallel finite element spaces
#ifdef MFEM_USE_MPI
    std::map<size_type, std::shared_ptr<FiniteElementSpace<true>>>
        parallel_fespaces;
#endif /* MFEM_USE_MPI */
    std::map<size_type, std::shared_ptr<FiniteElementSpace<false>>>
        sequential_fespaces;
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

  std::shared_ptr<FiniteElementSpace<true>>
  FiniteElementSpacesManager::getParallelFiniteElementSpace(
      Context& ctx, const size_type nc) const noexcept {
    return this->pimpl->getFiniteElementSpace<true>(ctx, nc);
  }  // end of getParallelFiniteElementSpace

  std::shared_ptr<FiniteElementSpace<false>>
  FiniteElementSpacesManager::getSequentialFiniteElementSpace(
      Context& ctx, const size_type nc) const noexcept {
    return this->pimpl->getFiniteElementSpace<false>(ctx, nc);
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

}  // end of namespace mfem_mgis