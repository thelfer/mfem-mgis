/*!
 * \file   MFEMMGIS/FiniteElementSpacesManager.hxx
 * \brief  This file declares the `FiniteElementSpacesManager` class
 * \author Thomas Helfer
 * \date   09/09/2026
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_HXX
#define LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  /*!
   * \brief This class manages similar finite element spaces, denoted as
   * siblings. They share the same finite element collection and are defined on
   * the same mesh or on one of its submeshes, but may have different vectorial
   * dimensions. The degrees of freedom of all those finite element spaces are
   * ordered in the same way (see the `FiniteElementSpaceOrdering` parameter).
   *
   * This class is designed to be lightweight, movable and copyable
   */
  struct MFEM_MGIS_EXPORT FiniteElementSpacesManager {
    //! \brief string associated to the `FiniteElementFamily` parameter
    static const char* const FiniteElementFamily;
    //! \brief string associated to the `FiniteElementOrder` parameter
    static const char* const FiniteElementOrder;
    /*!
     * \brief string associated to the `FiniteElementSpaceOrdering` parameter
     *
     * This parameter selects the ordering of the degrees of freedom of the
     * finite element spaces:
     *
     * - `byNODES` (default): all the values of the first component, then all
     *   the values of the second component, etc. (`XX...YY...ZZ...`).
     * - `byVDIM`: all the components of the first node, then all the
     *   components of the second node, etc. (`XYZXYZ...`).
     *
     * \note the `Elasticity` strategy of the `HypreBoomerAMG` preconditioner
     * requires the `byVDIM` ordering.
     */
    static const char* const FiniteElementSpaceOrdering;
    /*!
     * \return the list of parameters allowing to build a finite
     * element collection and to select the ordering of the degrees of freedom
     * of the finite element spaces.
     *
     * Those parameters are used when the mesh discretization is already built.
     */
    [[nodiscard]] static std::vector<std::string>
    getFiniteElementCollectionParametersList();
    /*!
     * \return the list of parameters allowing to build both a mesh and a finite
     * element space manager
     */
    [[nodiscard]] static std::vector<std::string> getParametersList();
    /*!
     * \brief constructor from parameters
     * \param[in, out] ctx: execution context
     * \param[in] parameters: parameters
     */
    FiniteElementSpacesManager(Context& ctx, const Parameters& parameters);
    /*!
     * \brief constructor from a mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] parameters: parameters
     */
    FiniteElementSpacesManager(Context& ctx,
                               const MeshDiscretization& m,
                               const Parameters& parameters);
    /*!
     * \brief constructor from a mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] c: finite element collection
     *
     * \note the degrees of freedom are ordered by nodes (`byNODES`)
     */
    FiniteElementSpacesManager(
        Context& ctx,
        const MeshDiscretization& m,
        std::shared_ptr<const FiniteElementCollection> c);
    //! \brief move constructor
    FiniteElementSpacesManager(FiniteElementSpacesManager&&) noexcept;
    //! \brief copy constructor
    FiniteElementSpacesManager(const FiniteElementSpacesManager&) noexcept;
    //! \return the mesh discretization
    MeshDiscretization getMeshDiscretization() const noexcept;
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept;
    //! \return the finite element collection
    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept;
    /*!
     * \brief structure used to create a finite element space on a submesh
     *
     * \see `getFiniteElementSpace` for details
     */
    struct GetFiniteElementSpaceOnSubMeshArguments {
      //! \brief location
      MeshDiscretization::Location location;
      /*!
       * \brief parameter used to identify the materials or the boundaries on
       * which the SubMesh is defined
       *
       * \see `MeshDiscretization::getSubMesh` for details
       */
      Parameter identifiers;
      /*!
       * \brief number of components (vectorial dimension) of the finite
       * element space (must be greater than or equal to 1).
       */
      size_type number_of_components;
    };  // end of struct GetFiniteElementSpaceOnSubMeshArguments
    /*!
     * \brief create a new finite element space on the whole mesh or reuse an
     * existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] nc: vectorial dimension
     *
     * \note if a finite element space is created, it is stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context& ctx, const size_type nc) const noexcept;
    /*!
     * \brief create a new finite element space or reuse an existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh on which the finite element space is defined
     * \param[in] nc: vectorial dimension
     *
     * \note the given mesh must be handled by the mesh discretization
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context& ctx,
                          const Mesh<parallel>& m,
                          const size_type nc) const noexcept;
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
     * mesh discretization.
     * \note if a finite element space is created, it is
     * stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(
        Context& ctx,
        const GetFiniteElementSpaceOnSubMeshArguments& args) const noexcept;
    /*!
     * \brief check if a finite element space is managed by this manager
     * \return if the given element space is also managed by this finite
     * element space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<true>& s) const noexcept;
    /*!
     * \brief check if a finite element space is managed by this manager
     * \return if the given element space is also managed by this finite
     * element space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<false>& s) const noexcept;

   private:
    /*!
     * \brief create a new parallel finite element space or reuse an existing
     * one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] nc: vectorial dimension
     */
    std::shared_ptr<FiniteElementSpace<true>> getParallelFiniteElementSpace(
        Context& ctx, const size_type nc) const noexcept;
    /*!
     * \brief create a new sequential finite element space or reuse an
     * existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] nc: vectorial dimension
     */
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<false>>
    getSequentialFiniteElementSpace(Context& ctx,
                                    const size_type nc) const noexcept;
    /*!
     * \brief create a new parallel finite element space or reuse an existing
     * one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] nc: vectorial dimension
     */
    std::shared_ptr<FiniteElementSpace<true>> getParallelFiniteElementSpace(
        Context& ctx, const Mesh<true>& m, const size_type nc) const noexcept;
    /*!
     * \brief create a new sequential finite element space or reuse an
     * existing one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] nc: vectorial dimension
     */
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<false>>
    getSequentialFiniteElementSpace(Context& ctx,
                                    const Mesh<false>& m,
                                    const size_type nc) const noexcept;
    /*!
     * \brief create a new parallel finite element space or reuse an existing
     * one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] args: arguments defining the finite element space
     *
     * \note if the list of materials identifiers contains the whole set of
     * material identifiers, the finite element space will be created on the
     * whole mesh and no submesh is created.
     *
     * \note if a sub mesh is created, it is stored internally by the underlying
     * mesh discretization.
     * \note if a finite element space is created, it is
     * stored internally.
     */
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<true>>
    getParallelFiniteElementSpace(
        Context& ctx,
        const GetFiniteElementSpaceOnSubMeshArguments& args) const noexcept;
    /*!
     * \brief create a new sequential finite element space or reuse an existing
     * one
     * \return the finite element space, a null pointer on failure
     * \param[in, out] ctx: execution context
     * \param[in] args: arguments defining the finite element space
     *
     * \note if the list of materials identifiers contains the whole set of
     * material identifiers, the finite element space will be created on the
     * whole mesh and no submesh is created.
     *
     * \note if a sub mesh is created, it is stored internally by the underlying
     * mesh discretization.
     * \note if a finite element space is created, it is
     * stored internally.
     */
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<false>>
    getSequentialFiniteElementSpace(
        Context& ctx,
        const GetFiniteElementSpaceOnSubMeshArguments& args) const noexcept;
    //! \internal structure to implement the PIMPL idiom
    struct Implementation;
    //! \brief pointer to the implementation
    std::shared_ptr<Implementation> pimpl;
  };  // end of struct FiniteElementSpacesManager;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementSpacesManager.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_HXX */
