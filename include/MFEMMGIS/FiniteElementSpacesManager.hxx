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
   * \brief This class manages similar finite elements spaces (same mesh, same
   * finite element collection, but different vectorial dimensions), denoted as
   * siblings.
   *
   * This class is designed to be lightweight, movable and copyable
   */
  struct MFEM_MGIS_EXPORT FiniteElementSpacesManager {
    //! \brief string associated to the `FiniteElementFamily` parameter
    static const char* const FiniteElementFamily;
    //! \brief string associated to the `FiniteElementOrder` parameter
    static const char* const FiniteElementOrder;
    /*!
     * \return the list of parameters allowing to build a finite
     * element collection.
     *
     * Those parameters are used when the mesh discretization is already built;
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
     * \param[in] ctx: execution context
     * \param[in] p: parameters
     */
    FiniteElementSpacesManager(Context& ctx, const Parameters&);
    /*!
     * \brief constructor from a mesh discretization
     * \param[in] ctx: execution context
     * \param[m] mesh discretization
     * \param[in] p: parameters
     */
    FiniteElementSpacesManager(Context& ctx,
                               const MeshDiscretization&,
                               const Parameters&);
    /*!
     * \brief constructor from a mesh discretization
     * \param[in] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] fec: finite element collection
     */
    FiniteElementSpacesManager(Context& ctx,
                               const MeshDiscretization&,
                               std::shared_ptr<const FiniteElementCollection>);
    //! \brief move constructor
    FiniteElementSpacesManager(FiniteElementSpacesManager&&) noexcept;
    //! \brief copy constructor
    FiniteElementSpacesManager(const FiniteElementSpacesManager&) noexcept;
    //! \brief return the mesh discretization
    MeshDiscretization getMeshDiscretization() const noexcept;
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept;
    //! \return the finite element collection
    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept;
    /*!
     * \brief create a new finite element space or reuse an existing one
     * \param[in] ctx: execution context
     * \param[in] nc: vectorial dimension
     *
     * \note if a finite element space is created, it is stored internally.
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context&, const size_type) const noexcept;
    /*!
     * \brief assign a suitable nodal finite element space to the underlying
     * mesh \param[in] ctx: execution context
     *
     * \note if a scalar finite element space has already been declared, it is
     * reused.
     */
    [[nodiscard]] bool setNodalFiniteElementSpace(Context&) const noexcept;
    /*!
     * \return if the given element space is also managed by this finite element
     * space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<true>& s) const noexcept;
    /*!
     * \return if the given element space is also managed by this finite element
     * space manager
     * \param[in] s: finite element space
     */
    [[nodiscard]] bool manages(
        const FiniteElementSpace<false>& s) const noexcept;

   private:
    /*!
     * \brief create a parallel finite element space
     * \param[in] ctx: execution context
     * \param[in] nc: vectorial dimension
     */
    std::shared_ptr<FiniteElementSpace<true>> getParallelFiniteElementSpace(
        Context&, const size_type) const noexcept;
    /*!
     * \brief create a sequential finite element space
     * \param[in] ctx: execution context
     * \param[in] nc: vectorial dimension
     */
    [[nodiscard]] std::shared_ptr<FiniteElementSpace<false>>
    getSequentialFiniteElementSpace(Context&, const size_type) const noexcept;
    //! \internal structure to implement the PIMPL idiom
    struct Implementation;
    //! \brief pointer to the implementation
    std::shared_ptr<Implementation> pimpl;
  };  // end of struct FiniteElementSpacesManager;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementSpacesManager.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_HXX */
