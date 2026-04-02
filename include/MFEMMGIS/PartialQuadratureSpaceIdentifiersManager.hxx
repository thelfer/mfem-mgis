/*!
 * \file   MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx
 * \brief  This file declares the `PartialQuadratureSpaceIdentifiersManager`
 * class
 * \author Thomas Helfer
 * \date   29/03/2026
 */

#ifndef LIB_MFEMMGIS_PARTIALQUADRATURESPACEIDENTIFIERSMANAGER_HXX
#define LIB_MFEMMGIS_PARTIALQUADRATURESPACEIDENTIFIERSMANAGER_HXX

#include <memory>
#include <unordered_map>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  /*!
   * \brief a class aiming at assigning an unique identifier to
   * a set of equivalent partial quadrature spaces.
   *
   * In short, two partial quadrature spaces are equivalent if they define the
   * same integration points on the same material id, see the `areEquivalent`
   * function for details.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureSpaceIdentifiersManager {
    /*!
     * \brief constructor
     * \param[in] m: mesh discretization
     */
    PartialQuadratureSpaceIdentifiersManager(
        const MeshDiscretization&) noexcept;
    /*!
     * \return the identifier associated with the given partial quadrature
     * space.
     *
     * \param[in, out] ctx: execution context
     * \param[in] s: partial quadrature space
     *
     * \note This identifier may be used for distinct material identifiers.
     */
    [[nodiscard]] std::optional<size_type> getIdentifier(
        Context&,
        const std::shared_ptr<const PartialQuadratureSpace>&) const noexcept;
    /*!
     * \brief return if two partial quadrature spaces are equivalent
     *
     * \param[in] qspace1: first partial quadrature space
     * \param[in] qspace2: second partial quadrature space
     *
     * \note if quadrature spaces are defined on two distinct mesh
     * discretizations, this method returns false
     */
    [[nodiscard]] bool areEquivalent(
        const std::shared_ptr<const PartialQuadratureSpace>&,
        const std::shared_ptr<const PartialQuadratureSpace>&) const noexcept;
    //! \brief destructor
    ~PartialQuadratureSpaceIdentifiersManager() noexcept;

   private:
    //
    PartialQuadratureSpaceIdentifiersManager() = delete;
    PartialQuadratureSpaceIdentifiersManager(
        PartialQuadratureSpaceIdentifiersManager&&) = delete;
    PartialQuadratureSpaceIdentifiersManager(
        const PartialQuadratureSpaceIdentifiersManager&&) = delete;
    PartialQuadratureSpaceIdentifiersManager& operator=(
        PartialQuadratureSpaceIdentifiersManager&&) = delete;
    PartialQuadratureSpaceIdentifiersManager& operator=(
        const PartialQuadratureSpaceIdentifiersManager&&) = delete;
    //! \brief a simple alias
    using IdentifiersMap =
        std::unordered_map<size_type,
                           std::vector<std::vector<
                               std::shared_ptr<const PartialQuadratureSpace>>>>;
    //! \brief mesh discretization on which the partial quadrature spaces are
    //! defined
    const MeshDiscretization mesh;
    /*!
     * \brief list of equivalent quadrature spaces, sorted by material
     * identifier
     *
     *
     */
    mutable IdentifiersMap ids;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARTIALQUADRATURESPACEIDENTIFIERSMANAGER_HXX */
