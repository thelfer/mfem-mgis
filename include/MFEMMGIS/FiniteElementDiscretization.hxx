/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.hxx
 * \brief  This file declares the `FiniteElementDiscretization` class
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "MFEMMGIS/Info.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;
  struct Parameters;

  /*!
   * \brief a simple class used to:
   * - handle the life time of the mesh and the finite element collection.
   * - create and handle a finite element space
   */
  struct MFEM_MGIS_EXPORT FiniteElementDiscretization : MeshDiscretization {
    //! \brief string associated to the `UnknownsSize` parameter
    static const char* const UnknownsSize;
    //! \brief report that no parallel finite element space is defined
    [[noreturn]] static void reportInvalidParallelFiniteElementSpace();
    //! \brief report that no sequential finite element space is defined
    [[noreturn]] static void reportInvalidSequentialFiniteElementSpace();
    //! \return the list of valid parameters
    static std::vector<std::string> getParametersList();

    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters
     */
    FiniteElementDiscretization(Context& ctx, const Parameters& params);
    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] m: finite element spaces manager
     * \param[in] params: parameters
     */
    FiniteElementDiscretization(Context& ctx,
                                const FiniteElementSpacesManager& m,
                                const Parameters& params);
    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] params: parameters
     */
    FiniteElementDiscretization(Context& ctx,
                                const MeshDiscretization& m,
                                const Parameters& params);
    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] params: parameters
     */
    FiniteElementDiscretization(Context& ctx,
                                std::shared_ptr<Mesh<true>> m,
                                const Parameters& params);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     * used. Supported families are:
     * - `H1`:
     * The default value is `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknowns
     */
    FiniteElementDiscretization(Context& ctx,
                                std::shared_ptr<Mesh<false>> m,
                                const Parameters& params);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] c: finite element collection
     * \param[in] d: size of the unknowns
     *
     * \note this method creates the finite element space.
     */
    FiniteElementDiscretization(
        Context& ctx,
        std::shared_ptr<Mesh<true>> m,
        std::shared_ptr<const FiniteElementCollection> c,
        const size_type d);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] d: size of the unknowns
     *
     * \note this method creates the finite element space.
     */
    FiniteElementDiscretization(
        Context& ctx,
        std::shared_ptr<Mesh<false>> m,
        std::shared_ptr<const FiniteElementCollection> c,
        const size_type d);
    /*!
     * \brief check if a finite element space is managed by the same finite
     * element spaces manager
     * \return if the given element space is also managed by the finite element
     * space manager
     * \param[in] s: finite element space
     */
    template <bool parallel>
    [[nodiscard]] bool isSibling(
        const FiniteElementSpace<parallel>& s) const noexcept;
    //! \return the underlying finite element space manager
    [[nodiscard]] FiniteElementSpacesManager getFiniteElementSpacesManager()
        const noexcept;
    //! \return the finite element space
    template <bool parallel>
    [[nodiscard]] FiniteElementSpace<parallel>& getFiniteElementSpace();
    //! \return the finite element space
    template <bool parallel>
    [[nodiscard]] const FiniteElementSpace<parallel>& getFiniteElementSpace()
        const;
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept;
    //! \return the finite element collection
    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept;
    //
    // expose MeshDiscretization's methods, even deprecated ones for backward
    // compatibility
    //
    using MeshDiscretization::getBoundariesIdentifiers;
    using MeshDiscretization::getBoundariesNames;
    using MeshDiscretization::getBoundaryIdentifier;
    using MeshDiscretization::getBoundaryName;
    using MeshDiscretization::getMaterialIdentifier;
    using MeshDiscretization::getMaterialName;
    using MeshDiscretization::getMaterialsIdentifiers;
    using MeshDiscretization::getMaterialsNames;
    using MeshDiscretization::setBoundariesNames;
    using MeshDiscretization::setMaterialsNames;
    //! \brief destructor
    ~FiniteElementDiscretization();

   private:
    //! \brief manager of the finite element spaces
    FiniteElementSpacesManager fespaces_manager;
#ifdef MFEM_USE_MPI
    //! \brief parallel finite element space
    std::shared_ptr<FiniteElementSpace<true>> parallel_fe_space;
#endif /* MFEM_USE_MPI */
    //! \brief sequential finite element space
    std::shared_ptr<FiniteElementSpace<false>> sequential_fe_space;
  };  // end of FiniteElementDiscretization

#ifdef MFEM_USE_MPI
  /*!
   * \return the underlying mesh
   * \param[in] s: finite element space
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const Mesh<true>& getMesh(
      const FiniteElementSpace<true>&) noexcept;
#endif /* MFEM_USE_MPI */
  /*!
   * \return the underlying mesh
   * \param[in] s: finite element space
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const Mesh<false>& getMesh(
      const FiniteElementSpace<false>&) noexcept;
  /*!
   * \brief return the number of components of the unknowns
   * \return the number of components (vectorial dimension) of the
   * underlying finite element space. \param[in] fed: finite element
   * discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getNumberOfComponents(
      const FiniteElementDiscretization& fed) noexcept;

  /*!
   * \brief return the number of unknowns
   * \return the total number of unknowns of the underlying
   * finite element space, including those required to handle ghost values or
   * hanging nodes.
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getVSize(
      const FiniteElementDiscretization& fed) noexcept;

  /*!
   * \brief return the number of true unknowns
   * \return the number of unknowns of the underlying finite element space,
   * excluding those required to handle ghost values or hanging nodes.
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getTrueVSize(
      const FiniteElementDiscretization& fed) noexcept;

  /*!
   * \brief display information about a finite element discretization
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] fed: finite element discretization
   * \return true on success
   */
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  getInformation<FiniteElementDiscretization>(
      Context& ctx,
      std::ostream& os,
      const FiniteElementDiscretization& fed) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX */
