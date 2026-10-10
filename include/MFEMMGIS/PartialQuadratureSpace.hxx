/*!
 * \file   PartialQuadratureSpace.hxx
 * \brief  This file declares the `PartialQuadratureSpace` class
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX

#include <map>
#include <memory>
#include <iosfwd>
#include <variant>
#include <functional>
#include <unordered_map>

#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/SpaceConcept.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Info.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief a space on quadrature points defined on a material or a boundary
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureSpace {
    /*!
     * \brief throw an exception in case of invalid element index
     * \param[in] l: location
     * \param[in] i: element number
     */
    [[noreturn]] static void treatInvalidElementIndex(
        const LocationIdentifier l, const size_type i);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] l: location identifier
     * \param[in] irs: function returning the integration rule for the
     * considered finite element.
     */
    PartialQuadratureSpace(const FiniteElementDiscretization &fed,
                           const LocationIdentifier &l,
                           const std::function<const mfem::IntegrationRule &(
                               const mfem::FiniteElement &,
                               const mfem::ElementTransformation &)> &irs);

#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] fespace: finite element space
     * \param[in] l: location identifier
     * \param[in] irs: function returning the integration rule for the
     * considered finite element.
     */
    PartialQuadratureSpace(const FiniteElementDiscretization &fed,
                           const FiniteElementSpace<true> &fespace,
                           const LocationIdentifier l,
                           const std::function<const mfem::IntegrationRule &(
                               const mfem::FiniteElement &,
                               const mfem::ElementTransformation &)> &irs);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] fespace: finite element space
     * \param[in] l: location identifier
     * \param[in] irs: function returning the integration rule for the
     * considered finite element.
     */
    PartialQuadratureSpace(const FiniteElementDiscretization &fed,
                           const FiniteElementSpace<false> &fespace,
                           const LocationIdentifier l,
                           const std::function<const mfem::IntegrationRule &(
                               const mfem::FiniteElement &,
                               const mfem::ElementTransformation &)> &irs);

    //! \brief move constructor
    PartialQuadratureSpace(PartialQuadratureSpace &&) noexcept = default;
    PartialQuadratureSpace(const PartialQuadratureSpace &) = delete;
    PartialQuadratureSpace &operator=(PartialQuadratureSpace &&) = delete;
    PartialQuadratureSpace &operator=(const PartialQuadratureSpace &) = delete;
    //! \return the material name or boundary name
    [[nodiscard]] std::string getLocationName() const noexcept;
    //! \return if the partial quadrature space is defined on a material
    [[nodiscard]] bool isDefinedOnAMaterial() const;
    //! \return if the partial quadrature space is defined on a boundary
    [[nodiscard]] bool isDefinedOnABoundary() const;
    //! \return the mesh discretization
    [[nodiscard]] const MeshDiscretization &getMeshDiscretization()
        const noexcept;
    //! \return the finite element discretization
    [[nodiscard]] const FiniteElementDiscretization &
    getFiniteElementDiscretization() const noexcept;
    /*!
     * \brief return the underlying mesh
     * \return the underlying mesh on which the partial quadrature space is
     * built.
     * \param[in, out] ctx: execution context
     */
    template <bool parallel>
    [[nodiscard]] OptionalReference<const Mesh<parallel>> getMesh(
        Context &ctx) const noexcept;
    /*!
     * \brief return the finite element space
     * \return the finite element space on which the partial quadrature space is
     * built.
     * \param[in, out] ctx: execution context
     */
    template <bool parallel>
    [[nodiscard]] OptionalReference<const FiniteElementSpace<parallel>>
    getFiniteElementSpace(Context &ctx) const noexcept;
    /*!
     * \return if one shall iterate on boundary elements
     *
     * The rationale behind this method is that partial quadrature spaces can be
     * defined on the main mesh or on submeshes.
     *
     * If a partial quadrature space is defined on a boundary, two cases may
     * happen:
     *
     * 1. it can be created on the main mesh or a submesh defined on elements.
     *    In this case, one shall iterate over boundary elements and use
     *    the MFEM API relative to boundary elements (`GetNBE`,
     *    `GetBoundaryElement`, etc.)
     * 2. if a submesh is defined on the boundaries of the main mesh,
     *    however, one shall iterate over the elements  and use
     *    the MFEM API relative to standard elements (`GetNE`,
     *    `GetElement`, etc.)
     */
    [[nodiscard]] bool shallUseBoundaryElementsAPI() const noexcept;
    /*!
     * \brief return the integration rule of an element
     * \return the integration rule associated with the given finite element and
     * element transformation
     * \param[in] e: finite element
     * \param[in] tr: element transformation
     */
    [[nodiscard]] const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const;
    //! \return the number of finite elements of the partial quadrature space
    [[nodiscard]] size_type getNumberOfElements() const noexcept;
    //! \return the number of integration points
    [[nodiscard]] size_type getNumberOfIntegrationPoints() const noexcept;
    /*!
     * \brief return the number of quadrature points for the given finite
     * element
     * \param[in] e: index of the finite element
     * \return the number of quadrature points
     */
    [[nodiscard]] size_type getNumberOfQuadraturePoints(
        const size_type e) const;
    /*!
     * \brief return the number of quadrature points for the given finite
     * element
     * \param[in, out] ctx: execution context
     * \param[in] e: index of the finite element
     * \return the number of quadrature points
     */
    [[nodiscard]] std::optional<size_type> getNumberOfQuadraturePoints(
        Context &ctx, const size_type e) const noexcept;
    /*!
     * \return the hash table associating global element numbers and
     * local offsets.
     */
    [[nodiscard]] const std::unordered_map<size_type, size_type> &getOffsets()
        const noexcept;
    /*!
     * \brief return the offset associated with an element
     * \param[in] i: element number (global numbering)
     * \return the offset
     */
    [[nodiscard]] size_type getOffset(const size_type i) const;
    //! \return the material or boundary identifier
    [[nodiscard]] LocationIdentifier getLocation() const noexcept;
    //! \brief destructor
    ~PartialQuadratureSpace();

   private:
    /*!
     * \brief internal method shared by constructors
     * \param[in] throwing: dummy attribute to indicate that this function may
     * throw an exception
     */
    template <bool parallel>
    void initialize(attributes::Throwing throwing);
    /*!
     * \return the underlying mesh
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<const Mesh<true>> getParallelMesh(
        Context &ctx) const noexcept;
    /*!
     * \return the underlying mesh
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<const Mesh<false>> getSequentialMesh(
        Context &ctx) const noexcept;
    /*!
     * \return the underlying finite element space
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<const FiniteElementSpace<true>>
    getParallelFiniteElementSpace(Context &ctx) const noexcept;
    /*!
     * \return the underlying finite element space
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<const FiniteElementSpace<false>>
    getSequentialFiniteElementSpace(Context &ctx) const noexcept;
    //! \brief underlying finite element discretization
    const FiniteElementDiscretization &fe_discretization;
#ifdef MFEM_USE_MPI
    //! \brief underlying parallel finite element space
    const FiniteElementSpace<true> *const parallel_fespace = nullptr;
#endif /* MFEM_USE_MPI */
    //! \brief underlying sequential finite element space
    const FiniteElementSpace<false> *const sequential_fespace = nullptr;
    /*!
     * \brief function returning the integration rule for the
     * considered finite element.
     */
    std::function<const mfem::IntegrationRule &(
        const mfem::FiniteElement &, const mfem::ElementTransformation &)>
        integration_rule_selector;
    //! \brief offsets associated with elements
    std::unordered_map<size_type,  // element number (global numbering)
                       size_type>  // offset
        offsets;
    //! \brief number of quadrature points associated with elements
    std::unordered_map<size_type,  // element number (global numbering)
                       size_type>  // number of quadrature points
        number_of_quadrature_points;
    //! \brief location in the main mesh
    LocationIdentifier location;
    //! \brief number of integration points
    size_type ng;
  };  // end of struct PartialQuadratureSpace

}  // end of namespace mfem_mgis

#ifdef MGIS_FUNCTION_SUPPORT

namespace mfem_mgis {

  /*!
   * \brief check if two quadrature spaces are equivalent
   * \param[in] s1: first quadrature space
   * \param[in] s2: second quadrature space
   * \return if two quadrature spaces are equivalent
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool areEquivalent(
      const PartialQuadratureSpace &s1,
      const PartialQuadratureSpace &s2) noexcept;

  /*!
   * \brief return the number of integration points
   * \param[in] s: partial quadrature space
   * \return the number of integration points
   *
   * \note this method is equivalent to `getNumberOfIntegrationPoints`
   * \note this is a requirement of mgis::function::SpaceConcept
   */
  [[nodiscard]] size_type getSpaceSize(const PartialQuadratureSpace &s);
  /*!
   * \brief return the number of quadrature points
   * \param[in] s: partial quadrature space
   * \return the number of quadrature points
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfIntegrationPoints`
   * \note this is a
   * requirement of mgis::function::QuadratureSpaceConcept
   */
  [[nodiscard]] size_type getNumberOfElements(const PartialQuadratureSpace &s);
  /*!
   * \brief return the number of finite elements of the space
   * \param[in] s: partial quadrature space
   * \return the number of finite elements
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfElements`
   * \note this is a requirement of mgis::function::QuadratureSpaceConcept
   */
  [[nodiscard]] size_type getNumberOfCells(const PartialQuadratureSpace &s);
  /*!
   * \brief return the number of quadrature points for the given finite
   * element
   * \param[in] s: partial quadrature space
   * \param[in] e: index of the finite element
   * \return the number of quadrature points
   */
  [[nodiscard]] size_type getNumberOfQuadraturePoints(
      const PartialQuadratureSpace &s, const size_type e);

}  // namespace mfem_mgis

namespace mgis::function {

  //! \brief specialisation for partial quadrature spaces
  template <>
  struct SpaceTraits<mfem_mgis::PartialQuadratureSpace> {
    /*!
     * \brief a simple alias
     *
     * \note this is a requirement of mgis::function::SpaceConcept
     */
    using size_type = mfem_mgis::size_type;
    /*!
     * \brief a simple alias
     *
     * \note this is a requirement of mgis::function::ElementSpaceConcept
     */
    using element_index_type = mfem_mgis::size_type;
    /*!
     * \brief boolean stating that the integration points are stored from 0 to
     * size()-1
     *
     * \note this is a requirement of mgis::function::LinearElementSpaceConcept
     */
    static constexpr auto linear_element_indexing = true;
    /*!
     * \brief a simple alias
     *
     * \note this is a requirement of mgis::function::QuadratureSpaceConcept
     */
    using cell_index_type = mfem_mgis::size_type;
    /*!
     * \brief a simple alias
     *
     * \note this is a requirement of mgis::function::QuadratureSpaceConcept
     */
    using quadrature_point_index_type = mfem_mgis::size_type;
  };

  static_assert(SpaceConcept<mfem_mgis::PartialQuadratureSpace>);
  static_assert(ElementSpaceConcept<mfem_mgis::PartialQuadratureSpace>);
  static_assert(LinearElementSpaceConcept<mfem_mgis::PartialQuadratureSpace>);
  static_assert(QuadratureSpaceConcept<mfem_mgis::PartialQuadratureSpace>);

}  // end of namespace mgis::function

#endif /* MGIS_FUNCTION_SUPPORT */

namespace mfem_mgis {

  /*!
   * \brief structure describing information about a partial quadrature space
   */
  struct PartialQuadratureSpaceInformation {
    //! \brief identifier of the underlying material or boundary
    size_type identifier;
    //! \brief name of the material or boundary
    std::string name;
    //! \brief number of cells (finite elements)
    size_type number_of_cells;
    //! \brief number of quadrature points
    size_type number_of_quadrature_points;
    /*!
     * \brief mapping giving for each geometric type in the partial quadrature
     * space the number of elements
     */
    std::map<mfem::Geometry::Type, size_type> number_of_cells_by_geometric_type;
    /*!
     * \brief mapping giving for each geometric type in the partial quadrature
     * space the number of quadrature points
     */
    std::map<mfem::Geometry::Type, size_type>
        number_of_quadrature_points_by_geometric_type;
#ifdef MFEM_USE_MPI
    //! \brief communicator
    MPI_Comm communicator = MPI_COMM_WORLD;
#endif /* MFEM_USE_MPI */
  };   // end of PartialQuadratureSpaceInformation

  /*!
   * \brief return local information about a partial quadrature space
   * \return information about the partial quadrature space on the current
   * process
   *
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<PartialQuadratureSpaceInformation>
  getLocalInformation(Context &ctx, const PartialQuadratureSpace &s) noexcept;
  /*!
   * \brief return global information about a partial quadrature space
   * \return information about the partial quadrature space, gathered from all
   * processes
   *
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<PartialQuadratureSpaceInformation> getInformation(
      Context &ctx, const PartialQuadratureSpace &s) noexcept;
  /*!
   * \brief write the given information about a partial quadrature space in
   * the output stream
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] info: information to be displayed
   * \return true on success
   */
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  getInformation<PartialQuadratureSpaceInformation>(
      Context &ctx,
      std::ostream &os,
      const PartialQuadratureSpaceInformation &info) noexcept;
  /*!
   * \brief write information, gathered from all processes in parallel, about
   * the partial quadrature space in the output stream
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] s: partial quadrature space
   * \return true on success
   */
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool getInformation<PartialQuadratureSpace>(
      Context &ctx, std::ostream &os, const PartialQuadratureSpace &s) noexcept;
  /*!
   * \brief synchronize information of all processes
   *
   * \param[in, out] ctx: execution context
   * \param[in] info: information to be shared
   * \return the information gathered from all processes
   */
  std::optional<PartialQuadratureSpaceInformation> synchronize(
      Context &ctx, const PartialQuadratureSpaceInformation &info) noexcept;

}  // end of  namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureSpace.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX */
