/*!
 * \file   PartialQuadratureSpace.hxx
 * \brief
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

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief a space on quadrature points defined on a material
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureSpace {
    /*!
     * \brief throw an exception in case of invalid offset
     * \param[in] id: material identifier
     * \param[in] i: element number
     */
    [[noreturn]] static void treatInvalidElementIndex(const size_type,
                                                      const size_type);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] irs: function returning the order of quadrature for the
     * considered finite element.
     */
    PartialQuadratureSpace(const FiniteElementDiscretization &,
                           const size_type,
                           const std::function<const mfem::IntegrationRule &(
                               const mfem::FiniteElement &,
                               const mfem::ElementTransformation &)> &);
    //! \return the finite element discretization
    const FiniteElementDiscretization &getFiniteElementDiscretization() const;
    /*!
     * \return the integration ruel associated with the given finite element and
     * element transformation
     * \param[in] e: finite element
     * \param[in] tr: element transformation
     */
    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &, const mfem::ElementTransformation &) const;
    /*!
     * \brief return the number of finite element associated with this material
     * identifier
     */
    size_type getNumberOfElements() const;
    //! \brief return the number of integration points
    size_type getNumberOfIntegrationPoints() const;
    /*!
     * \brief return the number of quadrature points for the given finite
     * element
     * \param[in] e: index of the finite element
     */
    [[nodiscard]] size_type getNumberOfQuadraturePoints(const size_type) const;
    /*!
     * \brief return the number of quadrature points for the given finite
     * element
     * \param[in, out] ctx: execution context
     * \param[in] e: index of the finite element
     */
    [[nodiscard]] std::optional<size_type> getNumberOfQuadraturePoints(
        Context &, const size_type) const noexcept;
    /*!
     * \brief return the hash table associating global element numbers and
     * local offsets.
     */
    const std::unordered_map<size_type, size_type> &getOffsets() const;
    /*!
     * \brief return the offset associated with an element
     * \param[in] i: element number (global numbering)
     */
    size_type getOffset(const size_type) const;
    //! \return the material id
    size_type getId() const;
    //! \brief destructor
    ~PartialQuadratureSpace();

   private:
    //! \brief underlying finite element discretization
    const FiniteElementDiscretization &fe_discretization;
    /*!
     * \brief function returning the order of quadrature for the
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
    //! \brief material identifier
    size_type id;
    //! \brief number of integration points
    size_type ng;
  };  // end of struct PartialQuadratureSpace

}  // end of namespace mfem_mgis

#ifdef MGIS_FUNCTION_SUPPORT

namespace mfem_mgis {

  /*!
   * \brief return the number of integration points
   *
   * \note this method is equivalent to `getNumberOfIntegrationPoints`
   * \note this is as requirement of mgis::function::SpaceConcept
   */
  constexpr bool areEquivalent(const PartialQuadratureSpace &,
                               const PartialQuadratureSpace &) noexcept;

  /*!
   * \brief return the number of integration points
   *
   * \note this method is equivalent to `getNumberOfIntegrationPoints`
   * \note this is as requirement of mgis::function::SpaceConcept
   */
  [[nodiscard]] size_type getSpaceSize(const PartialQuadratureSpace &);
  /*!
   * \brief return the number of quadrature points
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfIntegrationPoints`
   * \note this is as
   * requirement of mgis::function::QuadratureSpaceConcept
   */
  [[nodiscard]] size_type getNumberOfElements(const PartialQuadratureSpace &);
  /*!
   * \brief return the number of finite elements associated identifier
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfElements`
   * \note this is as requirement of mgis::function::QuadratureSpaceConcept
   */
  [[nodiscard]] size_type getNumberOfCells(const PartialQuadratureSpace &);
  /*!
   * \brief return the number of quadrature points for the given finite
   * element
   * \param[in] e: index of the finite element
   */
  [[nodiscard]] size_type getNumberOfQuadraturePoints(
      const PartialQuadratureSpace &, const size_type);

}  // namespace mfem_mgis

namespace mgis::function {

  template <>
  struct SpaceTraits<mfem_mgis::PartialQuadratureSpace> {
    /*!
     * \brief a simple alias
     *
     * \note this is as requirement of mgis::function::SpaceConcept
     */
    using size_type = mfem_mgis::size_type;
    /*!
     * \brief a simple alias
     *
     * \note this is as requirement of mgis::function::ElementSpaceConcept
     */
    using element_index_type = mfem_mgis::size_type;
    /*!
     * \brief boolean stating that the integration points are stored from 0 to
     * size()-1
     *
     * \note this is as requirement of mgis::function::LinearElementSpaceConcept
     */
    static constexpr auto linear_element_indexing = true;
    /*!
     * \brief a simple alias
     *
     * \note this is as requirement of mgis::function::QuadratureSpaceConcept
     */
    using cell_index_type = mfem_mgis::size_type;
    /*!
     * \brief a simple alias
     *
     * \note this is as requirement of mgis::function::QuadratureSpaceConcept
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
    //! \brief identifier of the underlying material
    size_type identifier;
    //! \brief name of the material, if defined
    std::string name;
    //! \brief number of cells (finite elements)
    size_type number_of_cells;
    //! \brief number of quadrature points
    size_type number_of_quadrature_points;
    /*!
     * \brief mapping giving for each geometric type in the partial quadrature
     * space the number of elements
     */
    std::map<mfem::Geometry::Type, size_type>
        number_of_cells_by_geometric_type;
    /*!
     * \brief mapping giving for each geometric type in the partial quadrature
     * space the number of quadrature points
     */
    std::map<mfem::Geometry::Type, size_type>
        number_of_quadrature_points_by_geometric_type;
  };  // end of PartialQuadratureSpaceInformation

  /*!
   * \return information about the partial quadrature space on the current
   * process
   *
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<PartialQuadratureSpaceInformation>
  getLocalInformation(Context &, const PartialQuadratureSpace &) noexcept;
  /*!
   * \return information about the partial quadrature space, gathered from all
   * processes
   *
   * \param[in, out] ctx: execution context
   * \param[in] s: partial quadrature space
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<PartialQuadratureSpaceInformation> getInformation(
      Context &, const PartialQuadratureSpace &) noexcept;
  /*!
   * \brief write information, gathered from all processes in parallel, about
   * the partial quadrature space in the output stream
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] info: information to be displayed
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool info(
      Context &,
      std::ostream &,
      const PartialQuadratureSpaceInformation &) noexcept;
  /*!
   * \brief write information, gathered from all processes in parallel, about
   * the partial quadrature space in the output stream
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] s: partial quadrature space
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool info(
      Context &, std::ostream &, const PartialQuadratureSpace &) noexcept;
  /*!
   * \brief synchronize information of all processes
   *
   * \param[in, out] ctx: execution context
   * \param[in] info: information to be shared
   */
  std::optional<PartialQuadratureSpaceInformation> synchronize(
      Context &ctx, const PartialQuadratureSpaceInformation &) noexcept;

}  // end of  namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureSpace.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX */
