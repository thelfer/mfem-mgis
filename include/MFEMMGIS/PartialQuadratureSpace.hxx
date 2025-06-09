/*!
 * \file   PartialQuadratureSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX

#include <memory>
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
    size_type getNumberOfQuadraturePoints(const size_type) const;
    /*!
     * \brief return the hash table associating global element numbers and local
     * offsets.
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
                       size_type>  // offset
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
  size_type getSpaceSize(const PartialQuadratureSpace &);
  /*!
   * \brief return the number of quadrature points
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfIntegrationPoints`
   * \note this is as
   * requirement of mgis::function::QuadratureSpaceConcept
   */
  size_type getNumberOfElements(const PartialQuadratureSpace &);

  /*!
   * \brief return the number of finite elements associated identifier
   *
   * \note this function calls the method
   * `PartialQuadratureSpace::getNumberOfElements`
   * \note this is as requirement of mgis::function::QuadratureSpaceConcept
   */
  size_type getNumberOfCells(const PartialQuadratureSpace &);
  /*!
   * \brief return the number of quadrature points for the given finite
   * element
   * \param[in] e: index of the finite element
   */
  size_type getNumberOfQuadraturePoints(const PartialQuadratureSpace &,
                                        const size_type);
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

#include "MFEMMGIS/PartialQuadratureSpace.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX */
