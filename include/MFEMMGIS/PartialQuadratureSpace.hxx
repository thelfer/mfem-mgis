/*!
 * \file   PartialQuadratureSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX

#include <functional>
#include <unordered_map>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief a space on quadrature points defined on a material
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureSpace {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] integration_rule_selector: function returning the order of quadrature for the
     *   considered finite element.
     */
    PartialQuadratureSpace(const FiniteElementDiscretization &,
                           const size_type,
                           const std::function<const mfem::IntegrationRule &(
                               const mfem::FiniteElement &,
                               const mfem::ElementTransformation &)> &);
    //! \brief return the offset associated with an element
    size_type getNumberOfElements() const;
    //! \brief return the number of integration points
    size_type getNumberOfIntegrationPoints() const;
    //! \brief return the offset associated with an element
    size_type getOffset(const size_type) const;
    //! \brief destructor
    ~PartialQuadratureSpace();

   private:
    //! \brief underlying finite element space
    const FiniteElementDiscretization &fe_discretization;
    //! \brief offsets associated with elements
    std::unordered_map<size_type,  // element number
                       size_type>  // offset
        offsets;
    //! \brief material identifier
    size_type id;
    //! \brief number of integration points
    size_type ng;
  };  // end of struct QuadratureSpace

}  // end of mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_HXX */
