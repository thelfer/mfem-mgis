/*!
 * \file   include/MFEMMGIS/PartialQuadratureSpace.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_IXX
#define LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_IXX

namespace mfem_mgis {

  inline size_type PartialQuadratureSpace::getNumberOfElements() const {
    return this->offsets.size();
  }  // end of getNumberOfElements

  inline size_type PartialQuadratureSpace::getNumberOfIntegrationPoints()
      const {
    return this->ng;
  }  // end of getNumberOfIntegrationPoints

  inline size_type PartialQuadratureSpace::getId() const {
    return this->id;
  }  // end of getId

  inline const std::unordered_map<size_type, size_type>&
  PartialQuadratureSpace::getOffsets() const {
    return this->offsets;
  }  // end of getOffsets

  inline size_type PartialQuadratureSpace::getOffset(const size_type i) const {
    const auto p = this->offsets.find(i);
    if (p == this->offsets.end()) {
      PartialQuadratureSpace::treatInvalidElementIndex(this->id, i);
    }
    return p->second;
  }  // end of getOffset

  inline size_type PartialQuadratureSpace::getNumberOfQuadraturePoints(
      const size_type e) const {
    const auto p = this->number_of_quadrature_points.find(e);
    if (p == this->number_of_quadrature_points.end()) {
      PartialQuadratureSpace::treatInvalidElementIndex(this->id, e);
    }
    return p->second;
  }  // end of getNumberOfQuadraturePoints

  inline size_type getSpaceSize(const PartialQuadratureSpace& s) {
    return s.getNumberOfIntegrationPoints();
  }  // end of getSpaceSize

  inline size_type getNumberOfCells(const PartialQuadratureSpace& s) {
    return s.getNumberOfElements();
  }  // end of getNumberOfCells

  inline size_type getNumberOfQuadraturePoints(const PartialQuadratureSpace& s,
                                               const size_type e) {
    return s.getNumberOfQuadraturePoints(e);
  }

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_IXX */
