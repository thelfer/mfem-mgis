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
  }  // end of PartialQuadratureSpace::getNumberOfElements

  inline size_type PartialQuadratureSpace::getNumberOfIntegrationPoints()
      const {
    return this->ng;
  }  // end of PartialQuadratureSpace::getNumberOfIntegrationPoints

  inline size_type PartialQuadratureSpace::getOffset(const size_type i) const {
    const auto p = this->offsets.find(i);
    if (p == this->offsets.end()) {
      PartialQuadratureSpace::treatInvalidOffset(this->id, i);
    }
    return p->second;
  }  // end of PartialQuadratureSpace::getOffset

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATURESPACE_IXX */
