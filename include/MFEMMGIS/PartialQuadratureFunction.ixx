/*!
 * \file   include/MFEMMGIS/PartialQuadratureFunction.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX

namespace mfem_mgis {

  inline size_type PartialQuadratureFunction::getDataOffset(const size_type o) const {
    return o * (this->data_stride) + this->data_begin;
  } // end of getDataOffset

  inline real& PartialQuadratureFunction::getIntegrationPointValue(const size_type o) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline const real& PartialQuadratureFunction::getIntegrationPointValue(
      const size_type o) const {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline mgis::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type o) {
    return mgis::span<real>(this->values.data() + this->getDataOffset(o),
                            this->data_size);
  }  // end of getIntegrationPointValues

  inline mgis::span<const real>
  PartialQuadratureFunction::getIntegrationPointValues(
      const size_type o) const {
    return mgis::span<const real>(this->values.data() + this->getDataOffset(o),
                                  this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline mgis::span<real, N>
  PartialQuadratureFunction::getIntegrationPointValues(const size_type o) {
    return mgis::span<real, N>(this->values.data() + this->getDataOffset(o),
                               this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline mgis::span<const real, N>
  PartialQuadratureFunction::getIntegrationPointValues(
      const size_type o) const {
    return mgis::span<const real, N>(
        this->values.data() + this->getDataOffset(o), this->data_size);
  }  // end of getIntegrationPointValues


} // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX */
