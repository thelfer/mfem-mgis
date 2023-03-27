/*!
 * \file   include/MFEMMGIS/PartialQuadratureFunction.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX

namespace mfem_mgis {

  inline size_type PartialQuadratureFunctionDataLayout::getDataStride() const {
    return this->data_stride;
  }  // end of getDataStride

  inline size_type PartialQuadratureFunctionDataLayout::getDataOffset() const {
    return this->data_begin;
  }  // end of getDataOffset

  inline size_type PartialQuadratureFunctionDataLayout::getNumberOfComponents()
      const {
    return this->data_size;
  }  // end of PartialQuadratureFunctionDataLayout::getNumberOfComponents

  inline size_type PartialQuadratureFunctionDataLayout::getDataOffset(
      const size_type o) const {
    return o * (this->data_stride) + this->data_begin;
  }  // end of getDataOffset

  inline const PartialQuadratureSpace&
  ImmutablePartialQuadratureFunctionView::getPartialQuadratureSpace() const {
    return *(this->qspace);
  }  // end of getPartialQuadratureSpace

  inline std::shared_ptr<const PartialQuadratureSpace>
  ImmutablePartialQuadratureFunctionView::getPartialQuadratureSpacePointer()
      const {
    return this->qspace;
  }  // end of getPartialQuadratureSpacePointer

  inline mgis::span<const real>
  ImmutablePartialQuadratureFunctionView::getValues() const {
    return this->immutable_values;
  }

  inline const real&
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type o) const {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline mgis::span<const real>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return mgis::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline mgis::span<const real, N>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return mgis::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  inline real& PartialQuadratureFunction::getIntegrationPointValue(
      const size_type o) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline mgis::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type o) {
    return mgis::span<real>(this->values.data() + this->getDataOffset(o),
                            this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline mgis::span<real, N>
  PartialQuadratureFunction::getIntegrationPointValues(const size_type o) {
    return mgis::span<real, N>(this->values.data() + this->getDataOffset(o),
                               this->data_size);
  }  // end of getIntegrationPointValues

  inline mgis::span<real> PartialQuadratureFunction::getValues() {
    return this->values;
  }

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX */
