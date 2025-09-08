/*!
 * \file   include/MFEMMGIS/PartialQuadratureFunction.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX

namespace mfem_mgis {

  inline bool PartialQuadratureFunctionDataLayout::isScalar() const noexcept {
    return this->getNumberOfComponents() == 1;
  }  // end of isScalar

  inline size_type PartialQuadratureFunctionDataLayout::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  inline size_type PartialQuadratureFunctionDataLayout::getDataOffset()
      const noexcept {
    return this->data_begin;
  }  // end of getDataOffset

  inline size_type PartialQuadratureFunctionDataLayout::getNumberOfComponents()
      const noexcept {
    return this->data_size;
  }  // end of PartialQuadratureFunctionDataLayout::getNumberOfComponents

  inline size_type PartialQuadratureFunctionDataLayout::getDataOffset(
      const size_type o) const noexcept {
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

  inline std::span<const real>
  ImmutablePartialQuadratureFunctionView::getValues() const {
    return this->immutable_values;
  }

  inline const real&
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type o) const {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline const real* ImmutablePartialQuadratureFunctionView::data(
      const size_type o) const {
    return this->immutable_values.data() + this->getDataOffset(o);
  }  // end of getIntegrationPointValues

  inline std::span<const real>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return std::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<const real, N>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return std::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  inline std::span<const real>
  ImmutablePartialQuadratureFunctionView::operator()(const size_type o) const {
    return this->getIntegrationPointValues(o);
  }

  inline std::span<const real>
  ImmutablePartialQuadratureFunctionView::operator()(const size_type e,
                                                     const size_type i) const {
    return this->getIntegrationPointValues(e, i);
  }

  inline real* PartialQuadratureFunction::data(const size_type o) {
    return this->mutable_values.data() + this->getDataOffset(o);
  }  // end of getIntegrationPointValues

  inline real& PartialQuadratureFunction::getIntegrationPointValue(
      const size_type o) {
    return *(this->mutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type o) {
    return std::span<real>(this->mutable_values.data() + this->getDataOffset(o),
                           this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<real, N>
  PartialQuadratureFunction::getIntegrationPointValues(const size_type o) {
    return std::span<real, N>(
        this->mutable_values.data() + this->getDataOffset(o), this->data_size);
  }  // end of getIntegrationPointValues

  inline std::span<real> PartialQuadratureFunction::operator()(
      const size_type o) {
    return this->getIntegrationPointValues(o);
  }

  inline std::span<real> PartialQuadratureFunction::operator()(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(e, i);
  }

  inline std::span<real> PartialQuadratureFunction::getValues() {
    return this->mutable_values;
  }

#ifdef MGIS_FUNCTION_SUPPORT

  constexpr bool check(AbstractErrorHandler&,
                       const ImmutablePartialQuadratureFunctionView&) noexcept {
    return true;
  }  // end of check

  constexpr void allocateWorkspace(
      ImmutablePartialQuadratureFunctionView&) noexcept {
  }  // end of allocateWorkspace

  inline mgis::size_type getNumberOfComponents(
      const ImmutablePartialQuadratureFunctionView& f) noexcept {
    return static_cast<mgis::size_type>(f.getNumberOfComponents());
  }  // end of getNumberOfComponents

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX */
