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

  inline PartialQuadratureFunctionView::PartialQuadratureFunctionView() =
      default;

  inline PartialQuadratureFunctionView::PartialQuadratureFunctionView(
      std::shared_ptr<const PartialQuadratureSpace> s,
      const size_type ds,
      const size_type db,
      const size_type dsize)
      : ImmutablePartialQuadratureFunctionView(s, ds, db, dsize) {
  }  // end of PartialQuadratureFunctionView

  inline PartialQuadratureFunctionView::PartialQuadratureFunctionView(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::span<real> v,
      const size_type db,
      const size_type ds)
      : ImmutablePartialQuadratureFunctionView(s, v, db, ds),
        mutable_values(v) {}  // end of PartialQuadratureFunctionView

  inline real* PartialQuadratureFunctionView::data(const size_type o) {
    return this->mutable_values.data() + this->getDataOffset(o);
  }  // end of getIntegrationPointValues

  inline real& PartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type o) {
    return *(this->mutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<real>
  PartialQuadratureFunctionView::getIntegrationPointValues(const size_type o) {
    return std::span<real>(this->mutable_values.data() + this->getDataOffset(o),
                           this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<real, N>
  PartialQuadratureFunctionView::getIntegrationPointValues(const size_type o) {
    return std::span<real, N>(
        this->mutable_values.data() + this->getDataOffset(o), this->data_size);
  }  // end of getIntegrationPointValues

  inline std::span<real> PartialQuadratureFunctionView::operator()(
      const size_type o) {
    return this->getIntegrationPointValues(o);
  }

  inline std::span<real> PartialQuadratureFunctionView::operator()(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(e, i);
  }

  inline std::span<real> PartialQuadratureFunctionView::getValues() {
    return this->mutable_values;
  }

  inline PartialQuadratureFunctionView PartialQuadratureFunction::view() {
    return *this;
  }  // end of view

  inline ImmutablePartialQuadratureFunctionView
  PartialQuadratureFunction::view() const {
    return *this;
  }  // end of view

#ifdef MGIS_FUNCTION_SUPPORT

  constexpr bool check(AbstractErrorHandler&,
                       const ImmutablePartialQuadratureFunctionView&) noexcept {
    return true;
  }  // end of check

  inline mgis::size_type getNumberOfComponents(
      const ImmutablePartialQuadratureFunctionView& f) noexcept {
    return static_cast<mgis::size_type>(f.getNumberOfComponents());
  }  // end of getNumberOfComponents

  inline PartialQuadratureFunctionView view(PartialQuadratureFunction& f) {
    return f.view();
  }  // end of view

  inline ImmutablePartialQuadratureFunctionView view(
      const PartialQuadratureFunction& f) {
    return f.view();
  }  // end of view

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_IXX */
