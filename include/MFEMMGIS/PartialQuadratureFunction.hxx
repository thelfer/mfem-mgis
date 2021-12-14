/*!
 * \file   PartialQuadratureFunction.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX

#include <memory>
#include <vector>
#include <limits>
#include <functional>
#include "MGIS/Span.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   *
   * The `ImmutablePartialQuadratureFunctionView` defines an immutable view
   * associated with a partial quadrature function on a memory region.
   *
   * This memory region may contain more data than the one associated with the
   * quadrature function as illustrated by the following figure:
   *
   * |---------------------------------------------------------------|
   * <-                         Raw data                            ->
   * |---------------------------------------|
   * <- Data of the first integration point-->
   * |      |---------------|                |
   *        <-function data->
   *        ^                                ^
   *        |                                |
   *    data_begin                           |
   *                                     data_size
   *
   * The size of the all data (including the one not related to the partial
   * quadrature function) associated with one integration point is called the
   * `data_stride` in the `ImmutablePartialQuadratureFunctionView` class.
   *
   * Inside the data associated with on integration point, the function data
   * starts at the offset given by `data_begin`.
   *
   * The size of the data hold by the function per integration point, i.e. th
   * number of components of the function is given by `data_size`.
   */
  struct MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: offset of data
     * \param[in] ds: size of the data per integration points
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace>,
        mgis::span<const real>,
        const size_type = 0,
        const size_type = std::numeric_limits<size_type>::max());
    //! \return the underlying quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type,
                                         const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    mgis::span<const real, N> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    mgis::span<const real> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<const real> getIntegrationPointValues(const size_type,
                                                     const size_type) const;
    //! \return a view to the function values
    mgis::span<const real> getValues() const;
    //! \return the number of components
    size_type getNumberOfComponents() const;
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    size_type getDataStride() const;
    //! \return the offset of the first element
    size_type getInitialDataOffset() const;
    //! \brief destructor
    ~ImmutablePartialQuadratureFunctionView();

   protected:
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] ds: data size
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view (stride)
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace>,
        const size_type ds,
        const size_type,
        const size_type);
    //! \brief underlying finite element space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] o: offset associated with the integration point
     */
    size_type getDataOffset(const size_type) const;
    //! \brief underlying values
    mgis::span<const real> immutable_values;
    //! \brief data stride
    size_type data_stride;
    /*!
     * \brief begin of the data (offset of the function with respect to the
     * beginning of data values)
     */
    size_type data_begin;
    //! \brief data size
    size_type data_size;
  };  // end of ImmutablePartialQuadratureFunctionView

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureFunction
      : ImmutablePartialQuadratureFunctionView {
    /*!
     * \brief evaluate a partial quadrature function at each integration point
     * \note if required, the current integration point can be retrieved using
     * the `GetIntPoint` of the element transformation
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::function<real(const mfem::FiniteElement&,
                           mfem::ElementTransformation&)>);
    /*!
     * \brief evaluate a spatial function in 2D
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::function<real(real, real)>);
    /*!
     * \brief evaluate a spatial function in 3D
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::function<real(real, real, real)>);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] size: size of the data stored per integration points.
     */
    PartialQuadratureFunction(std::shared_ptr<const PartialQuadratureSpace>,
                              const size_type = 1);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view
     */
    PartialQuadratureFunction(
        std::shared_ptr<const PartialQuadratureSpace>,
        mgis::span<real>,
        const size_type = 0,
        const size_type = std::numeric_limits<size_type>::max());
    //
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValue;
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValues;
    using ImmutablePartialQuadratureFunctionView::getValues;
    /*!
     * \brief return the value associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type);
    /*!
     * \brief return the value associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type, const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    mgis::span<real, N> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    mgis::span<real> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<real> getIntegrationPointValues(const size_type,
                                               const size_type);
    //! \return a view to the function values
    mgis::span<real> getValues();
    //! \brief destructor
    ~PartialQuadratureFunction();

   protected:
    //! \brief underlying values
    mgis::span<real> values;
    /*!
     * \brief storage for the values when the partial function holds the
     * values
     */
    std::vector<real> local_values_storage;
  };  // end of struct PartialQuadratureFunction

}  // namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureFunction.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX */
