/*!
 * \file   PartialQuadratureFunction.hxx
 * \brief  This file declares the `PartialQuadratureFunction` class and its
 * views
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX

#include <span>
#include <limits>
#include <memory>
#include <vector>
#include <optional>
#include <functional>

#include "MGIS/StorageMode.hxx"
#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  using mgis::StorageMode;

  /*!
   * \brief data structure used to specify a view from external data
   *
   * To be valid, the following conditions must hold:
   *
   * - data_begin must be non-negative
   * - data_size and data_stride must be strictly positive
   * - data_begin + data_size <= data_stride
   *
   */
  struct ViewSpecifications {
    /*!
     * \brief begin of the data (offset with respect to the
     * beginning of data values)
     */
    size_type data_begin;
    //! \brief data size
    size_type data_size;
    //! \brief data stride
    size_type data_stride;
  };

  /*!
   * \brief a simple data structure describing how the data of a partial
   * quadrature function is mapped in memory
   */
  struct PartialQuadratureFunctionDataLayout : ViewSpecifications {
    //! \brief default constructor
    PartialQuadratureFunctionDataLayout() = default;
    //! \brief move constructor
    PartialQuadratureFunctionDataLayout(PartialQuadratureFunctionDataLayout&&) =
        default;
    //! \brief copy constructor
    PartialQuadratureFunctionDataLayout(
        const PartialQuadratureFunctionDataLayout&) = default;
    //! \brief move assignment
    PartialQuadratureFunctionDataLayout& operator=(
        PartialQuadratureFunctionDataLayout&&) = default;
    //! \brief standard assignment
    PartialQuadratureFunctionDataLayout& operator=(
        const PartialQuadratureFunctionDataLayout&) = default;
    //! \return if the function is scalar
    bool isScalar() const noexcept;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;

    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    size_type getDataStride() const noexcept;
    //! \return the offset of the first element
    size_type getDataOffset() const noexcept;
    //! \brief destructor
    ~PartialQuadratureFunctionDataLayout() = default;

   protected:
    /*!
     * \brief compute the data offset of an integration point
     * \return the data offset associated with the given integration point.
     * \param[in] o: offset associated with the integration point
     */
    size_type getDataOffset(const size_type o) const noexcept;
  };  // end of struct PartialQuadratureFunctionDataLayout

  /*!
   * \brief immutable view of a partial quadrature function
   *
   * The `ImmutablePartialQuadratureFunctionView` defines an immutable view
   * associated with a partial quadrature function on a memory region.
   *
   * This memory region may contain more data than the one associated with the
   * quadrature function as illustrated by the following figure:
   *
   * \verbatim
   * |---------------------------------------------------------------|
   * <-                         Raw data                            ->
   * |---------------------------------------|
   * <- Data of the first integration point-->
   * |      |---------------|                |
   *        <-function data->
   *        ^                                ^
   *        |                                |
   *    data_begin                           |
   *                                    data_stride
   * \endverbatim
   *
   * The size of all the data (including the one not related to the partial
   * quadrature function) associated with one integration point is called the
   * `data_stride` in the `ImmutablePartialQuadratureFunctionView` class.
   *
   * Inside the data associated with one integration point, the function data
   * starts at the offset given by `data_begin`.
   *
   * The size of the data held by the function per integration point, i.e. the
   * number of components of the function is given by `data_size`.
   */
  struct MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView
      : PartialQuadratureFunctionDataLayout {
    /*!
     * \brief constructor from already allocated data
     *
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] specs: specifications
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::span<const real> v,
        const ViewSpecifications& specs);
    //! \brief move constructor
    ImmutablePartialQuadratureFunctionView(
        ImmutablePartialQuadratureFunctionView&&) noexcept;
    //! \brief copy constructor
    ImmutablePartialQuadratureFunctionView(
        const ImmutablePartialQuadratureFunctionView&) noexcept;
    //! \brief move assignment
    ImmutablePartialQuadratureFunctionView& operator=(
        ImmutablePartialQuadratureFunctionView&&) noexcept;
    //! \brief copy assignment
    ImmutablePartialQuadratureFunctionView& operator=(
        const ImmutablePartialQuadratureFunctionView&) noexcept;
    //! \return the underlying quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the underlying quadrature space
    std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return a pointer to the data
     */
    const real* data(const size_type o) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return a pointer to the data
     */
    const real* data(const size_type e, const size_type i) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the value
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */

    const real& getIntegrationPointValue(const size_type o) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the value
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type e,
                                         const size_type i) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    template <size_type N>
    std::span<const real, N> getIntegrationPointValues(const size_type o) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    std::span<const real> getIntegrationPointValues(const size_type o) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the values
     */
    std::span<const real> getIntegrationPointValues(const size_type e,
                                                    const size_type i) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    std::span<const real> operator()(const size_type o) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the values
     */
    std::span<const real> operator()(const size_type e,
                                     const size_type i) const;
    //! \return a view to the function values
    std::span<const real> getValues() const;
    /*!
     * \brief check the compatibility with the given view
     * \return if the current function has the same quadrature space and the
     * same number of components as the given view
     * \param[in] v: view
     */
    bool checkCompatibility(
        const ImmutablePartialQuadratureFunctionView& v) const;
    //! \brief destructor
    ~ImmutablePartialQuadratureFunctionView();

   protected:
    //! \brief default constructor
    ImmutablePartialQuadratureFunctionView();
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] specs: specifications
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace> s,
        const ViewSpecifications& specs);
    //! \brief underlying partial quadrature space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    //! \brief underlying values
    std::span<const real> immutable_values;
  };  // end of ImmutablePartialQuadratureFunctionView

  //! \brief mutable view of a partial quadrature function
  struct MFEM_MGIS_EXPORT PartialQuadratureFunctionView
      : ImmutablePartialQuadratureFunctionView {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] specs: specifications
     */
    PartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::span<real> v,
        const ViewSpecifications& specs);
    //! \brief move constructor
    PartialQuadratureFunctionView(PartialQuadratureFunctionView&&) noexcept;
    //! \brief copy constructor
    PartialQuadratureFunctionView(
        const PartialQuadratureFunctionView&) noexcept;
    //! \brief move assignment
    PartialQuadratureFunctionView& operator=(
        PartialQuadratureFunctionView&&) noexcept;
    //! \brief copy assignment
    PartialQuadratureFunctionView& operator=(
        const PartialQuadratureFunctionView&) noexcept;
    //
    using ImmutablePartialQuadratureFunctionView::data;
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValue;
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValues;
    using ImmutablePartialQuadratureFunctionView::getValues;
    using ImmutablePartialQuadratureFunctionView::operator();
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return a pointer to the data
     */
    real* data(const size_type o);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return a pointer to the data
     */
    real* data(const size_type e, const size_type i);
    /*!
     * \brief return the value associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the value
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type o);
    /*!
     * \brief return the value associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the value
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type e, const size_type i);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    template <size_type N>
    std::span<real, N> getIntegrationPointValues(const size_type o);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    std::span<real> getIntegrationPointValues(const size_type o);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the values
     */
    std::span<real> getIntegrationPointValues(const size_type e,
                                              const size_type i);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \return the values
     */
    std::span<real> operator()(const size_type o);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \return the values
     */
    std::span<real> operator()(const size_type e, const size_type i);
    //! \return a view to the function values
    std::span<real> getValues();

   protected:
    //! \brief default constructor
    PartialQuadratureFunctionView();
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] specs: specifications
     */
    PartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace> s,
        const ViewSpecifications& specs);
    //! \brief underlying values
    std::span<real> mutable_values;
  };  // end of PartialQuadratureFunctionView

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   *
   * A partial quadrature function is movable, but not copyable.
   * Most of the time, the partial quadrature function holds the memory,
   * but the `borrow` method allows to use externally allocated memory
   *
   * The main reason for this choice is to force usage of parallel algorithms.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureFunction
      : PartialQuadratureFunctionView {
    /*!
     * \brief evaluate a partial quadrature function at each integration point
     * \note if required, the current integration point can be retrieved using
     * the `GetIntPoint` of the element transformation
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     * \return the partial quadrature function
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::function<real(const mfem::FiniteElement&,
                           mfem::ElementTransformation&)> f);
    /*!
     * \brief evaluate a spatial function in 2D
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     * \return the partial quadrature function
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::function<real(real, real)> f);
    /*!
     * \brief evaluate a spatial function in 3D
     * \param[in] s: partial quadrature space
     * \param[in] f: function to be evaluated
     * \return the partial quadrature function
     */
    static std::shared_ptr<PartialQuadratureFunction> evaluate(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::function<real(real, real, real)> f);
    /*!
     * \brief copy a view
     * \return a newly created `PartialQuadratureFunction` which contains a copy
     * of the values of the source
     *
     * \param[in, out] ctx: execution context
     * \param[in] v: view being copied
     *
     * \note this method has been introduced to avoid defining a copy
     * constructor and an assignment operator in the
     * `PartialQuadratureFunction` class
     *
     * \note we strongly recommend not using this function, as it uses
     * `std::copy` to copy the values. It is much better to create another
     * version relying on the `assign` algorithm using the parallel programming
     * model you wish to use
     */
    [[nodiscard]] static std::optional<PartialQuadratureFunction> copy(
        Context& ctx, const ImmutablePartialQuadratureFunctionView& v) noexcept;
    /*!
     * \brief create a partial quadrature function on external memory
     * \return a partial quadrature function that does not manage its
     * values, but borrows them from an external memory
     *
     * \param[in, out] ctx: execution context
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] specs: specifications
     *
     * \pre v.size() must be equal to stride * getSpaceSize(*s)
     */
    [[nodiscard]] static std::optional<PartialQuadratureFunction> borrow(
        Context& ctx,
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::span<real> v,
        const ViewSpecifications& specs) noexcept;
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] nv: size of the data stored per integration point.
     */
    PartialQuadratureFunction(std::shared_ptr<const PartialQuadratureSpace> s,
                              const size_type nv = 1);
    /*!
     * \brief move constructor
     * \param[in, out] f: moved function
     * \note if the moved function holds the memory, the move constructor will
     * take ownership of the memory
     */
    PartialQuadratureFunction(PartialQuadratureFunction&& f);
    //! \return a view of the function
    PartialQuadratureFunctionView view();
    //! \return an immutable view of the function
    ImmutablePartialQuadratureFunctionView view() const;

    //! \brief destructor
    ~PartialQuadratureFunction();

   protected:
    /*!
     * \brief constructor
     * \param[in] v: view to be copied
     */
    PartialQuadratureFunction(const ImmutablePartialQuadratureFunctionView& v);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] sm: storage mode
     * \param[in] v: values
     * \param[in] specs: specifications
     */
    PartialQuadratureFunction(std::shared_ptr<const PartialQuadratureSpace> s,
                              const StorageMode sm,
                              std::span<real> v,
                              const ViewSpecifications& specs);
    /*!
     * \brief turns this function into a view to the given function
     * \param[in] f: function
     */
    void makeView(PartialQuadratureFunction& f);
    /*!
     * \brief copy the given function
     * \param[in] f: function
     */
    void copy(const ImmutablePartialQuadratureFunctionView& f);
    /*!
     * \brief copy values from an immutable view
     * \param[in] v: view
     * \note the execution is aborted if the view is not compatible
     */
    void copyValues(const ImmutablePartialQuadratureFunctionView& v);
    /*!
     * \brief storage for the values when the partial function holds the
     * values
     */
    std::vector<real> local_values_storage;
  };  // end of struct PartialQuadratureFunction

  /*!
   * \brief assign the values of an immutable view to the values of a mutable
   * view.
   *
   * \param[in, out] ctx: execution context
   * \param[in] f: mutable view
   * \param[in] v: immutable view
   * \return true on success, false on failure
   *
   * \pre both views must have the same number of components.
   * \note partial quadrature spaces of the views may be
   * different, we only require them to have the same identifier and the same
   * size.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool assign_values(
      Context& ctx,
      PartialQuadratureFunctionView f,
      const ImmutablePartialQuadratureFunctionView& v) noexcept;

  /*!
   * \brief update the partial quadrature function from the given grid function
   * \param[in, out] ctx: execution context
   * \param[out] dest: partial quadrature function to be updated
   * \param[in] src: grid function
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool update(
      Context& ctx,
      PartialQuadratureFunctionView& dest,
      const GridFunction<true>& src) noexcept;
  /*!
   * \brief update the partial quadrature function from the given grid function
   * \param[in, out] ctx: execution context
   * \param[out] dest: partial quadrature function to be updated
   * \param[in] src: grid function
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool update(
      Context& ctx,
      PartialQuadratureFunctionView& dest,
      const GridFunction<false>& src) noexcept;
  /*!
   * \brief create a grid function for the given functions
   * \return a grid function able to store the result of the given functions
   * \param[in, out] ctx: execution context
   * \param[in] fcts: functions
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  [[nodiscard]] std::unique_ptr<GridFunction<parallel>> makeGridFunction(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::unique_ptr<GridFunction<true>>
  makeGridFunction<true>(
      Context&, const std::vector<ImmutablePartialQuadratureFunctionView>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::unique_ptr<GridFunction<false>>
  makeGridFunction<false>(
      Context&, const std::vector<ImmutablePartialQuadratureFunctionView>&);

  /*!
   * \brief create a grid function for the given functions on the given mesh
   * \return a grid function able to store the result of the given functions
   * \param[in, out] ctx: execution context
   * \param[in] fcts: functions
   * \param[in] mesh: mesh on which the grid function is defined
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  std::unique_ptr<GridFunction<parallel>> makeGridFunction(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<parallel>& mesh);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<true>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<false>&);

  /*!
   * \brief create a grid function for the given functions on the given submesh
   * \return a grid function able to store the result of the given functions
   * \param[in, out] ctx: execution context
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  std::unique_ptr<GridFunction<parallel>> makeGridFunction(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<parallel>& mesh);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<true>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<false>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[out] f: grid function
   * \param[in] fcts: functions
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[out] f: grid function
   * \param[in] fcts: functions
   * \param[in] mesh: mesh on which the grid function is defined
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<parallel>& mesh);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<true>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<false>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[out] f: grid function
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<parallel>& mesh);

  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<true>&);

  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<false>&);

}  // namespace mfem_mgis

#ifdef MGIS_FUNCTION_SUPPORT

namespace mfem_mgis {

  //! \brief concept of an evaluator defined on a partial quadrature space
  template <typename EvaluatorType>
  concept QPEvaluatorConcept =
      ((mgis::function::EvaluatorConcept<EvaluatorType>)&&  //
       (requires(const EvaluatorType& e) {
         { getSpace(e) } -> std::same_as<const PartialQuadratureSpace&>;
       }));

  /*!
   * \brief perform consistency checks
   * \param[in, out] eh: error handler
   * \param[in] f: function
   * \return true on success
   */
  constexpr bool check(
      AbstractErrorHandler& eh,
      const ImmutablePartialQuadratureFunctionView& f) noexcept;

  //! \brief deleted, a partial quadrature function is not an evaluator
  bool check(AbstractErrorHandler&,
             const PartialQuadratureFunction&) noexcept = delete;

  /*!
   * \brief return the number of components
   * \param[in] f: function
   * \return the number of components
   */
  mgis::size_type getNumberOfComponents(
      const ImmutablePartialQuadratureFunctionView& f) noexcept;

  /*!
   * \brief return the quadrature space of a view
   * \param[in] f: view
   * \return the quadrature space
   */
  MFEM_MGIS_EXPORT const PartialQuadratureSpace& getSpace(
      const ImmutablePartialQuadratureFunctionView& f);

  /*!
   * \brief return the quadrature space of a function
   * \param[in] f: function
   * \return the quadrature space
   */
  MFEM_MGIS_EXPORT const PartialQuadratureSpace& getSpace(
      const PartialQuadratureFunction& f);

  /*!
   * \brief return a view of a function
   * \param[in] f: function
   * \return a view of the given function
   */
  PartialQuadratureFunctionView view(PartialQuadratureFunction& f);

  /*!
   * \brief return an immutable view of a function
   * \param[in] f: function
   * \return an immutable view of the given function
   */
  ImmutablePartialQuadratureFunctionView view(
      const PartialQuadratureFunction& f);

}  // namespace mfem_mgis

namespace mgis::function {

  template <>
  struct LightweightViewTraits<mfem_mgis::PartialQuadratureFunctionView>
      : std::true_type {};

  static_assert(
      EvaluatorConcept<mfem_mgis::ImmutablePartialQuadratureFunctionView>);
  static_assert(FunctionConcept<mfem_mgis::PartialQuadratureFunction>);
  static_assert(!EvaluatorConcept<mfem_mgis::PartialQuadratureFunction>);

}  // end of namespace mgis::function

#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/PartialQuadratureFunction.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX */
