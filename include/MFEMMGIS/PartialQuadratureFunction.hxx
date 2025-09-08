/*!
 * \file   PartialQuadratureFunction.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX

#include <span>
#include <limits>
#include <memory>
#include <vector>
#include <functional>

#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple data structure describing how the data of a partial
   * quadrature function is mapped in memory
   */
  struct PartialQuadratureFunctionDataLayout {
    // \brief default constructor
    PartialQuadratureFunctionDataLayout() = default;
    // \brief move constructor
    PartialQuadratureFunctionDataLayout(PartialQuadratureFunctionDataLayout&&) =
        default;
    // \brief copy constructor
    PartialQuadratureFunctionDataLayout(
        const PartialQuadratureFunctionDataLayout&) = default;
    // \brief move assignement
    PartialQuadratureFunctionDataLayout& operator=(
        PartialQuadratureFunctionDataLayout&&) = default;
    // \brief standard assignement
    PartialQuadratureFunctionDataLayout& operator=(
        const PartialQuadratureFunctionDataLayout&) = default;
    //! \return the number of components
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
     * \return the data offset associated with the given integration point.
     * \param[in] o: offset associated with the integration point
     */
    size_type getDataOffset(const size_type) const noexcept;
    //! \brief data stride
    size_type data_stride = size_type{};
    /*!
     * \brief begin of the data (offset of the function with respect to the
     * beginning of data values)
     */
    size_type data_begin = size_type{};
    //! \brief data size
    size_type data_size = size_type{};
  };  // end of struct PartialQuadratureFunctionDataLayout

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
  struct MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView
      : PartialQuadratureFunctionDataLayout {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: offset of data
     * \param[in] ds: size of the data per integration points
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::span<const real>,
        const size_type,
        const size_type);
    //! \return the underlying quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the underlying quadrature space
    std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    const real* data(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    const real* data(const size_type, const size_type) const;
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
    std::span<const real, N> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<const real> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<const real> getIntegrationPointValues(const size_type,
                                                    const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<const real> operator()(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<const real> operator()(const size_type, const size_type) const;
    //! \return a view to the function values
    std::span<const real> getValues() const;
    /*!
     * \return if the current function has the same quadrature space and the
     * same number of components than the given view
     * \param[in] v: view
     */
    bool checkCompatibility(
        const ImmutablePartialQuadratureFunctionView&) const;
    //! \brief destructor
    ~ImmutablePartialQuadratureFunctionView();

   protected:
    //! \brief default constructor
    ImmutablePartialQuadratureFunctionView();
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] ds: data size
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view (stride)
     */
    ImmutablePartialQuadratureFunctionView(
        std::shared_ptr<const PartialQuadratureSpace>,
        const size_type,
        const size_type,
        const size_type);
    //! \brief underlying finite element space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    //! \brief underlying values
    std::span<const real> immutable_values;
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
        std::span<real>,
        const size_type = 0,
        const size_type = std::numeric_limits<size_type>::max());
    /*!
     * \brief move constructor
     * \param[in] f: moved function
     * \param[in] local_copy: copy locally the function values if the moved
     * function does not holds them, i.e. is a view.
     * \note if the moved function holds the memory, the move constructor will
     * take ownership of the memory
     */
    PartialQuadratureFunction(PartialQuadratureFunction&&, const bool = false);
    //! \brief copy constructor
    PartialQuadratureFunction(const PartialQuadratureFunction&);
    /*!
     * \brief constructor from an immutable view
     *
     * \note: data are copied in a local array
     * \param[in] v: view
     */
    explicit PartialQuadratureFunction(
        const ImmutablePartialQuadratureFunctionView&);
    //! \brief assignement operator
    PartialQuadratureFunction& operator=(
        const ImmutablePartialQuadratureFunctionView&);
    //! \brief standard assignement operator
    PartialQuadratureFunction& operator=(const PartialQuadratureFunction&);
    //     //! \brief move assignement operator
    //     PartialQuadratureFunction& operator=(PartialQuadratureFunction&&);
    //
    using ImmutablePartialQuadratureFunctionView::data;
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValue;
    using ImmutablePartialQuadratureFunctionView::getIntegrationPointValues;
    using ImmutablePartialQuadratureFunctionView::getValues;
    using ImmutablePartialQuadratureFunctionView::operator();
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    real* data(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    real* data(const size_type, const size_type);
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
    std::span<real, N> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<real> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<real> getIntegrationPointValues(const size_type, const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<real> operator()(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<real> operator()(const size_type, const size_type);
    //! \return a view to the function values
    std::span<real> getValues();
    //! \brief destructor
    ~PartialQuadratureFunction();

   protected:
    /*!
     * \brief turns this function into a view to the given function
     * \param[in] f: function
     */
    void makeView(PartialQuadratureFunction&);
    /*!
     * \brief turns this function into a view to the given function
     * \param[in] f: function
     */
    void copy(const ImmutablePartialQuadratureFunctionView&);
    /*!
     * \brief copy values for an immutable view
     * \param[in] v: view
     * \note: no compatibility checks are perfomed
     */
    void copyValues(const ImmutablePartialQuadratureFunctionView&);
    //! \brief underlying values
    std::span<real> mutable_values;
    /*!
     * \brief storage for the values when the partial function holds the
     * values
     */
    std::vector<real> local_values_storage;
  };  // end of struct PartialQuadratureFunction

  /*!
   * \return a grid function able to store the result of the given functions
   * \param[in] p: underlying problem
   * \param[in] fcts: functions
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  std::pair<std::unique_ptr<FiniteElementSpace<parallel>>,
            std::unique_ptr<GridFunction<parallel>>>
  makeGridFunction(const std::vector<ImmutablePartialQuadratureFunctionView>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<true>>,
                             std::unique_ptr<GridFunction<true>>>
  makeGridFunction<true>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<false>>,
                             std::unique_ptr<GridFunction<false>>>
  makeGridFunction<false>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  /*!
   * \return a grid function able to store the result of the given functions
   * \param[in] p: underlying problem
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  std::pair<std::unique_ptr<FiniteElementSpace<parallel>>,
            std::unique_ptr<GridFunction<parallel>>>
  makeGridFunction(const std::vector<ImmutablePartialQuadratureFunctionView>&,
                   const Mesh<parallel>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<true>>,
                             std::unique_ptr<GridFunction<true>>>
  makeGridFunction<true>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<true>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<false>>,
                             std::unique_ptr<GridFunction<false>>>
  makeGridFunction<false>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<false>&);

  /*!
   * \return a grid function able to store the result of the given functions
   * \param[in] p: underlying problem
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the values of the grid function are computed by the
   * updateGridFunction function.
   */
  template <bool parallel>
  std::pair<std::unique_ptr<FiniteElementSpace<parallel>>,
            std::unique_ptr<GridFunction<parallel>>>
  makeGridFunction(const std::vector<ImmutablePartialQuadratureFunctionView>&,
                   const SubMesh<parallel>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<true>>,
                             std::unique_ptr<GridFunction<true>>>
  makeGridFunction<true>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<true>&);

  template <>
  MFEM_MGIS_EXPORT std::pair<std::unique_ptr<FiniteElementSpace<false>>,
                             std::unique_ptr<GridFunction<false>>>
  makeGridFunction<false>(
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<false>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[in] f: function
   * \param[in] fcts: functions
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[in] f: function
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<parallel>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<true>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const Mesh<false>&);

  /*!
   * \brief update a grid function using the values of the given functions
   * \param[in] f: function
   * \param[in] fcts: functions
   * \param[in] mesh: submesh on which the grid function is defined
   * \note the grid function must have been created by `makeGridFunction`
   */
  template <bool parallel>
  void updateGridFunction(
      GridFunction<parallel>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<parallel>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<true>&);

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const SubMesh<false>&);

}  // namespace mfem_mgis

#ifdef MGIS_FUNCTION_SUPPORT

namespace mfem_mgis {

  constexpr bool check(AbstractErrorHandler&,
                       const ImmutablePartialQuadratureFunctionView&) noexcept;

  constexpr void allocateWorkspace(
      ImmutablePartialQuadratureFunctionView&) noexcept;

  mgis::size_type getNumberOfComponents(
      const ImmutablePartialQuadratureFunctionView&) noexcept;

  MFEM_MGIS_EXPORT const PartialQuadratureSpace& getSpace(
      const ImmutablePartialQuadratureFunctionView&);

  MFEM_MGIS_EXPORT const PartialQuadratureSpace& getSpace(
      const PartialQuadratureFunction&);

  /*!
   * \brief this function is deleted to avoid
   * `PartialQuadratureFunction` to match `mgis::EvaluatorConcept`
   * as it is not a lightweight object
   */
  void allocateWorkspace(PartialQuadratureFunction&) = delete;

}  // namespace mfem_mgis

namespace mgis::function {

  static_assert(
      EvaluatorConcept<mfem_mgis::ImmutablePartialQuadratureFunctionView>);
  static_assert(FunctionConcept<mfem_mgis::PartialQuadratureFunction>);
  static_assert(!EvaluatorConcept<mfem_mgis::PartialQuadratureFunction>);

}  // end of namespace mgis::function

#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/PartialQuadratureFunction.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX */
