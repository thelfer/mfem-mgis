/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluator.hxx
 * \brief  This file declares the evaluators of the rotation matrix and of the
 * rotated gradients and thermodynamic forces
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATOR_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATOR_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Utilities/Buffer.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  //! \brief an evaluator returning the rotation matrix
  struct RotationMatrixQPEvaluator {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotationMatrixQPEvaluator(const Material& m);
    /*!
     * \brief perform consistency checks
     * \param[in, out] eh: error handler
     * \return true on success
     */
    [[nodiscard]] bool check(AbstractErrorHandler& eh) const;
    //! \return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    /*!
     * \brief access operator
     * \param[in] i: integration point index
     * \return the rotation matrix at the integration point
     */
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying material
    const Material& material;
  };  // end of RotationMatrixQPEvaluator

  /*!
   * \brief return the quadrature space of an evaluator
   * \param[in] e: evaluator
   * \return the quadrature space
   */
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotationMatrixQPEvaluator& e);

  /*!
   * \brief perform consistency checks
   * \param[in, out] eh: error handler
   * \param[in] e: evaluator
   * \return true on success
   */
  [[nodiscard]] bool check(AbstractErrorHandler& eh,
                           const RotationMatrixQPEvaluator& e);

  /*!
   * \brief return the number of components
   * \param[in] e: evaluator
   * \return the number of components
   */
  [[nodiscard]] constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixQPEvaluator& e) noexcept;

  /*!
   * \brief an evaluator returning the rotated thermodynamic forces
   */
  template <size_type ThermodynamicForcesSize = dynamic_extent>
  struct RotatedThermodynamicForcesMatrixQPEvaluator {
    /*!
     * \brief constructor
     * \param[in] m: material
     * \param[in] s: state considered
     */
    RotatedThermodynamicForcesMatrixQPEvaluator(
        const Material& m,
        const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
    //! \return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    /*!
     * \brief perform consistency checks
     * \param[in, out] ctx: error handler
     * \return true on success
     */
    [[nodiscard]] bool check(AbstractErrorHandler& ctx) const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief access operator
     * \param[in] i: integration point index
     * \return the rotated thermodynamic forces at the integration point
     */
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief thermodynamic forces
    std::span<real> thforces;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    mutable Buffer<ThermodynamicForcesSize> buffer;
  };  // end of
      // RotatedThermodynamicForcesMatrixQPEvaluator

  /*!
   * \brief return the quadrature space of an evaluator
   * \param[in] e: evaluator
   * \return the quadrature space
   */
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e);
  /*!
   * \brief perform consistency checks
   * \param[in, out] eh: error handler
   * \param[in] e: evaluator
   * \return true on success
   */
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] bool check(AbstractErrorHandler& eh,
                           const RotatedThermodynamicForcesMatrixQPEvaluator<
                               ThermodynamicForcesSize>& e);
  /*!
   * \brief return the number of components
   * \param[in] e: evaluator
   * \return the number of components
   */
  template <size_type ThermodynamicForcesSize>
  mgis::size_type getNumberOfComponents(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e) noexcept;

  /*!
   * \brief an evaluator returning the gradients rotated in the global frame
   */
  template <size_type GradientsSize = dynamic_extent>
  struct RotatedGradientsMatrixQPEvaluator {
    /*!
     * \brief constructor
     * \param[in] m: material
     * \param[in] s: state considered
     */
    RotatedGradientsMatrixQPEvaluator(
        const Material& m,
        const Material::StateSelection s = Material::END_OF_TIME_STEP);
    /*!
     * \brief perform consistency checks
     * \param[in, out] ctx: error handler
     * \return true on success
     */
    [[nodiscard]] bool check(AbstractErrorHandler& ctx) const;
    //! \return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief access operator
     * \param[in] i: integration point index
     * \return the rotated gradients at the integration point
     */
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief gradients
    std::span<real> gradients;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    mutable Buffer<GradientsSize> buffer;
  };  // end of RotatedGradientsMatrixQPEvaluator

  /*!
   * \brief return the quadrature space of an evaluator
   * \param[in] e: evaluator
   * \return the quadrature space
   */
  template <size_type GradientsSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e);
  /*!
   * \brief perform consistency checks
   * \param[in, out] eh: error handler
   * \param[in] e: evaluator
   * \return true on success
   */
  template <size_type GradientsSize>
  [[nodiscard]] bool check(
      AbstractErrorHandler& eh,
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e);
  /*!
   * \brief return the number of components
   * \param[in] e: evaluator
   * \return the number of components
   */
  template <size_type GradientsSize>
  mgis::size_type getNumberOfComponents(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e) noexcept;

  /*!
   * \brief check if the given evaluators have the same partial quadrature space
   * \param[in, out] ctx: execution context
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   * \return true on success
   */
  template <QPEvaluatorConcept EvaluatorType1,
            QPEvaluatorConcept EvaluatorType2>
  bool checkMatchingQuadratureSpaces(Context& ctx,
                                     const EvaluatorType1& e1,
                                     const EvaluatorType2& e2) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/QPEvaluator/QPEvaluator.ixx"

namespace mfem_mgis {

  static_assert(mgis::function::EvaluatorConcept<RotationMatrixQPEvaluator>);
  static_assert(!mgis::function::FunctionConcept<RotationMatrixQPEvaluator>);
  static_assert(mgis::function::EvaluatorConcept<
                RotatedThermodynamicForcesMatrixQPEvaluator<>>);
  static_assert(!mgis::function::FunctionConcept<
                RotatedThermodynamicForcesMatrixQPEvaluator<>>);
  static_assert(
      mgis::function::EvaluatorConcept<RotatedGradientsMatrixQPEvaluator<>>);
  static_assert(
      !mgis::function::FunctionConcept<RotatedGradientsMatrixQPEvaluator<>>);

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATOR_HXX */
