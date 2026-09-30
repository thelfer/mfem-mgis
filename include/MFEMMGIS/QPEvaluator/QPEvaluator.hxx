/*!
 * \file   MFEMMGIS/QPEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATOR_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATOR_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Utilities/Buffer.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
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
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler& eh) const;
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    // access operator
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying material
    const Material& material;
  };  // end of RotationMatrixQPEvaluator

  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotationMatrixQPEvaluator& e);

  [[nodiscard]] bool check(AbstractErrorHandler& eh,
                           const RotationMatrixQPEvaluator& e);

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
     */
    RotatedThermodynamicForcesMatrixQPEvaluator(
        const Material& m,
        const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler& ctx) const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // access operator
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

  //! \bref return the quadrature space
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e);
  //! \brief perform consistency checks
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] bool check(AbstractErrorHandler& eh,
                           const RotatedThermodynamicForcesMatrixQPEvaluator<
                               ThermodynamicForcesSize>& e);
  //! \brief return the number of components
  template <size_type ThermodynamicForcesSize>
  mgis::size_type getNumberOfComponents(
      const RotatedThermodynamicForcesMatrixQPEvaluator<
          ThermodynamicForcesSize>& e) noexcept;

  /*!
   * \brief an evaluator returning the gradients rotated in the material frame
   */
  template <size_type GradientsSize = dynamic_extent>
  struct RotatedGradientsMatrixQPEvaluator {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotatedGradientsMatrixQPEvaluator(
        const Material& m,
        const Material::StateSelection s = Material::END_OF_TIME_STEP);
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler& ctx) const;
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // \brief access operator
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief thermodynamic forces
    std::span<real> gradients;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    mutable Buffer<GradientsSize> buffer;
  };  // end of RotatedGradientsMatrixQPEvaluator

  //! \bref return the quadrature space
  template <size_type GradientsSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e);
  //! \brief perform consistency checks
  template <size_type GradientsSize>
  [[nodiscard]] bool check(
      AbstractErrorHandler& eh,
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e);
  //! \brief return the number of components
  template <size_type GradientsSize>
  mgis::size_type getNumberOfComponents(
      const RotatedGradientsMatrixQPEvaluator<GradientsSize>& e) noexcept;

  /*!
   * \brief check if the given evaluators have the same partial quadrature space
   * \param[in, out] ctx: execution context
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
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
