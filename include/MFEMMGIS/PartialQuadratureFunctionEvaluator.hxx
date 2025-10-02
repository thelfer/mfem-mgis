/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionEvaluator.hxx
 * \brief
b * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Buffer.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  //! \brief an evaluator returning the rotation matrix
  struct RotationMatrixPartialQuadratureFunctionEvalutor {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotationMatrixPartialQuadratureFunctionEvalutor(const Material&);
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler&) const;
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    // access operator
    auto operator()(const size_type) const;

   private:
    //! \brief underlying material
    const Material& material;
  };  // end of RotationMatrixPartialQuadratureFunctionEvalutor


  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotationMatrixPartialQuadratureFunctionEvalutor&);

  [[nodiscard]] bool check(
      AbstractErrorHandler&,
      const RotationMatrixPartialQuadratureFunctionEvalutor&);

  [[nodiscard]] constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixPartialQuadratureFunctionEvalutor&) noexcept;

  /*!
   * \brief an evaluator returning the rotated thermodynamic forces
   */
  template <size_type ThermodynamicForcesSize = dynamic_extent>
  struct RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor(
        const Material&,
        const Material::StateSelection = Material::END_OF_TIME_STEP);
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler&) const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // access operator
    auto operator()(const size_type) const;

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
      // RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor

  //! \bref return the quadrature space
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>&);
  //! \brief perform consistency checks
  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] bool check(
      AbstractErrorHandler&,
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>&);
  //! \brief return the number of components
  template <size_type ThermodynamicForcesSize>
  mgis::size_type getNumberOfComponents(
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>&) noexcept;

  /*!
   * \brief an evaluator returning the gradients rotated in the material frame
   */
  template <size_type GradientsSize = dynamic_extent>
  struct RotatedGradientsMatrixPartialQuadratureFunctionEvalutor {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotatedGradientsMatrixPartialQuadratureFunctionEvalutor(
        const Material&,
        const Material::StateSelection = Material::END_OF_TIME_STEP);
    //! \brief perform consistency checks
    [[nodiscard]] bool check(AbstractErrorHandler&) const;
    //! \brief return the underlying partial quadrature space
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // \brief access operator
    auto operator()(const size_type) const;

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief thermodynamic forces
    std::span<real> gradients;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    mutable Buffer<GradientsSize> buffer;
  };  // end of RotatedGradientsMatrixPartialQuadratureFunctionEvalutor

  //! \bref return the quadrature space
  template <size_type GradientsSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
          GradientsSize>&);
  //! \brief perform consistency checks
  template <size_type GradientsSize>
  [[nodiscard]] bool check(
      AbstractErrorHandler&,
      const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
          GradientsSize>&);
  //! \brief return the number of components
  template <size_type GradientsSize>
  mgis::size_type getNumberOfComponents(
      const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
          GradientsSize>&) noexcept;

  /*!
   * \brief check if the given evaluators have the same partial quadrature space
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   */
  template <PartialQuadratureFunctionEvaluatorConcept EvaluatorType1,
            PartialQuadratureFunctionEvaluatorConcept EvaluatorType2>
  void checkMatchingQuadratureSpaces(const EvaluatorType1&,
                                     const EvaluatorType2&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureFunctionEvaluator.ixx"

namespace mfem_mgis{

  static_assert(mgis::function::EvaluatorConcept<
                RotationMatrixPartialQuadratureFunctionEvalutor>);
  static_assert(!mgis::function::FunctionConcept<
                RotationMatrixPartialQuadratureFunctionEvalutor>);
  static_assert(
      mgis::function::EvaluatorConcept<
          RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<>>);
  static_assert(
      !mgis::function::FunctionConcept<
          RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<>>);
  static_assert(mgis::function::EvaluatorConcept<
                RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<>>);
  static_assert(!mgis::function::FunctionConcept<
                RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<>>);

}

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_HXX */
