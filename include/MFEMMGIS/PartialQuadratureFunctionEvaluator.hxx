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

  template <typename EvaluatorType>
  concept PartialQuadratureFunctionEvalutorConcept = requires(EvaluatorType& e) {
    e.allocateWorkspace();
  }
  &&requires(const EvaluatorType& e) {
    e.check();
    {e.getPartialQuadratureSpace()}->std::same_as<const PartialQuadratureSpace&>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
  }
  &&((requires(const EvaluatorType& e, size_type i) { e(i); }) ||
     (requires(EvaluatorType & e, size_type i) { e(i); }));

  /*!
   * \brief an evaluator returning the values of an immutable partial quadrature
   * function view as a fixed size span or a scalar
   * \tparam size of the returned value
   */
  template <size_type N>
  struct FixedSizedPartialQuadratureFunctionEvalutor {
    /*!
     * \brief constructor
     */
    FixedSizedPartialQuadratureFunctionEvalutor(
        const ImmutablePartialQuadratureFunctionView&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    // access operator
    auto operator()(const size_type i) const;

   private:
    //! \brief underlying partial quadrature space
    const ImmutablePartialQuadratureFunctionView& function;
  };  // end of FixedSizedPartialQuadratureFunctionEvalutor

  //! \brief an evaluator returning the rotation matrix
  struct RotationMatrixPartialQuadratureFunctionEvalutor {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    RotationMatrixPartialQuadratureFunctionEvalutor(const Material&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief return the underlying partial quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    // access operator
    auto operator()(const size_type) const;

   private:
    //! \brief underlying material
    const Material& material;
  };  // end of RotationMatrixPartialQuadratureFunctionEvalutor

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
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // access operator
    auto operator()(const size_type);

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief thermodynamic forces
    std::span<real> thforces;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    Buffer<ThermodynamicForcesSize> buffer;
  };  // end of
      // RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor

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
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const PartialQuadratureSpace& getPartialQuadratureSpace() const;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;
    // access operator
    auto operator()(const size_type);

   private:
    //! \brief underlying material
    const Material& material;
    //! \brief thermodynamic forces
    std::span<real> gradients;
    //! \brief time step stage
    const Material::StateSelection stage;
    //! \brief buffer
    Buffer<GradientsSize> buffer;
  };  // end of RotatedGradientsMatrixPartialQuadratureFunctionEvalutor

  /*!
   * \brief check if the given evaluators have the same partial quadrature space
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   */
  template <PartialQuadratureFunctionEvalutorConcept EvaluatorType1,
            PartialQuadratureFunctionEvalutorConcept EvaluatorType2>
  void checkMatchingQuadratureSpaces(const EvaluatorType1&,
                                     const EvaluatorType2&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureFunctionEvaluator.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_HXX */
