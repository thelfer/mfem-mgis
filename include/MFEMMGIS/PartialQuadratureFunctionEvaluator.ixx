/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_IXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_IXX

#include <iterator>
#include <algorithm>

namespace mfem_mgis::algorithm {

  /*!
   * \brief copy N values
   * \tparam N: number of values to be copied
   * \tparam InputIterator: input iterator type
   * \tparam OutputIterator: output iterator type
   * \param[in] p: iterator to the beginning of the values to be copied
   * \param[in] pe: iterator past the end of the values to be copied
   * \param[in] po: iterator to the output values
   */
  template <size_type N, typename InputIterator, typename OutputIterator>
  void copy(const InputIterator p,
            const InputIterator pe,
            OutputIterator po) requires(N > 0) {
    if constexpr ((std::random_access_iterator<InputIterator>)&&  //
                  (std::random_access_iterator<OutputIterator>)) {
      if constexpr (N > 9) {
        std::copy(p, pe, po);
      } else if constexpr (N == 1) {
        *po = *p;
      } else if constexpr (N == 2) {
        po[0] = p[0];
        po[1] = p[1];
      } else if constexpr (N == 3) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
      } else if constexpr (N == 4) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
        po[3] = p[3];
      } else {
        copy<N - 1>(++p, pe, ++po);
      }
    } else {
      std::copy(p, pe, po);
    }
  }  // end of copy

}  // end of namespace mfem_mgis::algorithm

namespace mfem_mgis {

  inline auto RotationMatrixPartialQuadratureFunctionEvalutor::operator()(
      const size_type i) const {
    return this->material.getRotationMatrixAtIntegrationPoint(i);
  }

  inline const PartialQuadratureSpace& getSpace(
      const RotationMatrixPartialQuadratureFunctionEvalutor& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  inline bool check(AbstractErrorHandler& eh,
                    const RotationMatrixPartialQuadratureFunctionEvalutor& e) {
    return e.check(eh);
  }  // end of check

  constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixPartialQuadratureFunctionEvalutor&) noexcept {
    return 9u;
  }  // end of getNumberOfComponents

  template <size_type ThermodynamicForcesSize>
  RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
      ThermodynamicForcesSize>::
      RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor(
          const Material& m, const Material::StateSelection s)
      : material(m),
        thforces(getStateManager(m, s).thermodynamic_forces),
        stage(s) {
    if constexpr (ThermodynamicForcesSize != dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto thsize = sm.thermodynamic_forces_stride;
      this->buffer.resize(thsize);
    }
  }  // end of RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor

  template <size_type ThermodynamicForcesSize>
  bool RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
      ThermodynamicForcesSize>::check(AbstractErrorHandler& ctx) const {
    if (this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage("material is not orthotropic");
    }
    if constexpr (ThermodynamicForcesSize != dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto thsize = sm.thermodynamic_forces_stride;
      if (ThermodynamicForcesSize != thsize) {
        return ctx.registerErrorMessage(
            "inconsistent number of components of the thermodynamic forces");
      }
    }
    return true;
  }

  template <size_type ThermodynamicForcesSize>
  const PartialQuadratureSpace&
  RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
      ThermodynamicForcesSize>::getPartialQuadratureSpace() const {
    return this->material.getPartialQuadratureSpace();
  }

  template <size_type ThermodynamicForcesSize>
  size_type RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
      ThermodynamicForcesSize>::getNumberOfComponents() const noexcept {
    if constexpr (ThermodynamicForcesSize == dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      return sm.thermodynamic_forces_stride;
    } else {
      return ThermodynamicForcesSize;
    }
  }

  template <size_type ThermodynamicForcesSize>
  auto RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
      ThermodynamicForcesSize>::operator()(const size_type i) const {
    const auto* const mf = [this, i] {
      if constexpr (ThermodynamicForcesSize == dynamic_extent) {
        return this->thforces.data() + i * (this->buffer.size());
      } else {
        return this->thforces.data() + i * ThermodynamicForcesSize;
      }
    }();
    const auto r = this->material.getRotationMatrixAtIntegrationPoint(i);
    this->material.b.rotate_thermodynamic_forces_ptr(this->buffer.data(), mf,
                                                     r.data());
    return makeSpan(this->buffer);
  }

  template <size_type ThermodynamicForcesSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  template <size_type ThermodynamicForcesSize>
  bool check(
      AbstractErrorHandler& eh,
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>& e) {
    return e.check(eh);
  }  // end of check

  template <size_type ThermodynamicForcesSize>
  inline mgis::size_type getNumberOfComponents(
      const RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor<
          ThermodynamicForcesSize>& e) noexcept {
    return e.getNumberOfComponents();
  }

  template <size_type GradientsSize>
  RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<GradientsSize>::
      RotatedGradientsMatrixPartialQuadratureFunctionEvalutor(
          const Material& m, const Material::StateSelection s)
      : material(m), gradients(getStateManager(m, s).gradients), stage(s) {
    if constexpr (GradientsSize == dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto gsize = sm.gradients_stride;
      this->buffer.resize(gsize);
    }
  }

  template <size_type GradientsSize>
  bool
  RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<GradientsSize>::check(
      AbstractErrorHandler& ctx) const {
    if (this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage("material is not orthotropic");
    }
    if constexpr (GradientsSize != dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      const auto gsize = sm.gradients_stride;
      if (GradientsSize != gsize) {
        return ctx.registerErrorMessage(
            "inconsistent number of components of the thermodynamic forces");
      }
    }
    return true;
  }  // end of check

  template <size_type GradientsSize>
  const PartialQuadratureSpace&
  RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
      GradientsSize>::getPartialQuadratureSpace() const {
    return this->material.getPartialQuadratureSpace();
  }

  template <size_type GradientsSize>
  size_type RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
      GradientsSize>::getNumberOfComponents() const noexcept {
    if constexpr (GradientsSize == dynamic_extent) {
      const auto& sm = getStateManager(this->material, this->stage);
      return sm.gradients_stride;
    } else {
      return GradientsSize;
    }
  }

  template <size_type GradientsSize>
  auto RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
      GradientsSize>::operator()(const size_type i) const {
    const auto* const mg = [this, i] {
      if constexpr (GradientsSize == dynamic_extent) {
        return this->gradients.data() + i * (this->buffer.size());
      } else {
        return this->gradients.data() + i * GradientsSize;
      }
    }();
    const auto r = this->material.getRotationMatrixAtIntegrationPoint(i);
    // transposed rotation matrix
    const auto tr = std::array<real, 9>{r[0], r[3], r[6],  //
                                        r[1], r[4], r[7],  //
                                        r[2], r[5], r[8]};
    this->material.b.rotate_gradients_ptr(this->buffer.data(), mg, tr.data());
    return makeSpan(this->buffer);
  }

  template <size_type GradientsSize>
  [[nodiscard]] const PartialQuadratureSpace& getSpace(
      const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
          GradientsSize>& e) {
    return e.getPartialQuadratureSpace();
  }  // end of getSpace

  template <size_type GradientsSize>
  bool check(AbstractErrorHandler& eh,
             const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
                 GradientsSize>& e) {
    return e.check(eh);
  }  // end of check

  template <size_type GradientsSize>
  inline void allocateWorkspace(
      RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<GradientsSize>&
          e) {
    e.allocateWorkspace();
  }

  template <size_type GradientsSize>
  inline mgis::size_type getNumberOfComponents(
      const RotatedGradientsMatrixPartialQuadratureFunctionEvalutor<
          GradientsSize>& e) noexcept {
    return e.getNumberOfComponents();
  }

  template <PartialQuadratureFunctionEvaluatorConcept EvaluatorType1,
            PartialQuadratureFunctionEvaluatorConcept EvaluatorType2>
  void checkMatchingQuadratureSpaces(const EvaluatorType1& e1,
                                     const EvaluatorType2& e2) {
    const auto& qspace1 = getSpace(e1);
    const auto& qspace2 = getSpace(e2);
    raise_if(&qspace1 != &qspace2, "unmatched quadrature spaces");
  }  // end of checkMatchingQuadratureSpaces

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATOR_IXX */
