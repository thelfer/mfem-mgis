/*!
 * \file   src/MechanicalPostProcessings.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   22/05/2025
 */

#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/Mechanics.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

#include "MFEMMGIS/MechanicalPostProcessings.hxx"

namespace mfem_mgis {

#ifdef MGIS_FUNCTION_SUPPORT

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStressForSmallStrainBehaviours_impl(
      Context& ctx, const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction seq(m.getPartialQuadratureSpacePointer());
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok = sig | as_stensor<N> | vmis | seq;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: computation of the von Mises "
          "stress failed");
    }
    return seq;
  }  // end of computeVonMisesEquivalentStressForSmallStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStressForSmallStrainBehaviours(
      Context & ctx, const Material& m, const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported modelling hypothesis");
  }  // end of computeVonMisesEquivalentStressForSmallStrainBehaviours

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl(
      Context& ctx, const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction seq(m.getPartialQuadratureSpacePointer());
    const auto F = getThermodynamicForce(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "Stress", s);
    const auto ok =
        pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) | vmis | seq;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: computation of the von Mises "
          "stress failed");
    }
    return seq;
  }  // end of computeVonMisesEquivalentStressForFiniteStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStressForFiniteStrainBehaviours(
      Context & ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported modelling hypothesis");
  }  // end of computeVonMisesEquivalentStressForFiniteStrainBehaviours

  std::optional<PartialQuadratureFunction> computeVonMisesEquivalentStress(
      Context & ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours(ctx, m, s);
    } else if (m.b.btype == mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours(ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported behaviour type");
  }  // end of computeVonMisesEquivalentStress

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeEigenStressesForSmallStrainBehaviours_impl(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction vp(m.getPartialQuadratureSpacePointer(), 3);
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok = sig | as_stensor<N> | eigen_values<> | vp;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: computation of the von Mises stress failed");
    }
    return vp;
  }  // end of computeEigenStressesForSmallStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeEigenStressesForSmallStrainBehaviours(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeEigenStressesForSmallStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeEigenStressesForSmallStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported modelling hypothesis");
  }  // end of computeEigenStressesForSmallStrainBehaviours

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeEigenStressesForFiniteStrainBehaviours_impl(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction vp(m.getPartialQuadratureSpacePointer(), 3);
    const auto F = getThermodynamicForce(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "Stress", s);
    const auto ok = pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) |
                    eigen_values<> | vp;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: computation of the von Mises stress failed");
    }
    return vp;
  }  // end of computeEigenStressesForFiniteStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeEigenStressesForFiniteStrainBehaviours(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeEigenStressesForFiniteStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeEigenStressesForFiniteStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported modelling hypothesis");
  }  // end of computeEigenStressesForFiniteStrainBehaviours

  std::optional<PartialQuadratureFunction> computeEigenStresses(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeEigenStressesForSmallStrainBehaviours(ctx, m, s);
    } else if (m.b.btype == mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeEigenStressesForFiniteStrainBehaviours(ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported behaviour type");
  }  // end of computeEigenStresses

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeFirstEigenStressForSmallStrainBehaviours_impl(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction s1(m.getPartialQuadratureSpacePointer());
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok =
        sig | as_stensor<N> | eigen_values<> | maximum_component | s1;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeFirstEigenStress: computation of the von Mises stress "
          "failed");
    }
    return s1;
  }  // end of computeFirstEigenStressForSmallStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeFirstEigenStressForSmallStrainBehaviours(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeFirstEigenStressForSmallStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeFirstEigenStressForSmallStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported modelling hypothesis");
  }  // end of computeFirstEigenStressForSmallStrainBehaviours

  template <unsigned short N>
  static std::optional<PartialQuadratureFunction>
  computeFirstEigenStressForFiniteStrainBehaviours_impl(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    using namespace mgis::function;
    PartialQuadratureFunction s1(m.getPartialQuadratureSpacePointer());
    const auto F = getThermodynamicForce(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "Stress", s);
    const auto ok = pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) |
                    eigen_values<> | maximum_component | s1;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeFirstEigenStress: computation of the von Mises stress "
          "failed");
    }
    return s1;
  }  // end of computeFirstEigenStressForFiniteStrainBehaviours

  static std::optional<PartialQuadratureFunction>
  computeFirstEigenStressForFiniteStrainBehaviours(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeFirstEigenStressForFiniteStrainBehaviours_impl<2>(
          ctx, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeFirstEigenStressForFiniteStrainBehaviours_impl<3>(
          ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported modelling hypothesis");
  }  // end of computeFirstEigenStressForFiniteStrainBehaviours

  std::optional<PartialQuadratureFunction> computeFirstEigenStress(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeFirstEigenStressForSmallStrainBehaviours(ctx, m, s);
    } else if (m.b.btype == mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeFirstEigenStressForFiniteStrainBehaviours(ctx, m, s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported behaviour type");
  }  // end of computeFirstEigenStress

  //   template <unsigned short N>
  //   static std::optional<PartialQuadratureFunction>
  //   computeCauchyStressInGlobalFrame_impl(Context& ctx,
  //                                         const Material& m,
  //                                         const Material::StateSelection s) {
  //     PartialQuadratureFunction sig(m.getPartialQuadratureSpacePointer(),
  //                                   tfel::math::StensorSize<N>::value);
  //     const auto F = getThermodynamicForce(m, "DeformationGradient", s);
  //     const auto pk1 = getThermodynamicForce(m, "Stress", s);
  //     const auto ok = pk1 | as_tensor<N> | from_pk1_to_cauchy(F |
  //     as_tensor<N>) |
  //                     rotate_backwards(R) | (sig | as_stensor<N>);
  //     if (!ok) {
  //       return ctx.registerErrorMessage(
  //           "computeFirstEigenStress: computation of the von Mises stress "
  //           "failed");
  //     }
  //     return sig;
  //   }
  //
  //   std::optional<PartialQuadratureFunction>
  //   computeCauchyStressInGlobalFrame(
  //       Context& ctx, const Material& m, const Material::StateSelection s) {
  //     using namespace mgis::function;
  //     if (m.b.btype !=
  //         mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
  //       return ctx.registerErrorMessage(
  //           "computeCauchyStressInGlobalFrame: not a finite strain
  //           behaviour");
  //     }
  //     if (m.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
  //       return ctx.registerErrorMessage(
  //           "computeCauchyStressInGlobalFrame: material is not orthotropic");
  //     }
  //
  //     if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
  //       return computeCauchyStressInGlobalFrame_impl<2>(ctx, m, s);
  //     } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
  //       return computeCauchyStressInGlobalFrame_impl<3>(ctx, m, s);
  //     }
  //     return ctx.registerErrorMessage(
  //         "computeCauchyStressInGlobalFrame: unsupported modelling
  //         hypothesis");
  //   }  // end of computeCauchyStressInGlobalFrame

#endif /* MGIS_FUNCTION_SUPPORT */

} // end of namespace mfem_mgis

