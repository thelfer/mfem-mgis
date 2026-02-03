/*!
 * \file   src/MechanicalPostProcessings.cxx
 * \brief
 * \author Thomas Helfer
 * \date   22/05/2025
 */

#include "MGIS/Function/Mechanics.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"

namespace mfem_mgis {

  template <unsigned short N>
  static bool computeVonMisesEquivalentStressForSmallStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& seq,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok = sig | as_stensor<N> | vmis | seq;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: computation of the von Mises "
          "stress failed");
    }
    return true;
  }  // end of computeVonMisesEquivalentStressForSmallStrainBehaviours

  static bool computeVonMisesEquivalentStressForSmallStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& seq,
      const Material& m,
      const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours_impl<2>(
          ctx, seq, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours_impl<3>(
          ctx, seq, m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported modelling hypothesis");
  }  // end of computeVonMisesEquivalentStressForSmallStrainBehaviours

  template <unsigned short N>
  static bool computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& seq,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto F = getGradient(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    const auto ok =
        pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) | vmis | seq;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: computation of the von Mises "
          "stress failed");
    }
    return true;
  }  // end of computeVonMisesEquivalentStressForFiniteStrainBehaviours

  static bool computeVonMisesEquivalentStressForFiniteStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& seq,
      const Material& m,
      const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl<2>(
          ctx, seq, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours_impl<3>(
          ctx, seq, m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported modelling hypothesis");
  }  // end of computeVonMisesEquivalentStressForFiniteStrainBehaviours

  std::optional<PartialQuadratureFunction> computeVonMisesEquivalentStress(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    PartialQuadratureFunction seq(m.getPartialQuadratureSpacePointer());
    if (!computeVonMisesEquivalentStress(ctx, seq, m, s)) {
      return {};
    }
    return seq;
  }  // end of computeVonMisesEquivalentStress

  bool computeVonMisesEquivalentStress(Context& ctx,
                                       PartialQuadratureFunction& seq,
                                       const Material& m,
                                       const Material::StateSelection s) {
    if (seq.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: quadrature function is not defined "
          "on the given material");
    }
    if (seq.getNumberOfComponents() != 1) {
      return ctx.registerErrorMessage(
          "computeVonMisesEquivalentStress: quadrature function is not scalar");
    }
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeVonMisesEquivalentStressForSmallStrainBehaviours(ctx, seq,
                                                                     m, s);
    } else if (m.b.btype ==
               mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeVonMisesEquivalentStressForFiniteStrainBehaviours(ctx, seq,
                                                                      m, s);
    }
    return ctx.registerErrorMessage(
        "computeVonMisesEquivalentStress: unsupported behaviour type");
  }  // end of computeVonMisesEquivalentStress

  template <unsigned short N>
  static bool computeEigenStressesForSmallStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& svp,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto sig = getThermodynamicForce(m, "Stress", s);
    auto svp_view = svp | as_tvector<3>;
    const auto ok = sig | as_stensor<N> | eigen_values<> | svp_view;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: computation of the von Mises stress failed");
    }
    return true;
  }  // end of computeEigenStressesForSmallStrainBehaviours

  static bool computeEigenStressesForSmallStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& svp,
      const Material& m,
      const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeEigenStressesForSmallStrainBehaviours_impl<2>(ctx, svp, m,
                                                                  s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeEigenStressesForSmallStrainBehaviours_impl<3>(ctx, svp, m,
                                                                  s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported modelling hypothesis");
  }  // end of computeEigenStressesForSmallStrainBehaviours

  template <unsigned short N>
  static bool computeEigenStressesForFiniteStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& svp,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto F = getGradient(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    auto svp_view = svp | as_tvector<3>;
    const auto ok = pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) |
                    eigen_values<> | svp_view;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: computation of the von Mises stress failed");
    }
    return true;
  }  // end of computeEigenStressesForFiniteStrainBehaviours

  static bool computeEigenStressesForFiniteStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& svp,
      const Material& m,
      const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeEigenStressesForFiniteStrainBehaviours_impl<2>(ctx, svp, m,
                                                                   s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeEigenStressesForFiniteStrainBehaviours_impl<3>(ctx, svp, m,
                                                                   s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported modelling hypothesis");
  }  // end of computeEigenStressesForFiniteStrainBehaviours

  bool computeEigenStresses(Context& ctx,
                            PartialQuadratureFunction& svp,
                            const Material& m,
                            const Material::StateSelection s) {
    if (svp.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not defined "
          "on the given material");
    }
    if (svp.getNumberOfComponents() != 3) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not scalar");
    }
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeEigenStressesForSmallStrainBehaviours(ctx, svp, m, s);
    } else if (m.b.btype ==
               mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeEigenStressesForFiniteStrainBehaviours(ctx, svp, m, s);
    }
    return ctx.registerErrorMessage(
        "computeEigenStresses: unsupported behaviour type");
  }  // end of computeEigenStresses

  std::optional<PartialQuadratureFunction> computeEigenStresses(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    PartialQuadratureFunction svp(m.getPartialQuadratureSpacePointer(), 3);
    if (!computeEigenStresses(ctx, svp, m, s)) {
      return {};
    }
    return svp;
  }  // end of computeEigenStresses

  template <unsigned short N>
  static bool computeFirstEigenStressForSmallStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& s1,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok =
        sig | as_stensor<N> | eigen_values<> | maximum_component | s1;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeFirstEigenStress: computation of the von Mises stress "
          "failed");
    }
    return true;
  }  // end of computeFirstEigenStressForSmallStrainBehaviours

  static bool computeFirstEigenStressForSmallStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& s1,
      const Material& m,
      const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeFirstEigenStressForSmallStrainBehaviours_impl<2>(ctx, s1, m,
                                                                     s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeFirstEigenStressForSmallStrainBehaviours_impl<3>(ctx, s1, m,
                                                                     s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported modelling hypothesis");
  }  // end of computeFirstEigenStressForSmallStrainBehaviours

  template <unsigned short N>
  static bool computeFirstEigenStressForFiniteStrainBehaviours_impl(
      Context& ctx,
      PartialQuadratureFunction& s1,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto F = getGradient(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    const auto ok = pk1 | as_tensor<N> | from_pk1_to_cauchy(F | as_tensor<N>) |
                    eigen_values<> | maximum_component | s1;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeFirstEigenStress: computation of the von Mises stress "
          "failed");
    }
    return true;
  }  // end of computeFirstEigenStressForFiniteStrainBehaviours

  static bool computeFirstEigenStressForFiniteStrainBehaviours(
      Context& ctx,
      PartialQuadratureFunction& s1,
      const Material& m,
      const Material::StateSelection s) {
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      return computeFirstEigenStressForFiniteStrainBehaviours_impl<2>(ctx, s1,
                                                                      m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeFirstEigenStressForFiniteStrainBehaviours_impl<3>(ctx, s1,
                                                                      m, s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported modelling hypothesis");
  }  // end of computeFirstEigenStressForFiniteStrainBehaviours

  bool computeFirstEigenStress(Context& ctx,
                               PartialQuadratureFunction& s1,
                               const Material& m,
                               const Material::StateSelection s) {
    if (s1.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not defined "
          "on the given material");
    }
    if (s1.getNumberOfComponents() != 1) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not scalar");
    }
    if (m.b.btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      return computeFirstEigenStressForSmallStrainBehaviours(ctx, s1, m, s);
    } else if (m.b.btype ==
               mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return computeFirstEigenStressForFiniteStrainBehaviours(ctx, s1, m, s);
    }
    return ctx.registerErrorMessage(
        "computeFirstEigenStress: unsupported behaviour type");
  }  // end of computeFirstEigenStress

  std::optional<PartialQuadratureFunction> computeFirstEigenStress(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    PartialQuadratureFunction s1(m.getPartialQuadratureSpacePointer());
    if (!computeFirstEigenStress(ctx, s1, m, s)) {
      return {};
    }
    return s1;
  }

  template <unsigned short N>
  static bool computeSmallStrainStressInGlobalFrame_impl(
      Context& ctx,
      PartialQuadratureFunction& rsig,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    if (rsig.getNumberOfComponents() !=
        tfel::math::StensorDimeToSize<N>::value) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: invalid quadrature function size");
    }
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto R = RotationMatrixEvaluator{m};
    auto sview = rsig | as_stensor<N>;
    const auto ok =
        sig | as_stensor<N> | rotate_backwards(R | as_tmatrix<3, 3>) | sview;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeStressInGlobalFrame: computation of the  "
          "stress in the global frame failed");
    }
    return true;
  }  // end of computeSmallStrainStressInGlobalFrame_impl

  template <unsigned short N>
  static bool computeFirstPiolaKirchhoffStressInGlobalFrame_impl(
      Context& ctx,
      PartialQuadratureFunction& rpk1,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    if (rpk1.getNumberOfComponents() !=
        tfel::math::TensorDimeToSize<N>::value) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: invalid quadrature function size");
    }
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    const auto R = RotationMatrixEvaluator{m};
    auto rpk1_view = rpk1 | as_tensor<N>;
    const auto ok =
        pk1 | as_tensor<N> | rotate_backwards(R | as_tmatrix<3, 3>) | rpk1_view;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeStressInGlobalFrame: computation of the  "
          "stress in the global frame failed");
    }
    return true;
  }  // end of computeFirstPiolaKirchhoffStressInGlobalFrame_impl

  bool computeStressInGlobalFrame(Context& ctx,
                                  PartialQuadratureFunction& rstress,
                                  const Material& m,
                                  const Material::StateSelection s) {
    if (rstress.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not defined "
          "on the given material");
    }
    if (m.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage(
          "computeStressInGlobalFrame: material is not orthotropic");
    }
    if (m.b.btype != mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      if ((m.b.hypothesis == Hypothesis::PLANESTRAIN) ||
          (m.b.hypothesis == Hypothesis::PLANESTRESS)) {
        return computeSmallStrainStressInGlobalFrame_impl<2>(ctx, rstress, m,
                                                             s);
      } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
        return computeSmallStrainStressInGlobalFrame_impl<3>(ctx, rstress, m,
                                                             s);
      } else {
        return ctx.registerErrorMessage(
            "computeStressInGlobalFrame: unsupported modelling hypothesis");
      }
    } else if (m.b.btype !=
               mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      if ((m.b.hypothesis == Hypothesis::PLANESTRAIN) ||
          (m.b.hypothesis == Hypothesis::PLANESTRESS)) {
        return computeFirstPiolaKirchhoffStressInGlobalFrame_impl<2>(
            ctx, rstress, m, s);
      } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
        return computeFirstPiolaKirchhoffStressInGlobalFrame_impl<3>(
            ctx, rstress, m, s);
      } else {
        return ctx.registerErrorMessage(
            "computeStressInGlobalFrame: unsupported modelling hypothesis");
      }
    }
    return ctx.registerErrorMessage(
        "computeStressInGlobalFrame: not a mechanical behaviour");
  }  // end of computeStressInGlobalFrame

  std::optional<PartialQuadratureFunction> computeStressInGlobalFrame(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    auto stress_size = size_type{};
    if (m.b.btype != mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      if ((m.b.hypothesis == Hypothesis::PLANESTRAIN) ||
          (m.b.hypothesis == Hypothesis::PLANESTRESS)) {
        stress_size = 4;
      } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
        stress_size = 6;
      } else {
        return ctx.registerErrorMessage(
            "computeStressInGlobalFrame: unsupported modelling hypothesis");
      }
    } else if (m.b.btype !=
               mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      if ((m.b.hypothesis == Hypothesis::PLANESTRAIN) ||
          (m.b.hypothesis == Hypothesis::PLANESTRESS)) {
        stress_size = 5;
      } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
        stress_size = 9;
      } else {
        return ctx.registerErrorMessage(
            "computeStressInGlobalFrame: unsupported modelling hypothesis");
      }
    } else {
      return ctx.registerErrorMessage(
          "computeStressInGlobalFrame: not a mechanical behaviour");
    }
    PartialQuadratureFunction rpk1(m.getPartialQuadratureSpacePointer(),
                                   stress_size);
    if (!computeStressInGlobalFrame(ctx, rpk1, m, s)) {
      return {};
    }
    return rpk1;
  }

  template <unsigned short N>
  static bool computeCauchyStressInGlobalFrame_impl(
      Context& ctx,
      PartialQuadratureFunction& sig,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto F = getGradient(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    const auto R = RotationMatrixEvaluator{m};
    auto sview = sig | as_stensor<N>;
    const auto ok = pk1 | as_tensor<N> |  //
                    from_pk1_to_cauchy(F | as_tensor<N>) |
                    rotate_backwards(R | as_tmatrix<3, 3>) | sview;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeCauchyStressInGlobalFrame: computation of the Cauchy "
          "stress in the global frame failed");
    }
    return true;
  }  // end of computeCauchyStressInGlobalFrame_impl

  bool computeCauchyStressInGlobalFrame(Context& ctx,
                                        PartialQuadratureFunction& sig,
                                        const Material& m,
                                        const Material::StateSelection s) {
    if (sig.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not defined "
          "on the given material");
    }
    if (m.b.btype !=
        mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage(
          "computeCauchyStressInGlobalFrame: not a finite strain behaviour");
    }
    if (m.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      return ctx.registerErrorMessage(
          "computeCauchyStressInGlobalFrame: material is not orthotropic");
    }
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      if (sig.getNumberOfComponents() != 4) {
        return ctx.registerErrorMessage(
            "computeCauchyStressInGlobalFrame: invalid quadrature function "
            "size");
      }
      return computeCauchyStressInGlobalFrame_impl<2>(ctx, sig, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      if (sig.getNumberOfComponents() != 6) {
        return ctx.registerErrorMessage(
            "computeCauchyStressInGlobalFrame: invalid quadrature function "
            "size");
      }
      return computeCauchyStressInGlobalFrame_impl<3>(ctx, sig, m, s);
    } else {
      return ctx.registerErrorMessage(
          "computeCauchyStressInGlobalFrame: unsupported modelling "
          "hypothesis");
    }
    return true;
  }  // end of computeCauchyStressInGlobalFrame

  std::optional<PartialQuadratureFunction> computeCauchyStressInGlobalFrame(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    auto stress_size = size_type{};
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      stress_size = 4;
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      stress_size = 6;
    } else {
      return ctx.registerErrorMessage(
          "computeCauchyStressInGlobalFrame: unsupported modelling "
          "hypothesis");
    }
    PartialQuadratureFunction sig(m.getPartialQuadratureSpacePointer(),
                                  stress_size);
    if (!computeCauchyStressInGlobalFrame(ctx, sig, m, s)) {
      return {};
    }
    return sig;
  }  // end of computeCauchyStressInGlobalFrame

  template <unsigned short N>
  static bool computeCauchyStress_impl(Context& ctx,
                                       PartialQuadratureFunction& sig,
                                       const Material& m,
                                       const Material::StateSelection s) {
    using namespace mgis::function;
    const auto F = getGradient(m, "DeformationGradient", s);
    const auto pk1 = getThermodynamicForce(m, "FirstPiolaKirchhoffStress", s);
    auto sview = sig | as_stensor<N>;
    const auto ok = pk1 | as_tensor<N> |  //
                    from_pk1_to_cauchy(F | as_tensor<N>) | sview;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeCauchyStress: computation of the Cauchy "
          "stress in the global frame failed");
    }
    return true;
  }  // end of computeCauchyStress_impl

  bool computeCauchyStress(Context& ctx,
                           PartialQuadratureFunction& sig,
                           const Material& m,
                           const Material::StateSelection s) {
    if (sig.getPartialQuadratureSpacePointer() !=
        m.getPartialQuadratureSpacePointer()) {
      return ctx.registerErrorMessage(
          "computeEigenStresses: quadrature function is not defined "
          "on the given material");
    }
    if (m.b.btype !=
        mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage(
          "computeCauchyStress: not a finite strain behaviour");
    }
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      if (sig.getNumberOfComponents() != 4) {
        return ctx.registerErrorMessage(
            "computeCauchyStress: invalid quadrature function "
            "size");
      }
      return computeCauchyStress_impl<2>(ctx, sig, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      if (sig.getNumberOfComponents() != 6) {
        return ctx.registerErrorMessage(
            "computeCauchyStress: invalid quadrature function "
            "size");
      }
      return computeCauchyStress_impl<3>(ctx, sig, m, s);
    } else {
      return ctx.registerErrorMessage(
          "computeCauchyStress: unsupported modelling "
          "hypothesis");
    }
    return true;
  }  // end of computeCauchyStress

  std::optional<PartialQuadratureFunction> computeCauchyStress(
      Context& ctx, const Material& m, const Material::StateSelection s) {
    auto stress_size = size_type{};
    if (m.b.hypothesis == Hypothesis::PLANESTRAIN) {
      stress_size = 4;
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      stress_size = 6;
    } else {
      return ctx.registerErrorMessage(
          "computeCauchyStress: unsupported modelling "
          "hypothesis");
    }
    PartialQuadratureFunction sig(m.getPartialQuadratureSpacePointer(),
                                  stress_size);
    if (!computeCauchyStress(ctx, sig, m, s)) {
      return {};
    }
    return sig;
  }  // end of computeCauchyStressInGlobalFrame

}  // end of namespace mfem_mgis
