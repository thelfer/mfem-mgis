/*!
 * \file MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   17/03/2026
 */

#ifndef LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX
#define LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"
#include "IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  template <Hypothesis H>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch;

  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::PLANESTRAIN> {
    using type =
        IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::PLANESTRESS> {
    using type =
        IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  template <>
  struct Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
      Hypothesis::TRIDIMENSIONAL> {
    using type =
        IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator;
  };

  template <Hypothesis H>
  using Faltus2026RegularizedIsotropicBehaviourIntegratorBase =
      typename Faltus2026RegularizedIsotropicBehaviourIntegratorBaseDispatch<
          H>::type;

  /*!
   * \brief a behaviour integrator which enhances a standard finite strain
   * behaviour with the regularization proposed by Faltus et al.
   *
   * This regularization adds a term penalizing the difference between the
   * deformation gradient \f$\underline{F}\f$ and its value
   * \f$\bar{\underline{F}}\f$ at the centroid of the element:
   *
   * \f[
   * W\left(\underline{F}, \bar{\underline{F}}\right) =
   * \alpha\,\left(\underline{F}-\bar{\underline{F}}\right)\,\colon\,
   * \left(\underline{F}-\bar{\underline{F}}\right)
   * \f]
   */
  template <Hypothesis H>
  struct Faltus2026RegularizedIsotropicBehaviourIntegrator final
      : Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H> {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     * \param[in] params: parameters defining the penalization coefficient
     */
    Faltus2026RegularizedIsotropicBehaviourIntegrator(
        const FiniteElementDiscretization& fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr,
        const Parameters& params);
    //
    void updateResidual(mfem::Vector& Fe,
                        const mfem::FiniteElement& e,
                        mfem::ElementTransformation& tr,
                        const mfem::Vector& u) override;
    void updateJacobian(mfem::DenseMatrix& Je,
                        const mfem::FiniteElement& e,
                        mfem::ElementTransformation& tr,
                        const mfem::Vector& u) override;
    [[nodiscard]] bool requiresCurrentSolutionForResidualAssembly()
        const noexcept override;
    [[nodiscard]] bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept override;

   private:
    //! \brief penalization coefficient
    const real alpha;
#ifndef MFEM_THREAD_SAFE
    /*!
     * \brief matrix used to store the derivatives of the shape functions at
     * the center of the element
     */
    mfem::DenseMatrix dshape0;
#endif
  };

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.ixx"

#endif /* LIB_MFEM_MGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_HXX */
