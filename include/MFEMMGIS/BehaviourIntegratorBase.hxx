/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX
#define LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX

#include <memory>
#include <string>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for behaviour integrators based on MFront
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegratorBase : AbstractBehaviourIntegrator,
                                                    Material {
    [[nodiscard]] const PartialQuadratureSpace& getPartialQuadratureSpace()
        const noexcept override;
    [[nodiscard]] real getTimeIncrement() const noexcept override;
    void setTimeIncrement(const real dt) override;
    [[nodiscard]] bool setup(Context& ctx,
                             const real t,
                             const real dt) noexcept override;
    [[nodiscard]] bool cleanup(Context& ctx) noexcept override;
    void revert() override;
    void update() override;
    [[nodiscard]] bool hasMaterial() const noexcept override;
    OptionalReference<Material> getMaterial(Context& ctx) noexcept override;
    OptionalReference<const Material> getMaterial(
        Context& ctx) const noexcept override;
    void setMacroscopicGradients(std::span<const real> g) override;
    [[nodiscard]] bool requiresCurrentSolutionForResidualAssembly()
        const noexcept override;
    [[nodiscard]] bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept override;
    [[nodiscard]] bool setMaterialProperty(
        Context& ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept override;
    [[nodiscard]] bool setExternalStateVariable(
        Context& ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept override;
    //! \brief destructor
    ~BehaviourIntegratorBase() override;

   protected:
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    BehaviourIntegratorBase(std::shared_ptr<const PartialQuadratureSpace> s,
                            std::unique_ptr<const Behaviour> b_ptr);
    /*!
     * \brief check that the behaviour is a standard finite strain behaviour
     * or a general behaviour whose only gradient is the deformation
     * gradient, whose only thermodynamic force is the first Piola-Kirchhoff
     * stress and whose only tangent operator block is the derivative of the
     * latter with respect to the former.
     */
    void checkIfAFiniteStrainBehaviourIsDeclared(
        attributes::Throwing throwing) const;
    /*!
     * \brief check if the behaviour has the expected symmetry.
     * \param[in] s: expected symmetry
     */
    void checkBehaviourSymmetry(attributes::Throwing throwing,
                                const Behaviour::Symmetry s) const;
    /*!
     * \brief check that the integrator hypothesis is the same than the
     * behaviour hypothesis.
     * \param[in] h: integrator' hypothesis
     */
    void checkHypothesis(attributes::Throwing throwing,
                         const Hypothesis h) const;
    /*!
     * \brief throw an exception stating that the behaviour type is not the
     * expected one.
     * \param[in] e: error message
     */
    [[noreturn]] void throwInvalidBehaviourType(attributes::Throwing throwing,
                                                const std::string& e) const;
    /*!
     * \brief throw an exception stating that the behaviour kinematic is not the
     * expected one.
     * \param[in] e: error message
     */
    [[noreturn]] void throwInvalidBehaviourKinematic(
        attributes::Throwing throwing, const std::string& e) const;
    /*!
     * \brief throw an exception stating that the behaviour symmetry is not the
     * expected one.
     * \param[in] e: error message
     */
    [[noreturn]] void throwInvalidBehaviourSymmetry(
        attributes::Throwing throwing, const std::string& e) const;
    /*!
     * \brief integrate the mechanical behaviour over the time step
     * If successful, the value of the stress, consistent tangent
     * operator and internal state variables are updated.
     * \return true if the integration is successful.
     * \param[in] ip: local integration point index
     * \param[in] it: integration type
     * \note this method shall be called after having set the gradients.
     */
    virtual bool performsLocalBehaviourIntegration(const size_type ip,
                                                   const IntegrationType it);
    /*!
     * \brief evaluators of the material properties at the beginning of the
     * time step.
     */
    std::map<std::string,
             std::shared_ptr<const AbstractQPEvaluator>,
             std::less<>>
        material_properties_evaluators_bts;
    /*!
     * \brief evaluators of the material properties at the end of the
     * time step.
     */
    std::map<std::string,
             std::shared_ptr<const AbstractQPEvaluator>,
             std::less<>>
        material_properties_evaluators_ets;
    /*!
     * \brief evaluators of the external state variables at the beginning of the
     * time step.
     */
    std::map<std::string,
             std::shared_ptr<const AbstractQPEvaluator>,
             std::less<>>
        external_state_variables_evaluators_bts;
    /*!
     * \brief evaluators of the external state variables at the end of the
     * time step.
     */
    std::map<std::string,
             std::shared_ptr<const AbstractQPEvaluator>,
             std::less<>>
        external_state_variables_evaluators_ets;
    //! \brief workspace
    struct {
      //! \brief array for material properties at the end of the time step
      std::vector<real> mps;
      //! \brief array for external state variables at the beginning of the time
      //! step
      std::vector<real> esvs0;
      //! \brief array for external state variables at the end of the time step
      std::vector<real> esvs1;
      /*!
       * \brief evaluators for the material properties at the end of the time
       * step
       */
      std::vector<std::tuple<size_type, size_type, const real*>> mps_evaluators;
      /*!
       * \brief evaluators for the external state variables at the beginning of
       * the time step
       */
      std::vector<std::tuple<size_type, size_type, const real*>>
          esvs0_evaluators;
      /*!
       * \brief evaluators for the external state variables at the end of
       * the time step
       */
      std::vector<std::tuple<size_type, size_type, const real*>>
          esvs1_evaluators;
      /*!
       * \brief partial quadrature functions resulting from the evaluations of
       * the evaluators associated with material properties at the beginning of
       * the time step
       */
      std::map<std::string, QPEvaluatorResult> pqfcts_mps_bts;
      /*!
       * \brief partial quadrature functions resulting from the evaluations of
       * the evaluators associated with material properties at the end of
       * the time step
       */
      std::map<std::string, QPEvaluatorResult> pqfcts_mps_ets;
      /*!
       * \brief partial quadrature functions resulting from the evaluations of
       * the evaluators associated with external state variables at the
       * beginning of the time step
       */
      std::map<std::string, QPEvaluatorResult> pqfcts_esvs_bts;
      /*!
       * \brief partial quadrature functions resulting from the evaluations of
       * the evaluators associated with external state variables at the end of
       * the time step
       */
      std::map<std::string, QPEvaluatorResult> pqfcts_esvs_ets;
    } wks;
    //! \brief time increment for the given time step
    real time_increment;
  };  // end of struct BehaviourIntegratorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX */
