/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX
#define LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for behaviour integrators
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegratorBase : BehaviourIntegrator,
                                                    Material {
    void setTimeIncrement(const real) override;
    void setup(const real, const real) override;
    void revert() override;
    void update() override;
    Material& getMaterial() override;
    const Material& getMaterial() const override;
    void setMacroscopicGradients(mgis::span<const real>) override;
    //! \brief destructor
    ~BehaviourIntegratorBase() override;

   protected:
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    BehaviourIntegratorBase(std::shared_ptr<const PartialQuadratureSpace>,
                            std::unique_ptr<const Behaviour>);
    /*!
     * \brief check that the integrator hypothesis is the same than the
     * behaviour hypothesis.
     * \param[in] h: integrator' hypothesis
     *
     * \throw an `std::runtime_error` if the hypotheses don't match.
     */
    void checkHypotheses(const Hypothesis) const;
    /*!
     * \brief throw an exception stating that the behaviour type is not the
     * expected one.
     * \param[in] mn: calling method name
     * \param[in] m: error message
     */
    [[noreturn]] void throwInvalidBehaviourType(const char* const,
                                                const char* const) const;
    /*!
     * \brief throw an exception stating that the behaviour kinematic is not the
     * expected one.
     * \param[in] mn: calling method name
     * \param[in] m: error message
     */
    [[noreturn]] void throwInvalidBehaviourKinematic(const char* const,
                                                     const char* const) const;
    /*!
     * \brief integrate the mechanical behaviour over the time step
     * If successful, the value of the stress, consistent tangent
     * operator and internal state variables are updated.
     *
     * \param[in] ip: local integration point index
     * \note this method shall be called after having set the gradients.
     */
    virtual void integrate(const size_type);
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
      std::vector<std::tuple<size_type, real*>> mps_evaluators;
      /*!
       * \brief evaluators for the external state variables at the beginning of
       * the time step
       */
      std::vector<std::tuple<size_type, real*>> esvs0_evaluators;
      /*!
       * \brief evaluators for the external state variables at the end of
       * the time step
       */
      std::vector<std::tuple<size_type, real*>> esvs1_evaluators;
    } wks;
    //! \brief time increment for the given time step
    real time_increment;
  };  // end of struct BehaviourIntegratorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX */
