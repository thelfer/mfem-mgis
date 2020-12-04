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
#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for behaviour integrators
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegratorBase : BehaviourIntegrator,
                                                    Material {
    void setTimeIncrement(const real) override;
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
     * \brief integrate the mechanical behaviour over the time step
     * If successful, the value of the stress, consistent tangent 
     * operator and internal state variables are updated.
     *
     * \param[in] ip: local integration point index
     * \note this method shall be called after having set the gradients.
     */
    virtual void integrate(const size_type);
    //! \brief time increment for the given time step
    real time_increment;
  };  // end of struct BehaviourIntegratorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATORBASE_HXX */
