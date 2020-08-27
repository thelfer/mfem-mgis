/*!
 * \file   MFEMMGIS/BehaviourIntegrator.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis{

  /*!
   * \brief 
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegrator {
    virtual void computeInnerForces(const mfem::FiniteElement &,
                                    mfem::ElementTransformation &,
                                    const mfem::Vector &,
                                    mfem::Vector &) = 0;

    virtual void computeStiffnessMatrix(const mfem::FiniteElement &,
                                     mfem::ElementTransformation &,
                                     const mfem::Vector &,
                                     mfem::DenseMatrix &) = 0;

    virtual ~BehaviourIntegrator();

   protected:
    /*!
     * \brief return a suitable integration rule for the given finite
     * element
     * and the finite element transformation.
     */
    virtual const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &,
        const mfem::ElementTransformation &) const = 0;
    /*!
     * \brief load a behaivour
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    static std::unique_ptr<mgis::behaviour::Behaviour> load(
        const std::string &, const std::string &);
  }; // end of struct BehaviourIntegrator

} // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX */
