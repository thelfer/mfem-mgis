/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorTraits.hxx
 * \brief  This file declares the `BehaviourIntegratorTraits` class
 * \author Thomas Helfer
 * \date   06/04/2021
 */

#ifndef LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX
#define LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief a traits class aimed at specializing the
   * `StandardBehaviourIntegratorCRTPBase` and
   * `FBarBehaviourIntegratorCRTPBase` classes.
   */
  template <typename BehaviourIntegrator>
  struct BehaviourIntegratorTraits {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 0u;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = false;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX */
