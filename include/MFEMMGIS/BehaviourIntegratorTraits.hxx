/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorTraits.hxx
 * \brief
 * \author Thomas Helfer
 * \date   06/04/2021
 */

#ifndef LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX
#define LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief a traits class aimed a specializing the
   * `StandardBehaviourIntegratorCRTPBase` class.
   */
  template <typename BehaviourIntegrator>
  struct BehaviourIntegratorTraits {
    //!
    static constexpr size_type unknownsSize = 0u;
    //! \brief
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    //! \brief
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = false;
    //! \brief
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORTRAITS_HXX */
