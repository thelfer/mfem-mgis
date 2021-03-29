/*!
 * \file   MGISForward.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_MGISFORWARD_HXX
#define LIB_MFEM_MGIS_MGISFORWARD_HXX

#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis {

  namespace behaviour {

    struct Behaviour;

  }  // end of namespace behaviour

}  // end of namespace mgis

namespace mfem_mgis {

  //! \brief a simple alias
  using Hypothesis = mgis::behaviour::Hypothesis;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MGISFORWARD_HXX */
