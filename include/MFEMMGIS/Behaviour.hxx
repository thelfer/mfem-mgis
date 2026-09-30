/*!
 * \file   include/MFEMMGIS/Behaviour.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOUR_HXX
#define LIB_MFEM_MGIS_BEHAVIOUR_HXX

#include <memory>
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MGISForward.hxx"

namespace mfem_mgis {

  //! \brief a simple alias
  using Behaviour = mgis::behaviour::Behaviour;

  /*!
   * \brief load a behaviour.
   *
   * Compared to the `mgis::behaviour::load` function, this function
   * handles specifically the case of finite strain behaviours.
   *
   * \param[in, out] ctx: execution context
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   * \return the loaded behaviour, a null pointer on failure
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::unique_ptr<Behaviour> load(
      Context& ctx,
      const std::string& l,
      const std::string& b,
      const Hypothesis h) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOUR_HXX */
