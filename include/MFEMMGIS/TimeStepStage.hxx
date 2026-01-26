/*!
 * \file   MFEMMGIS/TimeStepStage.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/01/2026
 */

#ifndef LIB_MFEMMGIS_TIMESTEPSTAGE_HXX
#define LIB_MFEMMGIS_TIMESTEPSTAGE_HXX

namespace mfem_mgis {

  /*!
   * \brief an enumeration describing either the beginning of the time step or
   * the end of the time step
   */
  enum struct TimeStepStage { BEGINNING_OF_TIME_STEP, END_OF_TIME_STEP };
  //! \brief variable associated with the beginning of the time step
  inline constexpr auto bts = TimeStepStage::BEGINNING_OF_TIME_STEP;
  //! \brief variable associated with the end of the time step
  inline constexpr auto ets = TimeStepStage::END_OF_TIME_STEP;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_TIMESTEPSTAGE_HXX */
