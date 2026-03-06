/*!
 * \file   TimeStep.hxx
 * \brief  This file introduces the TimeStep data structure
 * \author Thomas Helfer
 * \date   05/03/2026
 */

#ifndef LIB_MFEM_MGIS_TIMESTEP_HXX
#define LIB_MFEM_MGIS_TIMESTEP_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief data structure associated with at time step
   */
  struct TimeStep {
    //! \brief beginning of the time step
    real begin;
    //! \brief end of the time step
    real end;
    //! \brief time increment
    real dt;
  };  // end of struct TimeStep

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_TIMESTEP_HXX */
