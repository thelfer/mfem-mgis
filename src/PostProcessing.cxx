/*!
 * \file   src/PostProcessing.cxx
 * \brief
 * \author Thomas Helfer
 * \date   08/03/2021
 */

#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  PostProcessing<true>::~PostProcessing() = default;

#endif /* MFEM_USE_MPI */

  PostProcessing<false>::~PostProcessing() = default;

}  // end of namespace mfem_mgis
