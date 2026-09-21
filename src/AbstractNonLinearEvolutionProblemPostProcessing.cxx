/*!
 * \file   src/AbstractNonLinearEvolutionProblemPostProcessing.cxx
 * \brief
 * \author Thomas Helfer
 * \date   08/03/2021
 */

#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  AbstractNonLinearEvolutionProblemPostProcessing<
      true>::~AbstractNonLinearEvolutionProblemPostProcessing() = default;

#endif /* MFEM_USE_MPI */

  AbstractNonLinearEvolutionProblemPostProcessing<
      false>::~AbstractNonLinearEvolutionProblemPostProcessing() = default;

}  // end of namespace mfem_mgis
