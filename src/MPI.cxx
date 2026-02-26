/*!
 * \file   src/MPI.cxx
 * \brief
 * \author Thomas Helfer
 * \date   25/02/2026
 */

#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MPI.hxx"

namespace mfem_mgis {

  bool isTrueOnAllProcesses(const FiniteElementDiscretization& fed,
                            const bool b) noexcept {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto r = b;
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_C_BOOL, MPI_LAND,
                    fed.getFiniteElementSpace<true>().GetComm());
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return b;
  }  // end of isTrueOnAllProcesses

}  // end of namespace mfem_mgis
