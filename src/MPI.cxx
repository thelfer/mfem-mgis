/*!
 * \file   src/MPI.cxx
 * \brief  This file implements the functions declared in `MFEMMGIS/MPI.hxx`
 * \author Thomas Helfer
 * \date   25/02/2026
 */

#include "mfem/config/config.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/MPI.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI
  bool isTrueOnAllProcesses(const MPI_Comm& c, const bool b) noexcept {
    auto r = b;
    MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_CXX_BOOL, MPI_LAND, c);
    return r;
  }
#endif /* MFEM_USE_MPI */

  bool isTrueOnAllProcesses(const MeshDiscretization& m,
                            const bool b) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return isTrueOnAllProcesses(getMPICommunicator(m), b);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return b;
  }  // end of isTrueOnAllProcesses

}  // end of namespace mfem_mgis
