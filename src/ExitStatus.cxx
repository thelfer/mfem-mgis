/*!
 * \file   ExitStatus.cxx
 * \brief
 * \author Thomas Helfer
 * \date   15/03/2026
 */

#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/ExitStatus.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  void ExitStatus::synchronize(const MPI_Comm c) noexcept {
    MPI_Allreduce(MPI_IN_PLACE, &(this->status.value), 1, mpi_type<size_type>,
                  MPI_SUM, c);
  }  // end of synchronize

  ExitStatus synchronize(const ExitStatus s, const MPI_Comm c) noexcept {
    auto ns = s;
    ns.synchronize(c);
    return ns;
  }  // end of synchronize

#endif /* MFEM_USE_MPI */

}  // end of namespace mfem_mgis