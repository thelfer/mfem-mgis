/*!
 * \file   MFEMMGIS/MPI.ixx
 * \brief
 * \author Thomas Helfer
 * \date   06/02/2026
 */

#ifndef LIB_MFEM_MGIS_MPI_IXX
#define LIB_MFEM_MGIS_MPI_IXX

namespace mfem_mgis {

  template <typename T>
  bool isValidOnAllProcesses(const FiniteElementDiscretization& fed,
                             const T& v) noexcept {
    return isTrueOnAllProcesses(fed, isValid(v));
  }  // end of isValidOnAllProcesses

}  // end of namespace mfem_mgis

#include "MFEMMGIS/MPI.ixx"

#endif /* LIB_MFEM_MGIS_MPI_IXX */
