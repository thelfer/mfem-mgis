/*!
 * \file   include/MFEMMGIS/Config.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 */

#ifndef LIB_MFEM_MGIS_CONFIG_IXX
#define LIB_MFEM_MGIS_CONFIG_IXX

#include <utility>
#include "MGIS/Raise.hxx"
#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif /* MFEM_USE_MPI */

namespace mfem_mgis {

  template <typename Exception>
  inline void raise() {
    mgis::raise<Exception>();
  }  // end of raise

  template <typename Exception, typename... Args>
  inline void raise(Args&&... args) {
    mgis::raise<Exception>(std::forward<Args...>(args...));
  }  // end of raise

  template <typename Exception>
  void raise_if(const bool c) {
    if (c) {
      raise<Exception>();
    }
  }  // end of raise

  template <typename Exception, typename... Args>
  void raise_if(const bool c, Args&&... a) {
    if (c) {
      raise<Exception>(std::forward<Args...>(a...));
    }
  }  // end of raise

  inline int getMPIrank() {
    int rank = 0;
#ifdef MFEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
    return(0);
  } //end of getMPIrank

  inline int getMPIsize() {
    int size=1;
#ifdef MFEM_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif /* MFEM_USE_MPI */
    return(size);
  } //end of getMPIsize
  
}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_CONFIG_IXX */
