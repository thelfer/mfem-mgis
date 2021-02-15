/*!
 * \file   include/MFEMMGIS/Config.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   19/06/2018
 */

#ifndef LIB_MFEM_MGIS_CONFIG_HXX
#define LIB_MFEM_MGIS_CONFIG_HXX

#include "MGIS/Config.hxx"

#define MFEM_MGIS_VISIBILITY_LOCAL MGIS_VISIBILITY_LOCAL

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFEMMGIS_EXPORTS
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#else /* defined MFEMMGIS_EXPORTS */
#ifndef MFEM_MGIS_STATIC_BUILD
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_IMPORT
#else /* MFEM_MGIS_STATIC_BUILD */
#define MFEM_MGIS_EXPORT
#endif /* MFEM_MGIS_STATIC_BUILD */
#endif /* defined MFEMMGIS_EXPORTS */
#else /* defined _WIN32 || defined _WIN64 || defined __CYGWIN__ */
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

namespace mfem_mgis {

  //! a simple alias
  using size_type = int;

  //! alias to the numeric type used
  using real = mgis::real;

  /*!
   * \brief this function can be called to report that the parallel computation
   * are not supported.
   */
  MFEM_MGIS_EXPORT [[noreturn]] void reportUnsupportedParallelComputations();

  /*!
   * \brief this function can be called to report that the sequential computation
   * are not supported.
   */
  MFEM_MGIS_EXPORT [[noreturn]] void reportUnsupportedSerialComputations();

}  // end of namespace mfront_mgis

#endif /* LIB_MFEM_MGIS_CONFIG_HXX */
