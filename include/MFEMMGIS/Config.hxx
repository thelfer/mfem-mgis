/*!
 * \file   include/MFEMMGIS/Config.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   19/06/2018
 */

#ifndef LIB_MFEM_MGIS_CONFIG_HXX
#define LIB_MFEM_MGIS_CONFIG_HXX

#include <cstddef>

/*!
 * Macro extracted from :
 * "Why is the new C++ visibility support so useful?"
 * from http://gcc.gnu.org/wiki/Visibility
 */
#if defined _WIN32 || defined _WIN64 ||defined __CYGWIN__
#define MFEM_MGIS_VISIBILITY_IMPORT __declspec(dllimport)
#define MFEM_MGIS_VISIBILITY_EXPORT __declspec(dllexport)
#define MFEM_MGIS_VISIBILITY_LOCAL
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (! defined __INTEL_COMPILER)
#if __GNUC__ >= 4
#define MFEM_MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MFEM_MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MFEM_MGIS_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* __GNUC__ >= 4 */
#define MFEM_MGIS_VISIBILITY_IMPORT
#define MFEM_MGIS_VISIBILITY_EXPORT
#define MFEM_MGIS_VISIBILITY_LOCAL
#endif /* __GNUC__ >= 4 */
#elif defined __INTEL_COMPILER
#define MFEM_MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MFEM_MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MFEM_MGIS_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* defined __INTEL_COMPILER */
#define MFEM_MGIS_VISIBILITY_IMPORT
#define MFEM_MGIS_VISIBILITY_EXPORT
#define MFEM_MGIS_VISIBILITY_LOCAL
#endif /* defined __INTEL_COMPILER */
#endif /* defined _WIN32 || defined _WIN64 ||defined __CYGWIN__ */

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFEMMGIS_EXPORTS
#define MFEM_MGIS_EXPORT MFEM_MGIS_VISIBILITY_EXPORT
#else
#ifndef MFEM_MGIS_STATIC_BUILD
#define MFEM_MGIS_EXPORT MFEM_MGIS_VISIBILITY_IMPORT
#else
#define MFEM_MGIS_EXPORT
#endif
#endif
#else
#define MFEM_MGIS_EXPORT MFEM_MGIS_VISIBILITY_EXPORT
#endif /* */

namespace mfem_mgis {

  //! a simple alias
  using size_type = size_t;

  //! alias to the numeric type used
  using real = double;

}  // end of namespace mfront_mgis

#endif /* LIB_MFEM_MGIS_CONFIG_HXX */
