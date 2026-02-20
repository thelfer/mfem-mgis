/*!
 * \file   include/MFEMMGIS/Config.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 */

#ifndef LIB_MFEM_MGIS_CONFIG_HXX
#define LIB_MFEM_MGIS_CONFIG_HXX

#include <limits>
#include <cstdlib>
#include "MGIS/Config.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/InvalidResult.hxx"
#include "MGIS/Utilities/OptionalReference.hxx"

#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

#ifndef MGIS_BEHAVIOUR_API_VERSION
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

#if MGIS_BEHAVIOUR_API_VERSION != 1
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

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
#else  /* defined _WIN32 || defined _WIN64 || defined __CYGWIN__ */
#define MFEM_MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

namespace mfem_mgis {

  using mgis::AbstractErrorHandler;
  using mgis::Context;
  using mgis::isInvalid;
  using mgis::isValid;
  using mgis::OptionalReference;

  namespace attributes {
    //! \brief a simple alias
    using Throwing = ::mgis::attributes::ThrowingAttribute<true>;
    //! \brief a simple alias
    using Unsafe = ::mgis::attributes::UnsafeAttribute;
  }  // namespace attributes
  //
  inline constexpr auto unsafe = ::mgis::attributes::UnsafeAttribute{};
  inline constexpr auto throwing =
      ::mgis::attributes::ThrowingAttribute<true>{};

  //! \brief a simple alias
  using size_type = int;
  /*!
   * \brief constant used to determine if a function view has not a number of
   * components knwon at compile-time size
   */
  inline constexpr size_type dynamic_extent =
      std::numeric_limits<size_type>::max();
  //! \brief alias to the numeric type used
  using real = mgis::real;
  /*!
   * \brief this function can be called to report that the parallel computation
   * are not supported.
   */
  MFEM_MGIS_EXPORT [[noreturn]] void reportUnsupportedParallelComputations();
  //! \brief a simple alias
  using MainFunctionArguments = char**;
  /*!
   * \brief function that must be called to initialize `mfem-mgis`.
   * \param[in] argc: number of arguments
   * \param[in] argv: arguments
   *
   * In parallel, this function calls the `MPI_Init` function.
   * It is safe to call this function multiple time.
   */
  MFEM_MGIS_EXPORT void initialize(int&, MainFunctionArguments&);
  /*!
   * \brief function that must be called to end `mfem-mgis`
   * This call is optional if the code exits normally.
   */
  MFEM_MGIS_EXPORT void finalize();
  //! \return the MPI rank.
  MFEM_MGIS_EXPORT int getMPIrank();
  //! \return the MPI global communicator size.
  MFEM_MGIS_EXPORT int getMPIsize();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception's classes constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   */
  template <typename Exception = std::runtime_error>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception's classes constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception'
   * constructor.
   * \param[in] a: arguments passed to the exception' constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise(Args&&...);
  /*!
   * \brief raise an exception if the first argument is `true`.
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception'
   * constructor.
   * \param[in] b: condition to be checked. If `true`, an exception is
   * thrown.
   * \param[in] a: arguments passed to the exception' constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  MFEM_MGIS_VISIBILITY_LOCAL inline void raise_if(const bool, Args&&...);

  /*!
   * \brief function that must be called if one MPI process detect an fatal
   * error.
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const int = EXIT_FAILURE);
  /*!
   * \brief function that must be called if one MPI process detect an fatal
   * error.
   * \param[in] msg: message displayed of the calling process
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const char* const,
                                           const int = EXIT_FAILURE);
  //! \return if the usage of PETSc has been requested by the user.
  MFEM_MGIS_EXPORT bool usePETSc();
  //! \brief activate PETSc with the configuration file givien in parameter.
  MFEM_MGIS_EXPORT void setPETSc(const char*);
  //! \brief declare default options
  MFEM_MGIS_EXPORT void declareDefaultOptions(mfem::OptionsParser&);

  //! \brief get Output Stream
  MFEM_MGIS_EXPORT std::ostream& getOutputStream();
  //! \brief get Error Stream
  MFEM_MGIS_EXPORT std::ostream& getErrorStream();

}  // namespace mfem_mgis

#include "MFEMMGIS/Config.ixx"

#endif /* LIB_MFEM_MGIS_CONFIG_HXX */
