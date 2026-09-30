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
#include "MGIS/LogStream.hxx"
#include "MGIS/InvalidResult.hxx"
#include "MGIS/Utilities/Construct.hxx"
#include "MGIS/Utilities/Invoke.hxx"
#include "MGIS/Utilities/OptionalReference.hxx"

#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

#ifndef MGIS_BEHAVIOUR_API_VERSION
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

#if MGIS_BEHAVIOUR_API_VERSION != 1
#error "Incompatible version of MGIS"
#endif /* MGIS_BEHAVIOUR_API_VERSION */

//! \brief macro hiding a symbol from the shared library
#define MFEM_MGIS_VISIBILITY_LOCAL MGIS_VISIBILITY_LOCAL

//! \brief macro exporting a symbol from the shared library
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

  using mgis::abort;
  using mgis::AbstractErrorHandler;
  using mgis::areInvalid;
  using mgis::areValid;
  using mgis::construct;
  using mgis::Context;
  using mgis::getDefaultVerbosityLevel;
  using mgis::InvalidResult;
  using mgis::invoke;
  using mgis::isInvalid;
  using mgis::isValid;
  using mgis::make_shared;
  using mgis::make_shared_as;
  using mgis::make_unique;
  using mgis::make_unique_as;
  using mgis::OptionalReference;
  using mgis::registerExceptionInErrorBacktrace;
  using mgis::terminate;
  using mgis::verboseDebug;
  using mgis::verboseFull;
  using mgis::verboseLevel0;
  using mgis::verboseLevel1;
  using mgis::verboseLevel2;
  using mgis::verboseLevel3;
  using mgis::VerbosityLevel;

  using mgis::debug;
  using mgis::getDefaultLogStream;
  using mgis::setDefaultLogStream;
  using mgis::warning;

  namespace attributes {
    //! \brief a simple alias
    using Throwing = ::mgis::attributes::ThrowingAttribute<true>;
    //! \brief a simple alias
    using MayThrow = ::mgis::attributes::ThrowingAttribute<true>;
    //! \brief a simple alias
    using MayAbort = ::mgis::attributes::AbortingAttribute<true>;
    //! \brief a simple alias
    using Unsafe = ::mgis::attributes::UnsafeAttribute;
  }  // namespace attributes
  //! \brief tag marking a call as unsafe without precautions
  inline constexpr auto unsafe = ::mgis::attributes::UnsafeAttribute{};
  //! \brief deprecated, use may_throw
  [[deprecated]] inline constexpr auto throwing =
      ::mgis::attributes::ThrowingAttribute<true>{};
  //! \brief tag marking a call that may throw
  inline constexpr auto may_throw =
      ::mgis::attributes::ThrowingAttribute<true>{};
  //! \brief tag marking a call that may abort
  inline constexpr auto may_abort =
      ::mgis::attributes::AbortingAttribute<true>{};

  //! \brief a simple alias
  using size_type = int;
  /*!
   * \brief constant stating that a number of components is not known at
   * compile-time
   */
  inline constexpr size_type dynamic_extent =
      std::numeric_limits<size_type>::max();
  //! \brief alias to the numeric type used
  using real = mgis::real;
  /*!
   * \brief this function can be called to report that the parallel
   * computations are not supported.
   */
  MFEM_MGIS_EXPORT [[noreturn]] void reportUnsupportedParallelComputations();
  //! \brief a simple alias
  using MainFunctionArguments = char**;
  /*!
   * \brief function that must be called to initialize `mfem-mgis`.
   * \param[in, out] argc: number of arguments
   * \param[in, out] argv: arguments
   *
   * In parallel, this function calls the `MPI_Init` function.
   * It is safe to call this function multiple times.
   */
  MFEM_MGIS_EXPORT void initialize(int& argc, MainFunctionArguments& argv);
  /*!
   * \brief function that must be called to end `mfem-mgis`.
   * This call is optional if the code exits normally.
   */
  MFEM_MGIS_EXPORT void finalize();
  //! \return the MPI rank for the default communicator
  MFEM_MGIS_EXPORT [[deprecated]] int getMPIrank();
  //! \return the total number of MPI processes for the default communicator.
  MFEM_MGIS_EXPORT [[deprecated]] int getMPIsize();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception classes' constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   */
  template <typename Exception = std::runtime_error>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception classes' constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception's
   * constructor.
   * \param[in] a: arguments passed to the exception's constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  [[noreturn]] MFEM_MGIS_VISIBILITY_LOCAL void raise(Args&&... a);
  /*!
   * \brief raise an exception if the first argument is `true`.
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception's
   * constructor.
   * \param[in] b: condition to be checked. If `true`, an exception is
   * thrown.
   * \param[in] a: arguments passed to the exception's constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  MFEM_MGIS_VISIBILITY_LOCAL inline void raise_if(const bool b, Args&&... a);

  /*!
   * \brief function that must be called if one MPI process detects a fatal
   * error.
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const int error = EXIT_FAILURE);
  /*!
   * \brief function that must be called if one MPI process detects a fatal
   * error.
   * \param[in] msg: message displayed by the calling process
   * \param[in] error: exit status
   */
  MFEM_MGIS_EXPORT [[noreturn]] void abort(const char* const msg,
                                           const int error = EXIT_FAILURE);
  //! \return if the usage of PETSc has been requested by the user.
  MFEM_MGIS_EXPORT bool usePETSc();
  /*!
   * \brief activate PETSc with the configuration file given in parameter.
   * \param[in] petscrc_file: PETSc configuration file
   */
  MFEM_MGIS_EXPORT void setPETSc(const char* petscrc_file);
  /*!
   * \brief declare default options. Only the PETSc options are declared,
   * if PETSc is available.
   * \param[in, out] parser: options parser
   */
  MFEM_MGIS_EXPORT void declareDefaultOptions(mfem::OptionsParser& parser);

  //! \return the output stream
  MFEM_MGIS_EXPORT std::ostream& getOutputStream();
  //! \return the error stream
  MFEM_MGIS_EXPORT std::ostream& getErrorStream();

}  // namespace mfem_mgis

#include "MFEMMGIS/Config.ixx"

#endif /* LIB_MFEM_MGIS_CONFIG_HXX */
