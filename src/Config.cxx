/*!
 * \file   src/Config.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#include <iostream>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Config.hxx"
#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif /* MFEM_USE_MPI */

namespace mfem_mgis {

  /*!
   *
   */
  struct MGIS_VISIBILITY_LOCAL Finalizer {
    //! \return the unique instance of this class
    static Finalizer& get();
    //! \brief finalize the execution of the mfem-mgis
    void finalize();

   private:
    //! boolean stating if the finalize method has been called
    bool hasFinalizeBeenCalled = false;
    /*!
     * \brief constructor
     * \param[in] argc: number of command line arguments
     * \param[in] argv: command line arguments
     */
    Finalizer();
    //! \brief destructor
    ~Finalizer();
  };  // end of Finalizer

  Finalizer& Finalizer::get() { static Finalizer f;
    return f;
  }  // end of get

  Finalizer::Finalizer() = default;

  void Finalizer::finalize() {
    if (!this->hasFinalizeBeenCalled) {
#ifdef MFEM_USE_MPI
      MPI_Finalize();
      this->hasFinalizeBeenCalled = true;
#endif /* MFEM_USE_MPI */
    }
  }  // end of finalize

  Finalizer::~Finalizer() {
    this->finalize();
  }

  [[noreturn]] void
  reportUnsupportedParallelComputations() {
    mgis::raise(
        "reportUnsupportedParallelComputations: "
        "unsupported parallel computations");
  }  // end of reportUnsupportedParallelComputations

  static void exit_on_failure() {
    try {
      throw;
    } catch (std::exception& e) {
      std::cerr << e.what() << '\n';
    } catch (...) {
      std::cerr << "unknown exception thrown";
    }
    finalize();
    std::abort();
  } // end of exit_on_failure

#ifdef MFEM_USE_MPI

  void initialize(int& argc, MainFunctionArguments& argv){
    static bool first = true;
    if (first) {
      mgis::setExceptionHandler(exit_on_failure);
      MPI_Init(&argc, &argv);
      Finalizer::get();
      first = false;
    }
  }  // end of initialize

#else /* MFEM_USE_MPI */

  void initialize(int&, MainFunctionArguments& ){
    Finalizer::get();
  }  // end of initialize

#endif /* MFEM_USE_MPI */

  void finalize() { Finalizer::get().finalize(); }  // end of finalize

}  // namespace mfem_mgis
