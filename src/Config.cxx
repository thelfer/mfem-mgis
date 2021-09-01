/*!
 * \file   src/Config.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#include <cstring>
#include <iostream>
#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif /* MFEM_USE_MPI */
#include "mfem/general/optparser.hpp"
#include "mfem/general/communication.hpp"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /*MFEM_USE_PETSC */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief structure in charge of freeing ressources on exit.
   */
  struct MGIS_VISIBILITY_LOCAL Finalizer {
    //! \brief option used to select the PETSc configuration file
    static const char* const petsc_configuration_file_option;
    //! \return the unique instance of this class
    static Finalizer& get();
    //! \brief initialize the execution of the mfem-mgis
    void initialize(int&, MainFunctionArguments&);
    //! \return true if PETSc is used
    bool usePETSc() const;
    //! \brief activate PETSc and provide a configuration file
    void setPETSc(const char *petscrc_file);
    //! \brief finalize the execution of the mfem-mgis
    void finalize();
    //! \brief abort the process
    [[noreturn]] void abort(int error);

   private:
    //! \brief boolean stating if PETSc shall be used
    bool use_petsc = false;
    //! \brief boolean stating if the finalize method has been called
    bool pendingExit = false;
    /*!
     * \brief constructor
     * \param[in] argc: number of command line arguments
     * \param[in] argv: command line arguments
     */
    Finalizer();
    //! \brief destructor
    ~Finalizer();
  };  // end of Finalizer

  const char* const Finalizer::petsc_configuration_file_option =
      "--petsc-configuration-file";

  Finalizer::Finalizer() = default;

#ifdef MFEM_USE_PETSC
  void Finalizer::setPETSc(const char* petscrc_file) {
    this->use_petsc = true;
    mfem::MFEMInitializePetsc(nullptr, nullptr, petscrc_file, nullptr);
  }
#else /* MFEM_USE_PETSC */
  void Finalizer::setPETSc(const char*) {
    mgis::raise("Finalizer::setPETSc: PETSc is not supported");
  }
#endif /*MFEM_USE_PETSC*/

  void Finalizer::initialize(int& argc, MainFunctionArguments& argv) {
#ifdef MFEM_USE_PETSC
    const char* petscrc_file = nullptr;
#endif /* MFEM_USE_PETSC */
    for (const auto* a = argv; a != argv + argc; ++a) {
#ifdef MFEM_USE_PETSC
      if (std::strcmp(*a, "--use-petsc") == 0) {
       this->use_petsc=true;
      }
      if (std::strcmp(*a, petsc_configuration_file_option) == 0){
        if (petscrc_file != nullptr) {
          mgis::raise("initialize: PETSc configuration file already specified");
        }
        ++a;
        if (a == argv + argc) {
          mgis::raise("initialize: option values missing for --use-petsc");
        }
        petscrc_file = *a;
      }
#endif /* MFEM_USE_PETSC */
    }
#ifdef MFEM_USE_PETSC
    if (this->use_petsc) {
      if (petscrc_file== nullptr) {
        mgis::raise("initialize: no PETSc configuration file given");
      }
      mfem::MFEMInitializePetsc(nullptr, nullptr, petscrc_file, nullptr);
    }
#endif /* MFEM_USE_PETSC */
  }  // end of initialize

  Finalizer& Finalizer::get() {
    static Finalizer f;
    return f;
  }  // end of get

  bool Finalizer::usePETSc() const {
    return this->use_petsc;
  }  // end of usePETSc

  void Finalizer::finalize() {
    if (!this->pendingExit) {
#ifdef MFEM_USE_MPI
#ifdef MFEM_USE_PETSC
      if (this->use_petsc) {
        mfem::MFEMFinalizePetsc();
      }
#endif /* MFEM_USE_PETSC */
      MPI_Finalize();
      this->pendingExit = true;
#endif /* MFEM_USE_MPI */
    }
  }  // end of finalize

  [[noreturn]] void Finalizer::abort(int error) {
    if (!this->pendingExit) {
#ifdef MFEM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, error);
      this->pendingExit = true;
#endif /* MFEM_USE_MPI */
    }
    std::exit(error);
  }  // end of abort

  Finalizer::~Finalizer() { this->finalize(); }

  [[noreturn]] void reportUnsupportedParallelComputations() {
    raise(
        "reportUnsupportedParallelComputations: "
        "unsupported parallel computations");
  }  // end of reportUnsupportedParallelComputations

#ifdef MFEM_USE_MPI

  static void exit_on_failure() {
    try {
      throw;
    } catch (std::exception& e) {
      mfem_mgis::getErrorStream() << e.what() << '\n';
    } catch (...) {
      mfem_mgis::getErrorStream() << "unknown exception thrown";
    }
    abort();
    std::abort();
  }  // end of exit_on_failure

  void initialize(int& argc, MainFunctionArguments& argv) {
    static bool first = true;
    if (first) {
      mgis::setExceptionHandler(exit_on_failure);
      MPI_Init(&argc, &argv);
      if (getMPIrank() != 0) { mfem::out.Disable(); mfem::err.Disable(); } 
      Finalizer::get().initialize(argc, argv);
      first = false;
    }
  }  // end of initialize

#else /* MFEM_USE_MPI */

  void initialize(int&, MainFunctionArguments&) {
    Finalizer::get();
  }  // end of initialize

#endif /* MFEM_USE_MPI */

  void finalize() { Finalizer::get().finalize(); }  // end of finalize

  void abort(const int error) {
    Finalizer::get().abort(error);
    std::exit(error);
  }  // end of abort

  void abort(const char* const msg, const int error) {
    mfem_mgis::getErrorStream() << msg << '\n';
    Finalizer::get().abort(error);
    std::exit(error);
  }  // end of abort

  void declareDefaultOptions(mfem::OptionsParser& parser) {
    static_cast<void>(parser);
#ifdef MFEM_USE_PETSC
    static bool use_petsc = false;
    static const char* petscrc_file = "";
#endif /* MFEM_USE_PETSC */
#ifdef MFEM_USE_PETSC 
    parser.AddOption(&use_petsc, "-up", "--use-petsc", "-no-up", "--donot-use-petsc",
                     "Use or not PETSc to solve the nonlinear system.");
    parser.AddOption(&petscrc_file, "petscrc_file",
                     Finalizer::petsc_configuration_file_option,
                     "Path to the PETSc configuration file.");
#endif /* MFEM_USE_PETSC */
}

  bool usePETSc() { return Finalizer::get().usePETSc(); }  // end of usePETSc

  void setPETSc(const char* petscrc_file) { Finalizer::get().setPETSc(petscrc_file); }  // end of setPETSc

}  // namespace mfem_mgis
