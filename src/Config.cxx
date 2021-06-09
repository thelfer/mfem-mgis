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
    //! \return the unique instance of this class
    static Finalizer& get();
    //! \brief initialize the execution of the mfem-mgis
    void initialize(int&, MainFunctionArguments&);
    //! \return true if PETSc is used
    bool usePETSc() const;
    //! \brief initiate a parser for command line arguments
    std::shared_ptr<mfem::OptionsParser>  beginParser();
    //! \brief close the parser zone
    void endParser();
    //! \brief finalize the execution of the mfem-mgis
    void finalize();
    //! \brief abort the process
    [[noreturn]] void abort(int error);

   private:
    //! \brief boolean stating if PETSc shall be used
    bool use_petsc = false;
    //! \brief configuration file for PETSc
    const char* petscrc_file = "";
    //! \brief boolean stating if the finalize method has been called
    bool pendingExit = false;
    //! \brief all options passed one the command line 
    std::shared_ptr<mfem::OptionsParser> opt_parser;
    /*!
     * \brief constructor
     * \param[in] argc: number of command line arguments
     * \param[in] argv: command line arguments
     */
    Finalizer();
    //! \brief destructor
    ~Finalizer();
  };  // end of Finalizer

  Finalizer::Finalizer() = default;

  void Finalizer::initialize(int& argc, MainFunctionArguments& argv) {
    opt_parser = std::make_shared<mfem::OptionsParser>(argc, argv);
#ifdef MFEM_USE_PETSC
    opt_parser->AddOption(&this->use_petsc, "-up", "--use-petsc", "-no-up",
                  "--no-use-petsc", "Activate PETSC support.");
    opt_parser->AddOption(&this->petscrc_file, "-pcf", "--petsc-configuration-file",
		   "Petsc configuration file");
#endif /* MFEM_USE_PETSC */
  }  // end of initialize

  std::shared_ptr<mfem::OptionsParser> Finalizer::beginParser() {
    return opt_parser;
  }  // end of beginParser

  void Finalizer::endParser() {
#ifdef MFEM_USE_PETSC
    if (this->use_petsc) {
      MFEM_VERIFY(this->petscrc_file != "", 
		  "initialize: no PETSc configuration file given");
      mfem::MFEMInitializePetsc(nullptr, nullptr, this->petscrc_file, nullptr);
    }
#endif /* MFEM_USE_PETSC */
  }  // end of endParser

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
      MPI_Finalize();
#ifdef MFEM_USE_PETSC
      if (this->use_petsc) {
        mfem::MFEMFinalizePetsc();
      }
#endif /* MFEM_USE_PETSC */
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
      std::cerr << e.what() << '\n';
    } catch (...) {
      std::cerr << "unknown exception thrown";
    }
    abort();
    std::abort();
  }  // end of exit_on_failure

  void initialize(int& argc, MainFunctionArguments& argv) {
    static bool first = true;
    if (first) {
      mgis::setExceptionHandler(exit_on_failure);
      MPI_Init(&argc, &argv);
      Finalizer::get().initialize(argc, argv);
      first = false;
    }      
  }  // end of initialize

#else /* MFEM_USE_MPI */

  void initialize(int&, MainFunctionArguments&) {
    static bool first = true;
    if (first) {
      Finalizer::get().initialize(argc, argv);
      first = false;
    }      
  }  // end of initialize

#endif /* MFEM_USE_MPI */

  std::shared_ptr<mfem::OptionsParser> beginParser() {
    auto opt_parser = Finalizer::get().beginParser();
    MFEM_VERIFY(opt_parser.get() != nullptr, "beginParser called before initialize");
    return(opt_parser);
  }

  void endParser() {
    Finalizer::get().endParser();
  }

  void finalize() { Finalizer::get().finalize(); }  // end of finalize

  void abort(const int error) {
    Finalizer::get().abort(error);
    std::exit(error);
  }  // end of abort

  void abort(const char* const msg, const int error) {
    std::cerr << msg << '\n';
    Finalizer::get().abort(error);
    std::exit(error);
  }  // end of abort

  bool usePETSc() { return Finalizer::get().usePETSc(); }  // end of usePETSc

}  // namespace mfem_mgis
