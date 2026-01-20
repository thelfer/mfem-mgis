/*!
 * \file   src/LinearSolverHandler.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   20/01/2026
 */
#include "MFEMMGIS/LinearSolverHandler.hxx"

namespace mfem_mgis {

  bool isInvalid(const LinearSolverHandler& s) noexcept {
    return isInvalid(s.linear_solver);
  }  // end of isInvalid

} // end of namespace mfem_mgis
