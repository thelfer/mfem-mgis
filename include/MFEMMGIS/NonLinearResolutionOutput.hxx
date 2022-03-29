/*!
 * \file   include/MFEMMGIS/NonLinearResolutionOutput.hxx
 * \brief
 * \author Thomas Helfer
 * \date   09/12/2021
 */

#ifndef LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX
#define LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX

#include <limits>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief structure storing information about a non linear resolution
   */
  struct NonLinearResolutionOutput {
    /*!
     * \brief status of the resolution:
     * - true if the non linear resolution converged
     * - false if the non linear resolution did not converged
     */
    bool status = false;
    /*!
     * \brief initial norm of the residual
     * \note: this may not be defined by every solver (in particular by the ones
     * based on PETCs)
     */
    real initial_residual_norm = std::numeric_limits<real>::quiet_NaN();
    //! \brief final norm of the residual
    real final_residual_norm = std::numeric_limits<real>::quiet_NaN();
    //! \brief number of iterations
    size_type iterations = size_type{};
    //! \brief convertion operator to a boolean
    inline operator bool() const { return this->status; }
  };  // end of struct NonLinearResolutionOutput

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX */
