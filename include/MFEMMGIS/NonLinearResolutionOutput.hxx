/*!
 * \file   include/MFEMMGIS/NonLinearResolutionOutput.hxx
 * \brief  This file declares the `NonLinearResolutionOutput` class
 * \author Thomas Helfer
 * \date   09/12/2021
 */

#ifndef LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX
#define LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX

#include <limits>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/ComputeNextStateOutput.hxx"

namespace mfem_mgis {

  /*!
   * \brief structure storing information about a non linear resolution
   */
  struct NonLinearResolutionOutput {
    /*!
     * \brief status of the resolution:
     * - true if the non linear resolution converged
     * - false if the non linear resolution did not converge
     */
    bool status = false;
    /*!
     * \brief initial norm of the residual
     * \note this may not be defined by every solver (in particular by the ones
     * based on PETSc)
     */
    real initial_residual_norm = std::numeric_limits<real>::quiet_NaN();
    //! \brief final norm of the residual
    real final_residual_norm = std::numeric_limits<real>::quiet_NaN();
    //! \brief number of iterations
    size_type iterations = size_type{};
    //! \brief information collected during the iterations
    std::vector<Parameter> iterations_information;
    /*!
     * \brief conversion operator to a boolean
     * \return the status of the resolution
     */
    inline operator bool() const { return this->status; }
  };  // end of struct NonLinearResolutionOutput

  /*!
   * \brief convert a non linear resolution output
   * \return a set of parameters containing the same information as the given
   * non linear resolution output if the status is valid or an empty optional
   * otherwise.
   *
   * \param[in] output: non linear resolution output
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<ComputeNextStateOutput>
  convertToComputeNextStateOutput(
      const NonLinearResolutionOutput& output) noexcept;

  /*!
   * \brief check if the given non linear resolution output is invalid
   * \param[in] r: non linear resolution output
   * \return if the given non linear resolution output is invalid
   */
  [[nodiscard]] inline bool isInvalid(
      const NonLinearResolutionOutput& r) noexcept {
    return !r;
  }  // end of isInvalid

}  // end of namespace mfem_mgis

namespace mgis::internal {

  //! \brief specialisation for non linear resolution outputs
  template <>
  struct InvalidValueTraits<mfem_mgis::NonLinearResolutionOutput> {
    //! \brief tag indicating that this class is properly specialized
    static constexpr bool isSpecialized = true;
    //! \return an invalid non linear resolution output
    static mfem_mgis::NonLinearResolutionOutput getValue() noexcept {
      return {};
    }
  };

}  // end of namespace mgis::internal

#endif /* LIB_MFEMMGIS_NONLINEARRESOLUTIONOUTPUT_HXX */
