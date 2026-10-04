/*!
 * \file   MFEMMGIS/ComputeNextStateOutput.hxx
 * \brief  This file declares the ComputeNextStateOutput class
 * \date   28/08/2024
 */

#ifndef LIB_MFEM_MGIS_COMPUTE_NEXT_STATE_OUTPUT_HXX
#define LIB_MFEM_MGIS_COMPUTE_NEXT_STATE_OUTPUT_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  /*!
   * \brief a structure returned by the computeNextState method.
   *
   * This structure is mostly a wrapper around the `Parameters`
   * class.
   *
   * Coupling items are free to populate this structure
   * with any relevant information
   */
  struct MFEM_MGIS_EXPORT ComputeNextStateOutput : public Parameters {
    using Parameters::Parameters;
    using Parameters::operator=;
    //! \brief default constructor
    ComputeNextStateOutput() noexcept;
    //! \brief move constructor
    ComputeNextStateOutput(ComputeNextStateOutput &&) noexcept;
    //! \brief copy constructor
    ComputeNextStateOutput(const ComputeNextStateOutput &) noexcept;
    //! \brief move assignment
    ComputeNextStateOutput &operator=(ComputeNextStateOutput &&) noexcept;
    //! \brief copy assignment
    ComputeNextStateOutput &operator=(const ComputeNextStateOutput &) noexcept;
    //! \brief destructor
    ~ComputeNextStateOutput() noexcept;
  };  // end of ComputeNextStateOutput

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_COMPUTE_NEXT_STATE_OUTPUT_HXX */
