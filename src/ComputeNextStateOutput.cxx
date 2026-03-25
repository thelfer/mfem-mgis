/*!
 * \file   src/ComputeNextStateOutput.cxx
 * \brief  This file implements the ComputeNextStateOutput class
 * \date   28/08/2024
 */

#include "MFEMMGIS/ComputeNextStateOutput.hxx"

namespace mfem_mgis {

  ComputeNextStateOutput::ComputeNextStateOutput() noexcept = default;
  ComputeNextStateOutput::ComputeNextStateOutput(
      ComputeNextStateOutput &&) noexcept = default;
  ComputeNextStateOutput::ComputeNextStateOutput(
      const ComputeNextStateOutput &) noexcept = default;
  ComputeNextStateOutput &ComputeNextStateOutput::operator=(
      ComputeNextStateOutput &&) noexcept = default;
  ComputeNextStateOutput &ComputeNextStateOutput::operator=(
      const ComputeNextStateOutput &) noexcept = default;
  ComputeNextStateOutput::~ComputeNextStateOutput() noexcept = default;

}  // end of namespace mfem_mgis
