/*!
 * \file   src/PostProcessingBase.cxx
 * \brief  This file implements the `PostProcessingBase` class.
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing/PostProcessingBase.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  PostProcessingBase::getParametersDescription() noexcept {
    return {{"allTimeSteps",
             "boolean stating if the fields must be written at each time step "
             "or only at the end time steps marked by the user as  "
             "post-processing times (default option)"}};
  }  // end of getParametersDescription

  PostProcessingBase::PostProcessingBase(PhysicalSystem &ps,
                                         const Parameters &params,
                                         const bool defaultAllTimeStepsValue)
      : physicalSystem(ps),
        allTimeSteps(get_if<bool>(
            throwing, params, "allTimeSteps", defaultAllTimeStepsValue)) {}

  PhysicalSystem &PostProcessingBase::getPhysicalSystem() noexcept {
    return this->physicalSystem;
  }

  const PhysicalSystem &PostProcessingBase::getPhysicalSystem() const noexcept {
    return this->physicalSystem;
  }

  PostProcessingBase::~PostProcessingBase() noexcept = default;

}  // namespace mfem_mgis