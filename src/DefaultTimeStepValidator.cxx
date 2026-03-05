/*!
 * \file   MFEMMGIS/DefaultTimeStepValidator.cxx
 * \brief  This class implements the `DefaultTimeStepValidator` class
 * \date   04/12/2023
 */

#include "MFEMMGIS/DefaultTimeStepValidator.hxx"

namespace mfem_mgis {

  DefaultTimeStepValidator::DefaultTimeStepValidator() noexcept = default;

  std::optional<AbstractTimeStepValidator::Result>
  DefaultTimeStepValidator::validate(Context &ctx) const noexcept {
    return this->callExternalValidators(ctx);
  }  // end of validate

  DefaultTimeStepValidator::~DefaultTimeStepValidator() = default;

}  // end of namespace mfem_mgis