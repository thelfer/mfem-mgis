/*!
 * \file   MFEMMGIS/TimeStepValidatorBase.hxx
 * \brief  This class declares the `TimeStepValidatorBase` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEMMGIS_TIMESTEPVALIDATORBASE_HXX
#define LIB_MFEMMGIS_TIMESTEPVALIDATORBASE_HXX

#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractTimeStepValidator.hxx"

namespace mfem_mgis {

  //! \brief a common class for most time step validators
  struct MFEM_MGIS_EXPORT TimeStepValidatorBase : AbstractTimeStepValidator {
    //! \brief constructor
    TimeStepValidatorBase() noexcept;
    //
    void addValidator(const std::string_view,
                      const ExternalValidator &) noexcept override;
    void addValidator(const ExternalValidator &) noexcept override;
    //! \brief destructor
    ~TimeStepValidatorBase() override;

   protected:
    /*!
     * \brief call the external validators
     * \param[in] ctx: execution context
     */
    std::optional<Result> callExternalValidators(Context &) const noexcept;
    //! \brief registered external validators
    std::vector<std::pair<std::string, ExternalValidator>> externalValidators;
  };  // end of struct TimeStepValidatorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_TIMESTEPVALIDATORBASE_HXX */
