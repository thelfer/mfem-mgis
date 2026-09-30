/*!
 * \file   MFEMMGIS/TimeStepValidatorBase.hxx
 * \brief  This file declares the `TimeStepValidatorBase` class
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
    /*!
     * \brief add an external validator
     * \param[in] n: name of the external validator
     * \param[in] v: external validator
     */
    void addValidator(const std::string_view n,
                      const ExternalValidator& v) noexcept override;
    /*!
     * \brief add an external validator with a default name
     * \param[in] v: external validator
     */
    void addValidator(const ExternalValidator& v) noexcept override;
    //! \brief destructor
    ~TimeStepValidatorBase() override;

   protected:
    /*!
     * \brief call the external validators
     * \param[in, out] ctx: execution context
     * \return the combined result of the external validators on success
     */
    std::optional<Result> callExternalValidators(Context& ctx) const noexcept;
    //! \brief registered external validators
    std::vector<std::pair<std::string, ExternalValidator>> externalValidators;
  };  // end of struct TimeStepValidatorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_TIMESTEPVALIDATORBASE_HXX */
