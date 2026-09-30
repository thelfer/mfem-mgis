/*!
 * \file   MFEMMGIS/DefaultTimeStepValidator.hxx
 * \brief  This class declares the `DefaultTimeStepValidator` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEMMGIS_DEFAULTTIMESTEPVALIDATOR_HXX
#define LIB_MFEMMGIS_DEFAULTTIMESTEPVALIDATOR_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepValidatorBase.hxx"

namespace mfem_mgis {

  //! \brief the default time step validator
  struct MFEM_MGIS_EXPORT DefaultTimeStepValidator : TimeStepValidatorBase {
    //! \brief constructor
    DefaultTimeStepValidator() noexcept;
    /*!
     * \return the result of the external validators on success
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] std::optional<Result> validate(
        Context& ctx) const noexcept override;
    //! \brief destructor
    ~DefaultTimeStepValidator() override;
  };  // end of DefaultTimeStepValidator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_DEFAULTTIMESTEPVALIDATOR_HXX */
