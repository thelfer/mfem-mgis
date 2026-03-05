/*!
 * \file   MFEMMGIS/AbstractTimeStepValidator.hxx
 * \brief  This class declares the `AbstractTimeStepValidator` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEMMGIS_ABSTRACTTIMESTEPVALIDATOR_HXX
#define LIB_MFEMMGIS_ABSTRACTTIMESTEPVALIDATOR_HXX

#include <vector>
#include <limits>
#include <string>
#include <utility>
#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief class used to validate a new time step after
   * the convergence of the coupling scheme.
   */
  struct MFEM_MGIS_EXPORT AbstractTimeStepValidator {
    //! \brief a simple alias
    using ExternalValidator = std::function<std::pair<bool, real>()>;
    //! \brief structure returned by the validate method
    struct [[nodiscard]] Result {
      //! \brief boolean stating if the time step is valid or shall be rejected
      bool isValid = true;
      /*!
       * \brief the recommended time increment if the time step is rejected
       *
       * \note this output shall only be considered if the time step is rejected
       */
      real recommendedTimeIncrement = std::numeric_limits<real>::max();
      /*!
       * \brief reasons why the time step has been rejected.
       *
       * \note this output shall only be considered if the time step is rejected
       */
      std::vector<std::string> reasons;
    };
    /*!
     * \brief add an external validator
     * \param[in] n: name of the external validator
     * \param[in] v: external validator
     */
    virtual void addValidator(std::string_view,
                              const ExternalValidator &) noexcept = 0;
    /*!
     * \brief add an external validator
     * \param[in] v: external validator
     */
    virtual void addValidator(const ExternalValidator &) noexcept = 0;
    /*!
     * \return a pair on success. The first member states if the time step is
     * valid. The second member is an estimate of a better time step if the time
     * step is rejected.
     *
     * \param[in] ctx: execution context
     */
    [[nodiscard]] virtual std::optional<Result> validate(
        Context &) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractTimeStepValidator();
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTTIMESTEPVALIDATOR_HXX */
