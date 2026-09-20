/*!
 * \file   include/MFEMMGIS/Utilities/ParametersValidator.hxx
 * \brief  This files declares the `ParametersValidator` class
 * \author Thomas Helfer
 * \date   19/09/2026
 */

#ifndef LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_HXX
#define LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_HXX

#include <set>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  /*!
   * \brief an helper structure used to validate a dictionary of parameters.
   */
  struct MFEM_MGIS_EXPORT ParametersValidator {
    //! \brief arguments of the add methods
    struct AddArguments {
      /*!
       * \brief boolean stating if the added key or added keys are required
       *
       * \note that if a key has already been declared requried, this boolean
       * has no effect.
       */
      const bool required = false;
    };
    static constexpr AddArguments not_required =
        AddArguments{.required = false};
    //! \brief function used to validate a
    using ParameterValidator = std::function<bool(Context&, const Parameter&)>;
    /*!
     * \brief function use to report that the value associated with the given
     * key has not the expected type.
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     */
    static InvalidResult reportUnmatchedTypeError(Context&,
                                                  const std::string&) noexcept;
    /*!
     * \brief function use to report that the value associated with the given
     * key has not one of the expected types.
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     */
    static InvalidResult reportUnmatchedTypesError(Context&,
                                                   const std::string&) noexcept;
    //! \brief constructor
    ParametersValidator() noexcept;
    //! \brief move constructor
    ParametersValidator(ParametersValidator&&) noexcept;
    //! \brief copy constructor
    ParametersValidator(const ParametersValidator&) noexcept;
    //! \brief move assignement
    ParametersValidator& operator=(ParametersValidator&&) noexcept;
    //! \brief standard assignement
    ParametersValidator& operator=(const ParametersValidator&) noexcept;
    /*!
     * \brief add a new allowed key without description.
     *
     * \param[in] k: key
     * \param[in] is_required: state if the given key is required
     *
     * \note the given key already exists, nothing is done.
     */
    ParametersValidator& add(const std::string&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief add a new allowed key with its associated description.
     *
     * \param[in] keys: list of keys
     * \param[in] are_required: state if the given keys are required
     *
     * \note if one of the given key already exists, nothing is done.
     */
    ParametersValidator& add(const std::vector<std::string>&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief add a new allowed key with its associated description.
     *
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] is_required: state if the given key is required
     *
     * \note the given key already exists and the documentation is empty,
     * the given documentation overwrites it. Ortherwise, nothing is done.
     */
    ParametersValidator& add(const std::string&,
                             const std::string&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief add a dictionary associated some allowed keys with their
     * description.
     *
     * \param[in] m: keys and description
     * \param[in] are_required: state if the given keys are required
     *
     * \note the given key already exists and the documentation is empty,
     * the given documentation overwrites it. Ortherwise, nothing is done.
     */
    ParametersValidator& add(const std::map<std::string, std::string>&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief add an arbitrary parameter validator
     *
     * \param[in] k: name of the parameter
     * \param[in] f: validator
     * \param[in] opts: option used to declared the key
     *
     * \note the given name is automatically added to the allowed keys with an
     * empty description if it does not exist
     */
    ParametersValidator& add(const std::string&,
                             const ParameterValidator&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief add an arbitrary parameter validator
     *
     * \param[in] k: name of the parameter
     * \param[in] f: validator
     * \param[in] d: documentation
     * \param[in] opts: option used to declared the key
     *
     * \note the given name is automatically added to the allowed keys with an
     * empty description if it does not exist
     */
    ParametersValidator& add(const std::string&,
                             const std::string&,
                             const ParameterValidator&,
                             const AddArguments& = not_required) noexcept;
    /*!
     * \brief check that the data has one of types given as template arguments
     *
     * \tparam Types: list of allowed types
     * \param[in] k: key
     * \param[in] opts: option used to declared the key
     */
    template <typename... Types>
    ParametersValidator& add(const std::string&,
                             const AddArguments& = not_required) noexcept
        requires((sizeof...(Types) > 0) &&
                 (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data has one of types given as template arguments
     *
     * \tparam Types: list of allowed types
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] opts: option used to declared the key
     */
    template <typename... Types>
    ParametersValidator& add(const std::string&,
                             const std::string&,
                             const AddArguments& = not_required) noexcept
        requires((sizeof...(Types) > 0) &&
                 (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data associated with the given key is a strictly
     * positive integer
     *
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] opts: option used to declared the key
     */
    ParametersValidator& addStrictlyPositiveIntegerCheck(
        const std::string&,
        const std::string&,
        const AddArguments& = not_required) noexcept;
    /*!
     * \brief check that the data associated with the given key is a strictly
     * positive integer
     *
     * \param[in] k: key
     * \param[in] opts: option used to declared the key
     */
    ParametersValidator& addStrictlyPositiveIntegerCheck(
        const std::string&, const AddArguments& = not_required) noexcept;
    /*!
     * \brief declare a list of keys to be incompatible
     *
     * \param[in] k: keys
     *
     * \note given keys are automatically added to the authorized keys
     * \note none of those keys have to be declared required
     */
    ParametersValidator& addKeysIncompatibilityCheck(
        const std::vector<std::string>&) noexcept;
    /*!
     * \brief validate a dictionary of parameters
     *
     * \param[in] throwing: attribute indicating that errors are reported by
     * throwing an exception
     * \param[in] parameters: tested parameters
     *
     * \note This method shall be used in constructors or functions with the
     * throwing attributes when no context is available. If a context is
     * available, it is better to use the following pattern:
     *
     * \code{.cpp}
     * auto or_raise = ctx.getThrowingFailureHandler();
     * ...
     * validator.validate(ctx, parameters) | or_raise;
     * \endcode
     */
    void validate(attributes::Throwing, const Parameters&) const;
    /*!
     * \brief validate a dictionary of parameters
     * \param[in] throwing: attribute indicating that errors are reported by
     * throwing an exception
     * \param[in] parameters: tested parameters
     */
    [[nodiscard]] bool validate(Context&, const Parameters&) const noexcept;
    //! \brief destructor
    ~ParametersValidator();

   private:
    /*!
     * \brief add a key without documentation. If the key already exists, this
     * is a no-op.
     *
     * \param[in] k: key
     * \param[in] opts: option used to declared the key
     */
    void addKey(const std::string&, const AddArguments&) noexcept;
    //! \brief list of incompatible keys
    std::vector<std::vector<std::string>> incompatibilities;
    //! \brief dictionary associating a valid key and its description
    std::map<std::string, std::string> allowed_keys;
    //! \brief list of required keys
    std::set<std::string> required_keys;
    //! \brief validators, sorted by keywords
    std::map<std::string, std::vector<ParameterValidator>> validators;
  };

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Utilities/ParametersValidator.ixx"

#endif /* LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_HXX */
