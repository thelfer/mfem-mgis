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
   * \brief a helper structure used to validate a dictionary of parameters.
   * \note this class provides a convenient way to declare and validate
   * parameters in a structured manner.
   */
  struct MFEM_MGIS_EXPORT ParametersValidator {
    /*!
     * \brief arguments of the `add` methods
     */
    struct AddArguments {
      /*!
       * \brief boolean stating if the added key or added keys are required
       * \note if a key has already been declared required, this boolean
       * has no effect.
       */
      const bool required = false;
    };
    //! \brief arguments for keys which are not required
    static constexpr AddArguments not_required =
        AddArguments{.required = false};
    /*!
     * \brief function used to validate a parameter
     */
    using ParameterValidator = std::function<bool(Context&, const Parameter&)>;
    /*!
     * \brief report that the value associated with the given key does not
     * have the expected type.
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     * \return an invalid result
     */
    static InvalidResult reportUnmatchedTypeError(
        Context& ctx, const std::string& k) noexcept;
    /*!
     * \brief report that the value associated with the given key does not
     * have one of the expected types.
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     * \return an invalid result
     */
    static InvalidResult reportUnmatchedTypesError(
        Context& ctx, const std::string& k) noexcept;
    //! \brief default constructor
    ParametersValidator() noexcept;
    //! \brief move constructor
    ParametersValidator(ParametersValidator&&) noexcept;
    //! \brief copy constructor
    ParametersValidator(const ParametersValidator&) noexcept;
    //! \brief move assignment
    ParametersValidator& operator=(ParametersValidator&&) noexcept;
    //! \brief Standard assignment
    ParametersValidator& operator=(const ParametersValidator&) noexcept;
    /*!
     * \brief add a new allowed key without description.
     * \param[in] k: key
     * \param[in] opts: state if the given key is required
     * \return the modified validator
     * \note if the given key already exists, its description is left
     * unchanged.
     */
    ParametersValidator& add(const std::string& k,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief add a list of allowed keys without description.
     * \param[in] keys: list of keys
     * \param[in] opts: state if the given keys are required
     * \return the modified validator
     * \note the descriptions of the existing keys are left unchanged.
     */
    ParametersValidator& add(const std::vector<std::string>& keys,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief add a new allowed key with its associated description.
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] opts: state if the given key is required
     * \return the modified validator
     * \note if the given key already exists and the documentation is empty,
     * the given documentation overwrites it. Otherwise, the documentation is
     * left unchanged.
     */
    ParametersValidator& add(const std::string& k,
                             const std::string& d,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief add a dictionary associating some allowed keys with their
     * description.
     * \param[in] m: keys and associated descriptions
     * \param[in] opts: state if the given keys are required
     * \return the modified validator
     * \note if the given key already exists and the documentation is empty,
     * the given documentation overwrites it. Otherwise, the documentation is
     * left unchanged.
     */
    ParametersValidator& add(const std::map<std::string, std::string>& m,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief add an arbitrary parameter validator
     * \param[in] k: name of the parameter
     * \param[in] f: validator
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     * \note the given name is automatically added to the allowed keys with an
     * empty description if it does not exist
     */
    ParametersValidator& add(const std::string& k,
                             const ParameterValidator& f,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief add an arbitrary parameter validator
     * \param[in] k: name of the parameter
     * \param[in] d: documentation
     * \param[in] f: validator
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     * \note the given name is automatically added to the allowed keys if it
     * does not exist
     */
    ParametersValidator& add(const std::string& k,
                             const std::string& d,
                             const ParameterValidator& f,
                             const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief check that the data has one of the types given as template
     * arguments
     * \tparam Types: list of allowed types
     * \param[in] k: key
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    template <typename... Types>
    ParametersValidator& add(const std::string& k,
                             const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data has one of the types given as template
     * arguments
     * \tparam Types: list of allowed types
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    template <typename... Types>
    ParametersValidator& add(const std::string& k,
                             const std::string& d,
                             const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data has one of the types given as template
     * arguments
     * \tparam Types: list of allowed types
     * \param[in] keys: list of keys
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    template <typename... Types>
    ParametersValidator& add(const std::vector<std::string>& keys,
                             const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data has one of the types given as template
     * arguments
     * \tparam Types: list of allowed types
     * \param[in] m: dictionary of keys and descriptions
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    template <typename... Types>
    ParametersValidator& add(const std::map<std::string, std::string>& m,
                             const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief check that the data associated with the given key is a strictly
     * positive integer
     * \param[in] k: key
     * \param[in] d: description
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    ParametersValidator& addStrictlyPositiveIntegerCheck(
        const std::string& k,
        const std::string& d,
        const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief check that the data associated with the given key is a strictly
     * positive integer
     * \param[in] k: key
     * \param[in] opts: option used to declare the key
     * \return the modified validator
     */
    ParametersValidator& addStrictlyPositiveIntegerCheck(
        const std::string& k, const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief declare a list of keys to be incompatible
     * \tparam Types: list of allowed types
     * \param[in] keys: keys
     * \param[in] opts: option used to declare if one of the given keys is
     * required.
     * \return the modified validator
     * \note given keys are automatically added to the authorized keys
     * \note none of those keys have to be declared required
     */
    template <typename... Types>
    ParametersValidator& addIncompatibleParametersList(
        const std::vector<std::string>& keys,
        const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief declare a list of keys to be incompatible
     * \tparam Types: list of allowed types
     * \param[in] m: keys and associated descriptions
     * \param[in] opts: option used to declare if one of the given keys is
     * required.
     * \return the modified validator
     * \note given keys are automatically added to the authorized keys
     * \note none of those keys have to be declared required
     */
    template <typename... Types>
    ParametersValidator& addIncompatibleParametersList(
        const std::map<std::string, std::string>& m,
        const AddArguments& opts = not_required) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>)));
    /*!
     * \brief declare a list of keys to be incompatible
     * \param[in] keys: keys
     * \param[in] opts: option used to declare if one of the given keys is
     * required.
     * \return the modified validator
     * \note given keys are automatically added to the authorized keys
     * \note none of those keys have to be declared required
     */
    ParametersValidator& addIncompatibleParametersList(
        const std::vector<std::string>& keys,
        const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief declare a list of keys to be incompatible
     * \param[in] m: keys and associated descriptions
     * \param[in] opts: option used to declare if one of the given keys is
     * required.
     * \return the modified validator
     * \note given keys are automatically added to the authorized keys
     * \note none of those keys have to be declared required
     */
    ParametersValidator& addIncompatibleParametersList(
        const std::map<std::string, std::string>& m,
        const AddArguments& opts = not_required) noexcept;
    /*!
     * \brief validate a dictionary of parameters
     * \param[in] throwing: dummy attribute to indicate that this function may
     * throw an exception
     * \param[in] parameters: tested parameters
     * \note this method shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available. If a
     * context is available, it is better to use the following pattern:
     *
     * \code{.cpp}
     * auto or_raise = ctx.getThrowingFailureHandler();
     * validator.validate(ctx, parameters) | or_raise;
     * \endcode
     */
    void validate(attributes::Throwing throwing,
                  const Parameters& parameters) const;
    /*!
     * \brief validate a dictionary of parameters
     * \param[in, out] ctx: execution context
     * \param[in] m: tested parameters
     * \return true on success
     */
    [[nodiscard]] bool validate(Context& ctx,
                                const Parameters& m) const noexcept;
    /*!
     * \brief get the list of allowed parameters and their associated
     * documentation
     * \return the list of allowed parameters and their associated
     * documentation
     */
    [[nodiscard]] const std::map<std::string, std::string, std::less<>>&
    getAllowedParameters() const noexcept;
    /*!
     * \brief get the description of a parameter
     * \return the description of a parameter, empty on failure
     * \param[in, out] ctx: execution context
     * \param[in] k: name of the parameter
     */
    std::optional<std::string> getDescription(
        Context& ctx, std::string_view k) const noexcept;

    //! \brief destructor
    ~ParametersValidator();

   private:
    /*!
     * \brief add a key without documentation.
     * \param[in] k: key
     * \param[in] opts: option used to declare the key
     * \note if the key already exists, its description is left unchanged.
     */
    void addKey(const std::string& k, const AddArguments& opts) noexcept;
    //! \brief list of incompatible keys
    std::vector<std::vector<std::string>> incompatibilities;
    //! \brief dictionary associating a valid key and its description
    std::map<std::string, std::string, std::less<>> allowed_keys;
    //! \brief list of required keys
    std::set<std::string> required_keys;
    //! \brief list of sets of keys for which at least one key is required
    std::vector<std::vector<std::string>> required_keys_in_set;
    //! \brief validators, sorted by keywords
    std::map<std::string, std::vector<ParameterValidator>> validators;
  };

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Utilities/ParametersValidator.ixx"

#endif /* LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_HXX */
