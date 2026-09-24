/*!
 * \file   MFEMMGIS/MaterialQuantityProviderSearch.hxx
 * \brief  This header declares functions searching, among the behaviour
 * integrators defined on a given location, the one which provides a given
 * gradient, thermodynamic force or internal state variable.
 * \date   23/09/2026
 */

#ifndef LIB_MFEMMGIS_MATERIALQUANTITYPROVIDERSEARCH_HXX
#define LIB_MFEMMGIS_MATERIALQUANTITYPROVIDERSEARCH_HXX

#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  // forward declarations
  struct AbstractBehaviourIntegrator;
  struct AbstractNonLinearEvolutionProblem;

  /*!
   * \brief result of the search of the behaviour integrator providing a
   * given quantity (gradient, thermodynamic force or internal state
   * variable) on a given location.
   *
   * Several behaviour integrators can be defined on the same location. Each
   * of them may or may not be associated with a material. A quantity is
   * generally defined by only one of those materials. The search is
   * successful if, and only if, exactly one behaviour integrator associated
   * with a material provides the requested quantity.
   */
  struct MaterialQuantityProviderSearchResult {
    //! \brief outcome of the search
    enum struct Status {
      //! \brief a unique provider has been found
      SUCCESS,
      //! \brief no behaviour integrator provides the requested quantity
      NO_PROVIDER,
      //! \brief several behaviour integrators provide the requested quantity
      MULTIPLE_PROVIDERS
    };
    //! \brief behaviour integrator providing the quantity, if unique
    OptionalReference<const AbstractBehaviourIntegrator> behaviour_integrator;
    //! \brief outcome of the search
    Status status = Status::SUCCESS;
  };  // end of struct MaterialQuantityProviderSearchResult

  /*!
   * \return if the search failed, i.e. if no unique provider was found
   * \param[in] r: result of the search
   */
  [[nodiscard]] inline bool isInvalid(
      const MaterialQuantityProviderSearchResult& r) noexcept {
    return isInvalid(r.behaviour_integrator);
  }  // end of isInvalid

  /*!
   * \brief search the behaviour integrator providing the given gradient on
   * the given location
   *
   * \return the result of the search. An empty value is returned if an error
   * occured, in which case the error is reported in the execution context.
   * Not finding a unique provider is not considered as an error: this
   * information is reported by the status of the result.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] l: location
   * \param[in] n: name of the gradient
   *
   * \note behaviour integrators are currently only defined on materials, so
   * no provider can be found on a boundary.
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<MaterialQuantityProviderSearchResult>
  hasGradientProvider(Context&,
                      const AbstractNonLinearEvolutionProblem&,
                      const LocationIdentifier&,
                      std::string_view) noexcept;
  /*!
   * \brief search the behaviour integrator providing the given
   * thermodynamic force on the given location
   *
   * \return the result of the search. An empty value is returned if an error
   * occured, in which case the error is reported in the execution context.
   * Not finding a unique provider is not considered as an error: this
   * information is reported by the status of the result.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] l: location
   * \param[in] n: name of the thermodynamic force
   *
   * \note behaviour integrators are currently only defined on materials, so
   * no provider can be found on a boundary.
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<MaterialQuantityProviderSearchResult>
  hasThermodynamicForceProvider(Context&,
                                const AbstractNonLinearEvolutionProblem&,
                                const LocationIdentifier&,
                                std::string_view) noexcept;
  /*!
   * \brief search the behaviour integrator providing the given internal
   * state variable on the given location
   *
   * \return the result of the search. An empty value is returned if an error
   * occured, in which case the error is reported in the execution context.
   * Not finding a unique provider is not considered as an error: this
   * information is reported by the status of the result.
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] l: location
   * \param[in] n: name of the internal state variable
   *
   * \note behaviour integrators are currently only defined on materials, so
   * no provider can be found on a boundary.
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<MaterialQuantityProviderSearchResult>
  hasInternalStateVariableProvider(Context&,
                                   const AbstractNonLinearEvolutionProblem&,
                                   const LocationIdentifier&,
                                   std::string_view) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_MATERIALQUANTITYPROVIDERSEARCH_HXX */
