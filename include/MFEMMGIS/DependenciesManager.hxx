/*!
 * \file   MFEMMGIS/DependenciesManager.hxx
 * \brief  This file declares the `DependenciesManager` class
 * \author Thomas Helfer
 * \date   02/04/2026
 */

#ifndef LIB_MFEMMGIS_DEPENDENCIESMANAGER_HXX
#define LIB_MFEMMGIS_DEPENDENCIESMANAGER_HXX

#include <map>
#include <array>
#include <string>
#include <utility>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  struct QPDependency;
  struct Provider;
  struct QPEvaluatorsFactory;
  struct PartialQuadratureSpaceIdentifiersManager;

  /*!
   * \brief a class used to collect dependencies
   */
  struct MFEM_MGIS_EXPORT DependenciesManager {
    /*!
     * \brief constructor
     * \param[in] m: partial quadrature space identifiers
     */
    DependenciesManager(
        const PartialQuadratureSpaceIdentifiersManager &m) noexcept;
    /*!
     * \brief enumeration of filters usable as argument of the
     * analyseDependencies method: only required dependencies, only optional
     * dependencies or all dependencies
     */
    enum AnalyseDependenciesFilter {
      ONLY_REQUIRED,  //!< only required dependencies
      ONLY_OPTIONAL,  //!< only optional dependencies
      ALL             //!< all dependencies
    };
    /*!
     * \brief structure returned by the `analyseDependencies` method
     */
    struct DependenciesAnalysisOutput {
      //! \brief missing dependencies at the beginning of the time step
      std::vector<QPDependency> missingQPDependencies_bts;
      //! \brief missing dependencies at the end of the time step
      std::vector<QPDependency> missingQPDependencies_ets;
    };
    /*!
     * \return a description of the location of the given dependency at
     * integration point usable in an error message
     *
     * \param[in] d: dependency.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(const QPDependency &d,
                                              const TimeStepStage s) noexcept;
    /*!
     * \brief add a new dependency at integration points
     * \param[in, out] ctx: execution context
     * \param[in] s: time step stage
     * \param[in] d: dependency description
     * \return true on success
     */
    [[nodiscard]] bool declareDependency(Context &ctx,
                                         const TimeStepStage s,
                                         const QPDependency &d) noexcept;
    /*!
     * \brief set the provider of the dependency at integration points for the
     * given location with the given name
     *
     * \param[in, out] ctx: execution context
     * \param[in] pr: provider
     * \param[in] d: dependency
     * \param[in] qspace: partial quadrature space
     * \param[in] nc: number of components
     * \param[in] ts: time step stage
     * \return true on success
     *
     * \note the partial quadrature space must be passed to properly treat the
     * case when the dependency does not specify it.
     */
    [[nodiscard]] bool setProvider(
        Context &ctx,
        const Provider &pr,
        const QPDependency &d,
        std::shared_ptr<const PartialQuadratureSpace> qspace,
        const size_type nc,
        const TimeStepStage ts) noexcept;
    /*!
     * \brief analyse dependencies
     * \param[in] f: filter
     * \return the dependencies without provider selected by the filter
     */
    [[nodiscard]] DependenciesAnalysisOutput analyseDependencies(
        const AnalyseDependenciesFilter f =
            AnalyseDependenciesFilter::ONLY_REQUIRED) const noexcept;
    /*!
     * \brief resolve all dependencies
     * \param[in, out] ctx: execution context
     * \param[in, out] f: factory
     * \return true on success
     */
    [[nodiscard]] bool resolveDependencies(
        Context &ctx, QPEvaluatorsFactory &f) const noexcept;

   private:
    /*!
     * \return a local manager for dependencies at integration points for the
     * given time step stage
     *
     * \param[in] l: location
     * \param[in] s: time step stage
     */
    std::vector<QPDependency> &getLocalQPDependenciesManager(
        const LocationIdentifier l, const TimeStepStage s) noexcept;
    //! \brief partial quadrature space identifiers
    const PartialQuadratureSpaceIdentifiersManager &qids;
    //! \brief list of registered dependencies at integration points
    std::array<std::map<LocationIdentifier, std::vector<QPDependency>>, 2u>
        registeredQPDependencies;
  };  // end of DependenciesManager

  /*!
   * \return a pair containing the number of dependencies and a description of
   * those dependencies in the form of a list.
   * \param[in] a: output of the dependencies analysis
   */
  [[nodiscard]] std::pair<size_type, std::string> getDescription(
      const DependenciesManager::DependenciesAnalysisOutput &a) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_DEPENDENCIESMANAGER_HXX */
