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
        const PartialQuadratureSpaceIdentifiersManager &) noexcept;
    /*!
     * \brief enumeration of filters usable as argument of the
     * analyseDependencies method
     */
    enum AnalyseDependenciesFilter {
      ONLY_REQUIRED,
      ONLY_OPTIONAL,
      ALL
    };
    /*!
     * \brief structure return by the `analyseDependencies` method
     */
    struct DependenciesAnalysisOutput {
      //! \brief missing dependencies at the beginning of the time step
      std::vector<QPDependency> missingQPDependencies_bts;
      //! \brief missing dependencies at the beginning of the end of the time
      //! step
      std::vector<QPDependency> missingQPDependencies_ets;
    };
    /*!
     * \return a description of the location of the given dependency at
     * integration point usable in an error message
     *
     * \param[in] d: dependency.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(const QPDependency &,
                                              const TimeStepStage) noexcept;
    /*!
     * \brief add a new dependency at integration points
     * param[in] ctx: execution context
     * \param[in] s: time step stage
     * \param[in] d: dependency description
     */
    [[nodiscard]] bool declareDependency(Context &,
                                         const TimeStepStage,
                                         const QPDependency &) noexcept;
    /*!
     * \brief set the provider of the dependency at integration points for the
     * given location with the given name
     *
     * \param[in] ctx: execution context
     * \param[in] pr: provider
     * \param[in] d: dependency
     * \param[in] qid: quadrature id
     * \param[in] nc: number of components
     * \param[in] ts: time step stage
     *
     * \note the quadrature id must be passed to properly treat the case when
     * the dependency does not specify it.
     */
    [[nodiscard]] bool setProvider(
        Context &,
        const Provider &,
        const QPDependency &,
        std::shared_ptr<const PartialQuadratureSpace>,
        const size_type,
        const TimeStepStage) noexcept;
    /*!
     * \brief analyse dependencies
     * \param[in] f: filter
     */
    [[nodiscard]] DependenciesAnalysisOutput analyseDependencies(
        const AnalyseDependenciesFilter =
            AnalyseDependenciesFilter::ONLY_REQUIRED) const noexcept;

   private:
    /*!
     * \return a local manager for dependencies at integration points for the
     * given time step stage
     *
     * \param[in] m: material identifier
     * \param[in] s: time step stage
     */
    std::vector<QPDependency> &getLocalQPDependenciesManager(
        const size_type, const TimeStepStage) noexcept;
    //! \brief partial quadrature space identifiers
    const PartialQuadratureSpaceIdentifiersManager &qids;
    //! \brief list of registred dependencies at integration points
    std::array<std::map<size_type, std::vector<QPDependency>>, 2u>
        registeredQPDependencies;
  };  // end of DependenciesManager

  /*!
   * \return a pair containing the number of dependencies and a description of
   * thoses dependencies in the form of a list. \param[in] a: output of the
   * dependencies analysis
   */
  [[nodiscard]] std::pair<size_type, std::string> getDescription(
      const DependenciesManager::DependenciesAnalysisOutput &) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_DEPENDENCIESMANAGER_HXX */