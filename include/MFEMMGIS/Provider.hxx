/*!
 * \file   MFEMMGIS/Provider.hxx
 * \brief  This file declares the `Provider` class
 * \date   12/12/2022
 */

#ifndef LIB_MFEM_MGIS_PROFILER_HXX
#define LIB_MFEM_MGIS_PROFILER_HXX

#include <string>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ExitStatus.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  // forward declarations
  struct QPDependency;
  struct QPEvaluatorsFactory;
  struct DependenciesManager;

  //! \brief this class describe an object able to solve a dependency
  struct MFEM_MGIS_EXPORT Provider {
    /*!
     * \brief return the name of the provider
     *
     * \note as most providers are also coupling items which has a `getName`
     * method, we named this method `getIdentifier` rather than `getName` to
     * avoid conflicts.
     */
    [[nodiscard]] virtual std::string getIdentifier() const noexcept = 0;
    /*!
     * \brief analyse the given dependency. If a provider can resolve this
     * dependency, it  shall declare itself as the provider of the dependency.
     * It can then declare additional dependencies.
     *
     * \note since the given dependency may not have concrete  specifications,
     * they must not be checked in this method. The implementation may check
     * that the specifications of the dependencies are compatible with what
     * the provider may provide using the
     * `DependencyBase::matchesSpecifications` method, but this is automatically
     * done by the  `DependenciesManager::setProvider` method .
     *
     * \note The `reportInvalidResolveDependencyCall` can be called to
     * report an error
     *
     * \param[in, out] ctx: execution context
     * \param[in] dm: dependencies manager
     * \param[in] d: dependency
     * \param[in] ts: time step stage
     */
    [[nodiscard]] virtual bool analyseDependency(
        Context &,
        DependenciesManager &,
        const QPDependency &,
        const TimeStepStage) const noexcept = 0;
    /*!
     * \brief resolve the dependency
     *
     * \note since the given dependency must have been resolved by this
     * provider, checking the specifications of the dependencies shall not
     * be required. The implementations of `resolveDependencies` are free to
     * call `Dependencies::checkSpecifications` for a somehow paranoid check.
     *
     * \note The `reportInvalidResolveDependencyCall` can be called to
     * report an error
     *
     * \param[in] ctx: execution context
     * \param[in] f: evaluator factory
     * \param[in] d: dependency
     * \param[in] ts: time step stage
     */
    [[nodiscard]] virtual bool resolveDependency(
        Context &,
        QPEvaluatorsFactory &,
        const QPDependency &,
        const TimeStepStage) const noexcept = 0;
    //! \brief destructor
    virtual ~Provider() noexcept;

   protected:
    /*!
     * \brief method that can be called to report an invalid call to the
     * resolveDependency method
     * \param[in] ctx: execution context
     * \param[in] d: dependency
     */
    [[nodiscard]] static bool reportInvalidResolveDependencyCall(
        Context &, const QPDependency &) noexcept;
  };  // end of Provider

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PROFILER_HXX */
