/*!
 * \file   MFEMMGIS/CouplingSchemeBase.hxx
 * \brief  This file declares the `CouplingSchemeBase` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_COUPLING_SCHEME_BASE_HXX
#define LIB_MFEM_MGIS_COUPLING_SCHEME_BASE_HXX

#include <map>
#include <vector>
#include <string>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractCouplingScheme.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;

  //! \brief a base class for most coupling schemes
  struct MFEM_MGIS_EXPORT CouplingSchemeBase : AbstractCouplingScheme {
    //! \return a description of the parameters of this scheme
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    CouplingSchemeBase(const MeshDiscretization &) noexcept;
    //
    MeshDiscretization getMeshDiscretization() const noexcept override;
    [[nodiscard]] std::vector<std::string> getLocations()
        const noexcept override;
    [[nodiscard]] VerbosityLevel getVerbosityLevel()
        const noexcept override final;
    void setVerbosityLevel(const VerbosityLevel) noexcept override final;
    void setLogStream(std::shared_ptr<std::ostream>) noexcept override final;
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer() noexcept
        override final;
    [[nodiscard]] std::vector<const Provider *> getProviders() noexcept
        override;
    //    [[nodiscard]] bool add(Context &, const Parameters &) noexcept
    //    override;
    //     [[nodiscard]] bool addCouplingItem(Context &,
    //                                        std::string_view,
    //                                        std::string_view,
    //                                        const Parameters &) noexcept
    //                                        override;
    [[nodiscard]] bool addCouplingItem(
        Context &, std::shared_ptr<AbstractCouplingItem>) noexcept override;
    //     [[nodiscard]] bool addModel(Context &,
    //                                 std::string_view,
    //                                 const Parameters &) noexcept override;
    [[nodiscard]] bool addModel(
        Context &, std::shared_ptr<AbstractModel>) noexcept override;
    [[nodiscard]] bool addModel(
        Context &,
        std::shared_ptr<NonLinearEvolutionProblem>) noexcept override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    //     [[nodiscard]] bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept override;
    //     [[nodiscard]] bool initializeAfterResourcesAllocation(
    //         Context &) noexcept override;
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept override;
    std::optional<real> getNextTimeIncrement(
        Context &, const real, const real) const noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &, const real) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(Context &,
                                                  const TimeStep &,
                                                  const bool) noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
    //! \brief destructor
    ~CouplingSchemeBase() noexcept override;

   protected:
    /*!
     * \brief a short structure saving the state a Context
     * before its modification by the `updateState` function
     */
    struct [[nodiscard]] ContextState {
      //! \brief verbosity level
      VerbosityLevel verbosity_level;
      //! \brief log stream
      std::shared_ptr<std::ostream> log_stream;
    };  // end of ContextState
    /*!
     * \brief function updating a context to take into account
     * the settings of a coupling item.
     * \param[in,out] ctx: execution context
     * \param[in] i: coupling item
     * \return the state before the update
     */
    static ContextState update(Context &, AbstractCouplingItem &) noexcept;
    /*!
     * \brief restore the state of a context
     * \param[in,out] ctx: execution context
     * \param[in] s: context state
     */
    static void restore(Context &, const ContextState &) noexcept;
    //! \return a description of the coupling items
    [[nodiscard]] virtual std::string getCouplingItemsDescription()
        const noexcept;
    //! \brief mesh discretization
    MeshDiscretization mesh;
    //! \brief list of registered coupling items
    std::vector<std::shared_ptr<AbstractCouplingItem>> items;
    //! \brief the verbosity level associated with the coupling scheme
    std::optional<VerbosityLevel> verbosity_level;
    //! \brief a log stream associated with the coupling scheme
    std::shared_ptr<std::ostream> log_stream;
  };  // end of CouplingSchemeBase

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_COUPLING_SCHEME_BASE_HXX */
