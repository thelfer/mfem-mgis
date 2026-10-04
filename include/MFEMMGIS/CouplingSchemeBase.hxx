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
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    CouplingSchemeBase(Context &ctx, const MeshDiscretization &m) noexcept;

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] parameters: parameters
     */
    CouplingSchemeBase(Context &ctx,
                       const MeshDiscretization &m,
                       const Parameters &parameters);
    //
    MeshDiscretization getMeshDiscretization() const noexcept override;
    void setName(std::string_view n) noexcept override final;
    //! \return an empty list, a coupling scheme has no location
    [[nodiscard]] std::vector<std::string> getLocations()
        const noexcept override;
    //! \return the verbosity level of the scheme or the default one if unset
    [[nodiscard]] VerbosityLevel getVerbosityLevel()
        const noexcept override final;
    void setVerbosityLevel(const VerbosityLevel l) noexcept override final;
    void setLogStream(std::shared_ptr<std::ostream> s) noexcept override final;
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer() noexcept
        override final;
    /*!
     * \return the coupling items which are providers and, recursively, the
     * providers of the coupling items which are coupling schemes
     */
    [[nodiscard]] std::vector<const Provider *> getProviders() noexcept
        override;
    //    [[nodiscard]] bool add(Context &, const Parameters &) noexcept
    //    override;
    //     [[nodiscard]] bool addCouplingItem(Context &,
    //                                        std::string_view,
    //                                        std::string_view,
    //                                        const Parameters &) noexcept
    //                                        override;
    /*!
     * \brief add a new coupling item
     * \param[in, out] ctx: execution context
     * \param[in] i: coupling item, defined on the mesh of the scheme
     * \return true on success
     */
    [[nodiscard]] bool addCouplingItem(
        Context &ctx,
        std::shared_ptr<AbstractCouplingItem> i) noexcept override;
    //     [[nodiscard]] bool addModel(Context &,
    //                                 std::string_view,
    //                                 const Parameters &) noexcept override;
    [[nodiscard]] bool addModel(
        Context &ctx, std::shared_ptr<AbstractModel> m) noexcept override;
    [[nodiscard]] bool addModel(
        Context &ctx,
        std::shared_ptr<NonLinearEvolutionProblem> m) noexcept override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    //     [[nodiscard]] bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept override;
    //     [[nodiscard]] bool initializeAfterResourcesAllocation(
    //         Context &) noexcept override;
    /*!
     * \brief call `performInitializationTaksAtTheBeginningOfTheTimeStep` on
     * all coupling items
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return true on success
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \return the minimum of `te - t` and of the time increments proposed by
     * the coupling items
     * \param[in, out] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    std::optional<real> getNextTimeIncrement(
        Context &ctx, const real t, const real te) const noexcept override;
    /*!
     * \brief execute the initial post-processings of all coupling items
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept override;
    /*!
     * \brief execute the post-processings of all coupling items
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] b: boolean stating if the time at the end of the time
     * step is a post-processing time
     * \return true on success
     */
    [[nodiscard]] bool executePostProcessingTasks(
        Context &ctx, const TimeStep &ts, const bool b) noexcept override;
    /*!
     * \brief update all coupling items
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    /*!
     * \brief revert all coupling items
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~CouplingSchemeBase() noexcept override;

   protected:
    /*!
     * \brief a short structure saving the state of a Context
     * before its modification by the `update` function
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
     * \param[in, out] ctx: execution context
     * \param[in] m: coupling item
     * \return the state before the update
     */
    static ContextState update(Context &ctx, AbstractCouplingItem &m) noexcept;
    /*!
     * \brief restore the state of a context
     * \param[in, out] ctx: execution context
     * \param[in] s: context state
     */
    static void restore(Context &ctx, const ContextState &s) noexcept;
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
    //! \brief name of the coupling scheme, specified externally
    std::optional<std::string> name;
  };  // end of CouplingSchemeBase

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_COUPLING_SCHEME_BASE_HXX */
