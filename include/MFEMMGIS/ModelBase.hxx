/*!
 * \file   MFEMMGIS/ModelBase.hxx
 * \brief  This file declares the  `ModelBase` class
 * \date   15/11/2022
 */

#ifndef LIB_MFEM_MGIS_MODEL_BASE_HXX
#define LIB_MFEM_MGIS_MODEL_BASE_HXX

#include <map>
#include <vector>
#include <string>
#include <optional>
#include <functional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractModel.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;

  //! \brief a base class for most model
  struct MFEM_MGIS_EXPORT ModelBase : AbstractModel {
    //! \return a description of the parameters of this model
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    ModelBase(Context &ctx, const MeshDiscretization &m) noexcept;

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] parameters: parameters
     */
    ModelBase(Context &ctx,
              const MeshDiscretization &m,
              const Parameters &parameters);
    //
    [[nodiscard]] std::string getIdentifier() const noexcept override final;
    MeshDiscretization getMeshDiscretization() const noexcept override;
    //! \return the verbosity level of the model or the default one if unset
    [[nodiscard]] VerbosityLevel getVerbosityLevel()
        const noexcept override final;
    void setName(std::string_view n) noexcept override final;
    void setVerbosityLevel(const VerbosityLevel l) noexcept override final;
    void setLogStream(std::shared_ptr<std::ostream> s) noexcept override final;
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer() noexcept
        override final;
    //! \return an empty list of locations
    [[nodiscard]] std::vector<std::string> getLocations()
        const noexcept override;
    /*!
     * \return a description of the model
     * \param[in, out] ctx: execution context
     * \param[in] b: boolean being the default value for information requests
     * \param[in] parameters: information requests. Supported requests are
     * `HelpOptions`, `ShortDescription`, `DetailedDescription`,
     * `UnknownFields`, `StateVariables` and `Dependencies`.
     */
    [[nodiscard]] std::optional<std::string> describe(
        Context &ctx,
        const bool b,
        const Parameters &parameters) const noexcept override;
    //! \return an empty list
    [[nodiscard]] std::vector<std::string> getAvailablePostProcessings()
        const noexcept override;
    /*!
     * \brief report that no post-processing named `n` is available
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the post-processing
     * \param[in] params: parameters defining the post-processing
     * \return false
     */
    [[nodiscard]] bool addPostProcessing(
        Context &ctx,
        std::string_view n,
        const Parameters &params) noexcept override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    /*!
     * \brief do nothing: by default, a model provides no dependency
     * \param[in, out] ctx: execution context
     * \param[in, out] dm: dependencies manager
     * \param[in] d: dependency
     * \param[in] ts: time step stage
     * \return true
     */
    [[nodiscard]] bool analyseDependency(
        Context &ctx,
        DependenciesManager &dm,
        const QPDependency &d,
        const TimeStepStage ts) const noexcept override;
    /*!
     * \brief report an invalid call: by default, a model provides no
     * dependency
     * \param[in, out] ctx: execution context
     * \param[in, out] f: evaluators factory
     * \param[in] d: dependency
     * \param[in] ts: time step stage
     * \return false
     */
    [[nodiscard]] bool resolveDependency(
        Context &ctx,
        QPEvaluatorsFactory &f,
        const QPDependency &d,
        const TimeStepStage ts) const noexcept override;
    //     [[nodiscard]] bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept override;
    //     [[nodiscard]] bool initializeAfterResourcesAllocation(Context &)
    //     noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return true
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time
     * \return true
     */
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept override;
    /*!
     * \brief execute the registered post-processings
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] b: boolean stating if the time at the end of the time
     * step is a post-processing time
     * \return true on success
     */
    [[nodiscard]] bool executePostProcessingTasks(
        Context &ctx, const TimeStep &ts, const bool b) noexcept override;
    /*!
     * \return the time remaining until the end of the temporal sequence,
     * minimum over all MPI processes
     * \param[in, out] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    std::optional<real> getNextTimeIncrement(
        Context &ctx, const real t, const real te) const noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return `ExitStatus::success` and an empty output
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \return true
     */
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \return true
     */
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~ModelBase() noexcept override;

   protected:
    //! \return a detailed description of the model
    [[nodiscard]] virtual std::string getDetailedDescription() const noexcept;
    //! \return a description of the unknown fields
    [[nodiscard]] virtual std::string getUnknownFieldsDescription()
        const noexcept;
    //! \return a description of the state variables
    [[nodiscard]] virtual std::string getStateVariablesDescription()
        const noexcept;
    //! \return a description of the dependencies
    [[nodiscard]] virtual std::string getDependenciesDescription()
        const noexcept;
    /*!
     * \brief add a post-processing, see `executePostProcessingTasks`
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(
        std::function<bool(Context &, bool)> p) noexcept;

    //! \brief name of the model, specified externally
    std::optional<std::string> name;

   private:
    //! \brief mesh discretization
    MeshDiscretization mesh;
    //! \brief list of registered post-processings
    std::vector<std::function<bool(Context &, bool)>> postProcessings;
    //! \brief the verbosity level associated with the model
    std::optional<VerbosityLevel> verbosityLevel;
    //! \brief a log stream associated with the model
    std::shared_ptr<std::ostream> logStream;
  };  // end of class ModelBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MODEL_BASE_HXX */
