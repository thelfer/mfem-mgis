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
     * \param[in,out] ctx: execution context
     * \param[in] m: mesh
     */
    ModelBase(Context &ctx, const MeshDiscretization &m) noexcept;

    /*!
     * \brief constructor
     * \param[in,out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] parameters: parameters
     */
    ModelBase(Context &ctx,
              const MeshDiscretization &m,
              const Parameters &parameters);
    //
    [[nodiscard]] std::string getIdentifier() const noexcept override final;
    MeshDiscretization getMeshDiscretization() const noexcept override;
    [[nodiscard]] VerbosityLevel getVerbosityLevel()
        const noexcept override final;
    void setName(std::string_view n) noexcept override final;
    void setVerbosityLevel(const VerbosityLevel l) noexcept override final;
    void setLogStream(std::shared_ptr<std::ostream> s) noexcept override final;
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer() noexcept
        override final;
    [[nodiscard]] std::vector<std::string> getLocations()
        const noexcept override;
    [[nodiscard]] std::optional<std::string> describe(
        Context &ctx,
        const bool b,
        const Parameters &parameters) const noexcept override;
    [[nodiscard]] std::vector<std::string> getAvailablePostProcessings()
        const noexcept override;
    [[nodiscard]] bool addPostProcessing(
        Context &ctx,
        std::string_view n,
        const Parameters &params) noexcept override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    [[nodiscard]] bool analyseDependency(
        Context &ctx,
        DependenciesManager &dm,
        const QPDependency &d,
        const TimeStepStage ts) const noexcept override;
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
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(
        Context &ctx, const TimeStep &ts, const bool b) noexcept override;
    std::optional<real> getNextTimeIncrement(
        Context &ctx, const real t, const real te) const noexcept override;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~ModelBase() noexcept override;

   protected:
    // \brief return a detailed description of the model
    [[nodiscard]] virtual std::string getDetailedDescription() const noexcept;
    // \brief return a description of the unknown fields
    [[nodiscard]] virtual std::string getUnknownFieldsDescription()
        const noexcept;
    // \brief return a description of the state variables
    [[nodiscard]] virtual std::string getStateVariablesDescription()
        const noexcept;
    // \brief return a description of the dependencies
    [[nodiscard]] virtual std::string getDependenciesDescription()
        const noexcept;
    //! \brief add a post-processing (see executePostProceccing for details)
    virtual void addPostProcessing(
        std::function<bool(Context &, bool)> p) noexcept;

    //! \brief name of the model, specified externally
    std::optional<std::string> name;

   private:
    //! \brief mesh discretization
    MeshDiscretization mesh;
    //! \brief list of registred post-processings
    std::vector<std::function<bool(Context &, bool)>> postProcessings;
    //! \brief the verbosity level associated with the model
    std::optional<VerbosityLevel> verbosityLevel;
    //! \brief a log stream associated with the model
    std::shared_ptr<std::ostream> logStream;
  };  // end of class ModelBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MODEL_BASE_HXX */
