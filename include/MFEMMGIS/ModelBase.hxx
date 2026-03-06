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
     * \param[in] m: mesh
     */
    ModelBase(const MeshDiscretization &) noexcept;
    //
    [[nodiscard]] std::string getIdentifier() const noexcept override final;
    MeshDiscretization getMeshDiscretization() const noexcept override;
    [[nodiscard]] VerbosityLevel getVerbosityLevel()
        const noexcept override final;
    void setVerbosityLevel(const VerbosityLevel) noexcept override final;
    void setLogStream(std::shared_ptr<std::ostream>) noexcept override final;
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer() noexcept
        override final;
    [[nodiscard]] std::vector<std::string> getLocations()
        const noexcept override;
    [[nodiscard]] std::optional<std::string> describe(
        Context &, const bool, const Parameters &) const noexcept override;
    [[nodiscard]] std::vector<std::string> getAvailablePostProcessings()
        const noexcept override;
    [[nodiscard]] bool addPostProcessing(Context &,
                                         std::string_view,
                                         const Parameters &) noexcept override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    //     [[nodiscard]] bool analyseDependency(Context &,
    //                                  DependenciesManager &,k
    //                                  const ValueDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool resolveDependency(Context &,
    //                                  ValueEvaluatorsFactory &,
    //                                  const ValueDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool analyseDependency(Context &,
    //                                  DependenciesManager &,
    //                                  const NodalDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool resolveDependency(Context &,
    //                                  NodalEvaluatorsFactory &,
    //                                  const NodalDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool analyseDependency(Context &,
    //                                  DependenciesManager &,
    //                                  const IPDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool resolveDependency(Context &,
    //                                  IPEvaluatorsFactory &,
    //                                  const IPDependency &,
    //                                  const TimeStepStage) const noexcept
    //                                  override;
    //     [[nodiscard]] bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept override;
    //     [[nodiscard]] bool initializeAfterResourcesAllocation(Context &)
    //     noexcept override;
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &, const real) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(Context &,
                                                  const TimeStep &,
                                                  const bool) noexcept override;
    std::optional<real> getNextTimeIncrement(
        Context &, const real, const real) const noexcept override;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
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
        std::function<bool(Context &, bool)>) noexcept;

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