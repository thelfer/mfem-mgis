/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluatorsFactory.hxx
 * \brief  This file declares the QPEvaluatorsFactory class
 * \author Thomas Helfer
 * \date   20/10/2022
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORFACTORY_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORFACTORY_HXX

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  /*!
   * This class is an abstract factory for quadrature point evaluators
   * evaluators, i.e. which provides methods to:
   *
   * 1. registrer generators of quadrature point evaluators stored by name,
   * material identifier and equivalent partial quadrature space.
   * 2. handle the generation of evaluators.
   */
  struct MFEM_MGIS_EVALUATOR QPEvaluatorsFactory {
    /*!
     * \brief structure return by the `analyseDependencies` method
     */
    struct LocalDependenciesAnalysisOutput : AdvancedExitStatus {
      //! inherinting constructors
      using AdvancedExitStatus::AdvancedExitStatus;
      //! \brief missing dependencies
      std::vector<QPEvaluatorDescription> missingIPDependencies;
    };
    //! \brief a simple alias
    using Generator = std::shared_ptr<AbstractQPEvaluatorGenerator>;
    /*!
     * \return a description of the given location usable in an error message
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(const size_type,
                                              const TimeStepStage) noexcept;
    /*!
     * \return a description of the given location usable in an error message
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage) noexcept;
    /*!
     * \brief default constructor
     * \param[in] m: model
     */
    QPEvaluatorsFactory(ResourcesManager &) noexcept;
    //! \brief return the underlying resources manager
    ResourcesManager &getResourcesManager() noexcept;
    //! \brief return the underlying resources manager
    const ResourcesManager &getResourcesManager() const noexcept;
    /*!
     * \return if a generator associated with the given mesh set, quadrature id
     * and time set stage has been declared. \param[in] m: mesh set on which the
     * evaluator is defined. \param[in] qid: quadrature id for which the
     * evaluator is defined. \param[in] s: stage in the time step \param[in] n:
     * name of the evaluator
     */
    [[nodiscard]] bool containsGenerator(
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &) const noexcept;
    /*!
     * \brief register a new evaluator generator
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] g: generator
     */
    [[nodiscard]] bool registerGenerator(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const Generator &) noexcept;
    /*!
     * \brief register a new generator for an uniform scalar evaluator
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] v: value
     */
    [[nodiscard]] bool registerGenerator(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const Real) noexcept;
    /*!
     * \brief register a new generator for an uniform evaluator
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] v: value
     */
    [[nodiscard]] bool registerGenerator(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const ConstMatrixRef<> &) noexcept;
    /*!
     * \brief register a new evaluator from a field on integration points
     * param[in] ctx: execution context.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] f: field
     */
    template <size_type nRows, size_type nCols>
    [[nodiscard]] bool registerGenerator(
        Context &,
        const TimeStepStage,
        const std::string &,
        const IPField<nRows, nCols> &) noexcept;
    /*!
     * \brief register a new evaluator from a view on a field on integration
     * points param[in] ctx: execution context. \param[in] s: stage in the time
     * step \param[in] n: name of the evaluator \param[in] f: field \param[in]
     * bs: block specifications defining the view
     */
    template <size_type nRows, size_type nCols>
    [[nodiscard]] bool registerGenerator(
        Context &,
        const TimeStepStage,
        const std::string &,
        const IPField<nRows, nCols> &,
        const IPFieldViewStorageSpecifications &) noexcept;
    /*!
     * \brief register a new evaluator from a field on integration points
     * param[in] ctx: execution context.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] f: field
     */
    [[nodiscard]] bool registerGenerator(Context &,
                                         const TimeStepStage,
                                         const std::string &,
                                         const ConstIPFieldView &) noexcept;
    /*!
     * \brief register a new evaluator for a function of an evaluator
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] f: function
     * \param[in] a: name of the evaluator
     */
    [[nodiscard]] bool registerGenerator(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const std::function<Real(const Real)> &,
        const std::string &) noexcept;
    /*!
     * \brief register a new evaluator for a function of two evaluators
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     * \param[in] f: function
     * \param[in] a1: name of the first argument
     * \param[in] a2: name of the second argument
     */
    [[nodiscard]] bool registerGenerator(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const std::function<Real(const Real, const Real)> &,
        const std::string &,
        const std::string &) noexcept;
    /*!
     * \brief check if the dependencies of the given evaluator are met
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] d: description of the dependency
     *
     * The result is to be understood as follows:
     *
     * - one or more missing dependencies is considered as a recoverable error
     * - if inconsistent dependencies are detected, an unrecoverable error is
     * returned
     * - if a cyclic dependency is detected, an unrecoverable error is returned
     *
     * If an unrecoverable error is returned, no missing dependencies are
     * reported by convention.
     */
    LocalDependenciesAnalysisOutput analyseDependencies(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const QPEvaluatorDescription &) const noexcept;
    /*!
     * \brief check if the dependencies of the given set of evaluators are met
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] deps: dependencies are met
     *
     * The result is to be understood as follows:
     *
     * - one or more missing dependencies is considered as a recoverable error
     * - if inconsistent dependencies are detected, an unrecoverable error is
     * returned
     * - if a cyclic dependency is detected, an unrecoverable error is returned
     *
     * If an unrecoverable error is returned, no missing dependencies are
     * reported by convention.
     */
    LocalDependenciesAnalysisOutput analyseDependencies(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::vector<QPEvaluatorDescription> &) const noexcept;
    /*!
     * \brief check if the dependencies of the given evaluator are met
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     *
     * The result is to be understood as follows:
     *
     * - one or more missing dependencies is considered as a recoverable error
     * - if no generator named 'n' is found, an unrecoverable error is returned
     * - if inconsistent dependencies are detected, an unrecoverable error is
     * returned
     * - if a cyclic dependency is detected, an unrecoverable error is returned
     *
     * If an unrecoverable error is returned, no missing dependencies are
     * reported by convention.
     */
    LocalDependenciesAnalysisOutput analyseDependencies(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &) const noexcept;
    /*!
     * \brief check if the dependencies of the given set of evaluators are met
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] names: names of the evaluators
     *
     * The result is to be understood as follows:
     *
     * - one or more missing dependencies is considered as a recoverable error
     * - if one of the given names is not associated with a generator, an
     * unrecoverable error is returned
     * - if inconsistent dependencies are detected, an unrecoverable error is
     * returned
     * - if a cyclic dependency is detected, an unrecoverable error is returned
     *
     * If an unrecoverable error is returned, no missing dependencies are
     * reported by convention.
     */
    LocalDependenciesAnalysisOutput analyseDependencies(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::vector<std::string> &) const noexcept;
    /*!
     * \brief generate an evaluator on the given mesh set and quadrature id.
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] d: description of the evaluator
     *
     * \return a shared pointer to the evaluator. If the evaluation failed, the
     * shared pointer is empty.
     */
    [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const QPEvaluatorDescription &) const noexcept;
    /*!
     * \brief generate an evaluator on the given mesh set and quadrature id.
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     *
     * \return a shared pointer to the evaluator. If the evaluation failed, the
     * shared pointer is empty.
     */
    [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &) const noexcept;
    /*!
     * \brief generate a set of evaluators on the given mesh set and quadrature
     * id. param[in] ctx: execution context. \param[in] m: mesh set on which the
     * evaluator is defined. \param[in] qid: quadrature id for which the
     * evaluator is defined. \param[in] s: stage in the time step \param[in]
     * descriptions: description of the evaluators
     */
    [[nodiscard]] std::optional<
        std::vector<std::pair<size_type, std::shared_ptr<AbstractQPEvaluator>>>>
    generate(Context &,
             const size_type,
             const std::shared_ptr<const PartialQuadratureSpace>,
             const TimeStepStage,
             const std::vector<QPEvaluatorDescription> &) const noexcept;
    /*!
     * \brief generate a set of evaluators on the given mesh set and quadrature
     * id. param[in] ctx: execution context. \param[in] m: mesh set on which the
     * evaluator is defined. \param[in] qid: quadrature id for which the
     * evaluator is defined. \param[in] s: stage in the time step \param[in]
     * names: names of the evaluators
     */
    [[nodiscard]] std::optional<
        std::vector<std::pair<size_type, std::shared_ptr<AbstractQPEvaluator>>>>
    generate(Context &,
             const size_type,
             const std::shared_ptr<const PartialQuadratureSpace>,
             const TimeStepStage,
             const std::vector<std::string> &) const noexcept;
    /*!
     * \brief generate an evaluator on the given mesh set and quadrature id.
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] d: description of the evaluator
     *
     * \return a shared pointer to the evaluator. If the evaluation failed, the
     * shared pointer is empty.
     *
     * \note the quadrature is unspecified. This call can only work if the
     * evaluator is avaiable for only one quadrature.
     */
    [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
        Context &,
        const size_type,
        const TimeStepStage,
        const QPEvaluatorDescription &) const noexcept;
    /*!
     * \brief generate an evaluator on the given mesh set and quadrature id.
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     *
     * \return a shared pointer to the evaluator. If the evaluation failed, the
     * shared pointer is empty.
     *
     * \note the quadrature is unspecified. This call can only work if the
     * evaluator is avaiable for only one quadrature.
     */
    [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
        Context &,
        const size_type,
        const TimeStepStage,
        const std::string &) const noexcept;
    /*!
     * \brief generate a set of evaluators on the given mesh set and quadrature
     * id. param[in] ctx: execution context. \param[in] m: mesh set on which the
     * evaluator is defined. \param[in] s: stage in the time step \param[in]
     * descriptions: description of the evaluators
     *
     * \note the quadrature is unspecified. This call can only work if the
     * evaluator for each mesh set is avaiable for only one quadrature.
     */
    [[nodiscard]] std::optional<
        std::vector<std::pair<size_type, std::shared_ptr<AbstractQPEvaluator>>>>
    generate(Context &,
             const size_type,
             const TimeStepStage,
             const std::vector<QPEvaluatorDescription> &) const noexcept;
    /*!
     * \brief generate a set of evaluators on the given mesh set and quadrature
     * id. param[in] ctx: execution context. \param[in] m: mesh set on which the
     * evaluator is defined. \param[in] s: stage in the time step \param[in]
     * names: names of the evaluators
     *
     * \note the quadrature is unspecified. This call can only work if the
     * evaluator for each mesh set is avaiable for only one quadrature.
     */
    [[nodiscard]] std::optional<
        std::vector<std::pair<size_type, std::shared_ptr<AbstractQPEvaluator>>>>
    generate(Context &,
             const size_type,
             const TimeStepStage,
             const std::vector<std::string> &) const noexcept;
    //! \return the list of registered generators
    [[nodiscard]] std::string getRegisteredGeneratorsList() const noexcept;
    //! \brief destructor
    ~QPEvaluatorsFactory() noexcept;

   private:
    // disabling default constructors and assignement operators
    QPEvaluatorsFactory() = delete;
    QPEvaluatorsFactory(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory(const QPEvaluatorsFactory &) = delete;
    QPEvaluatorsFactory &operator=(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory &operator=(const QPEvaluatorsFactory &) = delete;
    //! \brief a simple alias
    using GeneratorsContainer =
        std::map<const MeshSet *,
                 std::map<std::shared_ptr<const PartialQuadratureSpace>,
                          std::map<std::string, Generator>>>;
    /*!
     * \brief perform basics checks when defining a view on an ip field
     * \param[in] ectx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     * \param[in] n: name of the generator to be registered
     * \param[in] fn: name of the source field
     * \param[in] r: numbers of rows of the source field
     * \param[in] f: numbers of columns of the source field
     * \param[in] bs: block specifications defining the view
     */
    static [[nodiscard]] bool checkIPFieldViewSpecifications(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &,
        const std::string &,
        const size_type,
        const size_type,
        const IPFieldViewStorageSpecifications &) noexcept;
    /*!
     * \brief check if:
     * - the given mesh set is defined on the same model
     * - the given mesh contains only elements of one geometric type
     *
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     */
    [[nodiscard]] bool checkMeshSetConsistency(Context &,
                                               const size_type) const noexcept;
    /*!
     * \brief check if:
     * - the given mesh set is defined on the same model
     * - the given mesh contains only elements of one geometric type
     * - the given mesh has a geometric support compatible with the given
     * quadrature
     *
     * param[in] ctx: execution context.
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     */
    [[nodiscard]] bool checkMeshSetConsistency(
        Context &,
        const size_type,
        const std::shared_ptr<const PartialQuadratureSpace>) const noexcept;
    /*!
     * \return generators associated with the given stage in the time step
     * \param[in] s: stage in the time step
     */
    GeneratorsContainer &getGeneratorsContainer(const TimeStepStage) noexcept;
    /*!
     * \return fields associated with the given stage in the time step
     * \param[in] s: stage in the time step
     */
    const GeneratorsContainer &getGeneratorsContainer(
        const TimeStepStage) const noexcept;
    /*!
     * \return the list of registered generators for the given time step stage
     * \param[in] s: stage in the time step
     */
    [[nodiscard]] std::string getRegisteredGeneratorsList(
        const TimeStepStage) const noexcept;
    /*!
     * \brief check if the dependencies of the given evaluator are met
     * \param[in] m: registered generators
     * \param[in] n: name of the evaluator
     * \param[in] previous_dependencies: list of the previous iterations
     */
    LocalDependenciesAnalysisOutput analyseDependencies(
        Context &,
        const std::map<std::string, Generator> &,
        const std::string &n,
        const std::vector<QPEvaluatorDescription> &) const noexcept;
    /*!
     * \brief list of registered generators computing values of the beginning of
     * the time step, sorted by mesh sets, quadrature id's and name
     */
    GeneratorsContainer generatorsAtTheBeginningOfTheTimeStep;
    /*!
     * \brief list of registered generators computing values of the end of the
     * time step, sorted by mesh sets, quadrature id's and name
     */
    GeneratorsContainer generatorsAtTheEndOfTheTimeStep;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORFACTORY_HXX */
