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
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  /*!
   * This class is an abstract factory for quadrature point evaluators
   * evaluators, i.e. which provides methods to:
   *
   * 1. registrer generators of quadrature point evaluators stored by name,
   * material identifier and equivalent partial quadrature space.
   * 2. handle the generation of evaluators.
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorsFactory {
    /*!
     * \brief structure return by the `analyseDependencies` method
     */
    struct LocalDependenciesAnalysisOutput {
      //! \brief missing dependencies
      std::vector<QPEvaluatorDescription> missingIPDependencies;
    };
    //! \brief a simple alias
    using Generator = std::shared_ptr<AbstractQPEvaluatorGenerator>;
    /*!
     * \return a description of the given location usable in an error message
     * \param[in] m: mesh set on which the evaluator is defined.
     * \param[in] qid: quadrature id for which the evaluator is defined.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage) noexcept;
    //! \brief default constructor
    QPEvaluatorsFactory() noexcept;
    // disabling default constructors and assignement operators
    QPEvaluatorsFactory(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory(const QPEvaluatorsFactory &) = delete;
    QPEvaluatorsFactory &operator=(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory &operator=(const QPEvaluatorsFactory &) = delete;
    /*!
     * \return if a generator associated with the given partial quadrature
     * space, quadrature id and time set stage has been declared.
     *
     * \param[in] qspace: partial quadrature space
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     */
    [[nodiscard]] bool containsGenerator(
        const std::shared_ptr<const PartialQuadratureSpace>,
        const TimeStepStage,
        const std::string &) const noexcept;
    //     /*!
    //      * \brief register a new evaluator generator
    //      * \param[in] ctx: execution context.
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      * \param[in] g: generator
    //      */
    //     [[nodiscard]] bool registerGenerator(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const std::string &,
    //         const Generator &) noexcept;
    //     /*!
    //      * \brief generate an evaluator on the given mesh set and quadrature
    //      id.
    //      *
    //      * \param[in] ctx: execution context.
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] d: description of the evaluator
    //      *
    //      * \return a shared pointer to the evaluator. If the evaluation
    //      failed, the
    //      * shared pointer is empty.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const QPEvaluatorDescription &) const noexcept;
    //     /*!
    //      * \brief generate an evaluator on the given mesh set and quadrature
    //      id.
    //      *
    //      * \param[in] ctx: execution context.
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      *
    //      * \return a shared pointer to the evaluator. If the evaluation
    //      failed, the
    //      * shared pointer is empty.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const std::string &) const noexcept;
    //     /*!
    //      * \brief generate an evaluator on the given mesh set and quadrature
    //      id.
    //      *
    //      * \param[in] ctx: execution context.
    //      * \param[in] m: mesh set on which the evaluator is defined.
    //      * \param[in] s: stage in the time step
    //      * \param[in] d: description of the evaluator
    //      *
    //      * \return a shared pointer to the evaluator. If the evaluation
    //      failed, the
    //      * shared pointer is empty.
    //      *
    //      * \note the quadrature is unspecified. This call can only work if
    //      the
    //      * evaluator is avaiable for only one quadrature.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const size_type,
    //         const TimeStepStage,
    //         const QPEvaluatorDescription &) const noexcept;
    //     /*!
    //      * \brief generate an evaluator on the given mesh set and quadrature
    //      id.
    //      * param[in] ctx: execution context.
    //      * \param[in] m: mesh set on which the evaluator is defined.
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      *
    //      * \return a shared pointer to the evaluator. If the evaluation
    //      failed, the
    //      * shared pointer is empty.
    //      *
    //      * \note the quadrature is unspecified. This call can only work if
    //      the
    //      * evaluator is avaiable for only one quadrature.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const size_type,
    //         const TimeStepStage,
    //         const std::string &) const noexcept;
    //     /*!
    //      * \brief generate a set of evaluators on the given mesh set and
    //      quadrature
    //      * id.
    //      * \param[in] ctx: execution context.
    //      * \param[in] m: mesh set on which the evaluator is defined.
    //      * \param[in] s: stage in the time step
    //      * \param[in] descriptions: description of the evaluators
    //      *
    //      * \note the quadrature is unspecified. This call can only work if
    //      the
    //      * evaluator for each mesh set is avaiable for only one quadrature.
    //      */
    //     [[nodiscard]] std::optional<
    //         std::vector<std::pair<size_type,
    //         std::shared_ptr<AbstractQPEvaluator>>>>
    //     generate(Context &,
    //              const size_type,
    //              const TimeStepStage,
    //              const std::vector<QPEvaluatorDescription> &) const noexcept;
    //     /*!
    //      * \brief generate a set of evaluators on the given mesh set and
    //      quadrature
    //      * id.
    //      * \param[in] ctx: execution context.
    //      * \param[in] m: mesh set on which the evaluator is defined.
    //      * \param[in] s: stage in the time step
    //      * \param[in] names: names of the evaluators
    //      *
    //      * \note the quadrature is unspecified. This call can only work if
    //      the
    //      * evaluator for each mesh set is avaialble for only one quadrature.
    //      */
    //     [[nodiscard]] std::optional<
    //         std::vector<std::pair<size_type,
    //         std::shared_ptr<AbstractQPEvaluator>>>>
    //     generate(Context &,
    //              const size_type,
    //              const TimeStepStage,
    //              const std::vector<std::string> &) const noexcept;
    //! \return the list of registered generators
    [[nodiscard]] std::string getRegisteredGeneratorsList() const noexcept;
    //! \brief destructor
    ~QPEvaluatorsFactory() noexcept;

   private:
    /*!
     * \brief a simple alias
     */
    using GeneratorsContainer =
        std::map<LocationIdentifier,
                 std::map<std::shared_ptr<const PartialQuadratureSpace>,
                          std::map<std::string, Generator>>>;
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
    //     /*!
    //      * \brief check if the dependencies of the given evaluator are met
    //      * \param[in] m: registered generators
    //      * \param[in] n: name of the evaluator
    //      * \param[in] previous_dependencies: list of the previous iterations
    //      */
    //     LocalDependenciesAnalysisOutput analyseDependencies(
    //         Context &,
    //         const std::map<std::string, Generator> &,
    //         const std::string &n,
    //         const std::vector<QPEvaluatorDescription> &) const noexcept;
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
