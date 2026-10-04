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
   * \brief abstract factory for quadrature point evaluators
   *
   * This class is an abstract factory for quadrature point evaluators,
   * i.e. which provides methods to:
   *
   * 1. register generators of quadrature point evaluators stored by name,
   * material identifier and equivalent partial quadrature space.
   * 2. handle the generation of evaluators.
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorsFactory {
    /*!
     * \brief structure returned by the `analyseDependencies` method
     */
    struct LocalDependenciesAnalysisOutput {
      //! \brief missing dependencies
      std::vector<QPEvaluatorDescription> missingIPDependencies;
    };
    //! \brief a simple alias
    using Generator = std::shared_ptr<AbstractQPEvaluatorGenerator>;
    /*!
     * \brief describe a location
     * \return a description of the given location usable in an error message
     * \param[in] qspace: partial quadrature space on which the evaluator is
     * defined.
     * \param[in] s: stage in the time step
     */
    static std::string getLocationDescription(
        const std::shared_ptr<const PartialQuadratureSpace> qspace,
        const TimeStepStage s) noexcept;
    //! \brief default constructor
    QPEvaluatorsFactory() noexcept;
    // disabling copy and move constructors and assignment operators
    QPEvaluatorsFactory(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory(const QPEvaluatorsFactory &) = delete;
    QPEvaluatorsFactory &operator=(QPEvaluatorsFactory &&) = delete;
    QPEvaluatorsFactory &operator=(const QPEvaluatorsFactory &) = delete;
    /*!
     * \brief check if a generator has been declared
     * \return if a generator associated with the given partial quadrature
     * space, time step stage and name has been declared.
     *
     * \param[in] qspace: partial quadrature space
     * \param[in] s: stage in the time step
     * \param[in] n: name of the evaluator
     */
    [[nodiscard]] bool containsGenerator(
        const std::shared_ptr<const PartialQuadratureSpace> qspace,
        const TimeStepStage s,
        const std::string &n) const noexcept;
    //     /*
    //      * \brief register a new evaluator generator
    //      * \param[in, out] ctx: execution context
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      * \param[in] g: generator
    //      * \return true on success
    //      */
    //     [[nodiscard]] bool registerGenerator(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const std::string &,
    //         const Generator &) noexcept;
    //     /*
    //      * \brief generate an evaluator on the given partial quadrature space
    //      * \param[in, out] ctx: execution context
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] d: description of the evaluator
    //      * \return the evaluator, a null pointer on failure
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const QPEvaluatorDescription &) const noexcept;
    //     /*
    //      * \brief generate an evaluator on the given partial quadrature space
    //      * \param[in, out] ctx: execution context
    //      * \param[in] qspace: partial quadrature space
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      * \return the evaluator, a null pointer on failure
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const std::shared_ptr<const PartialQuadratureSpace>,
    //         const TimeStepStage,
    //         const std::string &) const noexcept;
    //     /*
    //      * \brief generate an evaluator on the given mesh set
    //      * \param[in, out] ctx: execution context
    //      * \param[in] m: mesh set on which the evaluator is defined
    //      * \param[in] s: stage in the time step
    //      * \param[in] d: description of the evaluator
    //      * \return the evaluator, a null pointer on failure
    //      *
    //      * \note the quadrature is unspecified. This call only works if the
    //      * evaluator is available for only one quadrature.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const size_type,
    //         const TimeStepStage,
    //         const QPEvaluatorDescription &) const noexcept;
    //     /*
    //      * \brief generate an evaluator on the given mesh set
    //      * \param[in, out] ctx: execution context
    //      * \param[in] m: mesh set on which the evaluator is defined
    //      * \param[in] s: stage in the time step
    //      * \param[in] n: name of the evaluator
    //      * \return the evaluator, a null pointer on failure
    //      *
    //      * \note the quadrature is unspecified. This call only works if the
    //      * evaluator is available for only one quadrature.
    //      */
    //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> generate(
    //         Context &,
    //         const size_type,
    //         const TimeStepStage,
    //         const std::string &) const noexcept;
    //     /*
    //      * \brief generate a set of evaluators on the given mesh set
    //      * \param[in, out] ctx: execution context
    //      * \param[in] m: mesh set on which the evaluators are defined
    //      * \param[in] s: stage in the time step
    //      * \param[in] descriptions: descriptions of the evaluators
    //      * \return the evaluators and their offsets, nothing on failure
    //      *
    //      * \note the quadrature is unspecified. This call only works if each
    //      * evaluator is available for only one quadrature.
    //      */
    //     [[nodiscard]] std::optional<
    //         std::vector<std::pair<size_type,
    //         std::shared_ptr<AbstractQPEvaluator>>>>
    //     generate(Context &,
    //              const size_type,
    //              const TimeStepStage,
    //              const std::vector<QPEvaluatorDescription> &) const noexcept;
    //     /*
    //      * \brief generate a set of evaluators on the given mesh set
    //      * \param[in, out] ctx: execution context
    //      * \param[in] m: mesh set on which the evaluators are defined
    //      * \param[in] s: stage in the time step
    //      * \param[in] names: names of the evaluators
    //      * \return the evaluators and their offsets, nothing on failure
    //      *
    //      * \note the quadrature is unspecified. This call only works if each
    //      * evaluator is available for only one quadrature.
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
    GeneratorsContainer &getGeneratorsContainer(const TimeStepStage s) noexcept;
    /*!
     * \return generators associated with the given stage in the time step
     * \param[in] s: stage in the time step
     */
    const GeneratorsContainer &getGeneratorsContainer(
        const TimeStepStage s) const noexcept;
    /*!
     * \return the list of registered generators for the given time step stage
     * \param[in] s: stage in the time step
     */
    [[nodiscard]] std::string getRegisteredGeneratorsList(
        const TimeStepStage s) const noexcept;
    //     /*
    //      * \brief check if the dependencies of the given evaluator are met
    //      * \param[in, out] ctx: execution context
    //      * \param[in] m: registered generators
    //      * \param[in] n: name of the evaluator
    //      * \param[in] previous_dependencies: evaluators already visited,
    //      * used to detect cycles
    //      * \return the missing dependencies
    //      */
    //     LocalDependenciesAnalysisOutput analyseDependencies(
    //         Context &,
    //         const std::map<std::string, Generator> &,
    //         const std::string &n,
    //         const std::vector<QPEvaluatorDescription> &) const noexcept;
    /*!
     * \brief list of registered generators computing values at the beginning of
     * the time step, sorted by location, partial quadrature space and name
     */
    GeneratorsContainer generatorsAtTheBeginningOfTheTimeStep;
    /*!
     * \brief list of registered generators computing values at the end of the
     * time step, sorted by location, partial quadrature space and name
     */
    GeneratorsContainer generatorsAtTheEndOfTheTimeStep;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORFACTORY_HXX */
