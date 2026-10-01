/*!
 * \file   MFEMMGIS/QPEvaluator/AbstractQPEvaluatorGenerator.hxx
 * \brief  This file declares the `AbstractQPEvaluatorGenerator` class
 * \author Thomas Helfer
 * \date   20/10/2022
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATORGENERATOR_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATORGENERATOR_HXX

#include <memory>
#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct AbstractQPEvaluator;
  struct QPEvaluatorsFactory;
  struct QPEvaluatorDescription;
  struct PartialQuadratureSpace;
  struct AbstractQPEvaluatorGenerator;

  /*!
   * \brief this class is meant to generate integration
   * point value evaluators.
   *
   * Such generator is meant to be used as follows:
   *
   * 1. They are registered in an instance of `QPEvaluatorsFactory`
   * 2. They are called through the `generate` method of the
   *    `QPEvaluatorsFactory` class
   */
  struct MFEM_MGIS_EXPORT AbstractQPEvaluatorGenerator {
    //! \returns the number of components of the evaluator
    [[nodiscard]] virtual size_type getNumberOfComponents() const = 0;
    //! \return the dependencies of the generated evaluator
    [[nodiscard]] virtual std::vector<QPEvaluatorDescription> getDependencies()
        const = 0;
    /*!
     * \brief generate an evaluator
     * \return the evaluator generated
     * \param[in, out] ctx: execution context
     * \param[in] f: factory used to retrieve the dependencies of the generated
     * evaluator.
     * \param[in] qspace: partial quadrature space
     * \param[in] s: stage in the time step
     */
    [[nodiscard]] virtual std::shared_ptr<AbstractQPEvaluator> operator()(
        Context& ctx,
        const QPEvaluatorsFactory& f,
        const PartialQuadratureSpace& qspace,
        const TimeStepStage s) const = 0;
    //! \brief destructor
    virtual ~AbstractQPEvaluatorGenerator();
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATORGENERATOR_HXX */