/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluatorDescription.hxx
 * \brief  This file declares the `QPEvaluatorDescription` class
 * \author Thomas Helfer
 * \date   04/11/2022
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORDESCRIPTION_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORDESCRIPTION_HXX

#include <string>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief this class is meant to describe an QPEvaluator.
   * A typicaly usage of this class is to describe a set of evaluators to be
   * passed to the `QPEvaluatorsSet` class.
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorDescription {
    /*!
     * \brief constructor
     * \param[in] n: name of the dependency
     * \param[in] nc: expected number of components
     * \param[in] b: flag stating if this description shall not be used to
     * generate an evaluator. See the description of the
     * `shallNotBeUsedInEvaluatorsGeneration` method.
     */
    QPEvaluatorDescription(std::string, const size_type, const bool = false);
    //! \return the name of the dependency
    [[nodiscard]] std::string getName() const noexcept;
    //! \return the expected number of components to be computed by the
    //! dependency
    [[nodiscard]] size_type getNumberOfComponents() const noexcept;
    /*!
     * \return if this description shall be used to generate an evaluator
     *
     * This flag has been introduced because some formulations may want to
     * evaluate some external state variables internally. A typical example is
     * the temperature which is evaluated internally by the
     * `ImplicitHeatTransferSolidMaterialFormulation` and passed as the first
     * external state variable to the heat transfer behavior (the fact that the
     * temperature is passed as an external state variable is an `MFront`
     * convention).
     */
    bool shallNotBeUsedInEvaluatorsGeneration() const noexcept;
    //! \brief destructor
    ~QPEvaluatorDescription() noexcept;

   private:
    //! \brief name of the dependency
    const std::string name;
    //! \brief expected number of rows
    const size_type number_of_components;
    //! \brief flag stating if an evaluator is just a placeholder
    const bool shall_not_be_used_in_evaluators_generation;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORDESCRIPTION_HXX */
