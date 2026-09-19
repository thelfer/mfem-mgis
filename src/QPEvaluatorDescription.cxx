/*!
 * \file   MFEMMGIS/QPEvaluatorDescription.cxx
 * \brief  This file implements the methods of the `QPEvaluatorDescription`
 * class \author Thomas Helfer \date   04/11/2022
 */

#include "MFEMMGIS/QPEvaluator/QPEvaluatorDescription.hxx"

namespace mfem_mgis {

  QPEvaluatorDescription::QPEvaluatorDescription(std::string n,
                                                 const size_type nc,
                                                 const bool b)
      : name(n),
        number_of_components(nc),
        shall_not_be_used_in_evaluators_generation(b) {
    if (n.empty()) {
      raise("invalid name");
    }
    if (this->number_of_components < 1) {
      raise("invalid number of components");
    }
  }  // end of QPEvaluatorDescription

  std::string QPEvaluatorDescription::getName() const noexcept {
    return this->name;
  }  // end of getName

  size_type QPEvaluatorDescription::getNumberOfComponents() const noexcept {
    return this->number_of_components;
  }  // end of getNumRows

  bool QPEvaluatorDescription::shallNotBeUsedInEvaluatorsGeneration()
      const noexcept {
    return this->shall_not_be_used_in_evaluators_generation;
  }  // end of shallNotBeUsedInEvaluatorsGeneration

  QPEvaluatorDescription::~QPEvaluatorDescription() noexcept = default;

}  // namespace mfem_mgis