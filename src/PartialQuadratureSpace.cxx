/*!
 * \file   src/PartialQuadratureSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <iostream>
#include <iterator>
#include <cmath>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  PartialQuadratureSpace::PartialQuadratureSpace(
      const mfem::FiniteElementSpace& fs,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>&
          integration_rule_selector)
      : fespace(fs), id(m) {
    this->ng = size_type{};
    for (size_type i = 0; i != this->fespace.GetNE(); ++i) {
      if (this->fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(this->fespace.GetFE(i));
      const auto& tr = *(this->fespace.GetElementTransformation(i));
      this->offsets[i] = this->ng;
      const auto& ir = integration_rule_selector(fe, tr);
      this->ng += ir.GetNPoints();
    }
  }  // end of PartialQuadratureSpace::PartialQuadratureSpace

  size_type PartialQuadratureSpace::getNumberOfElements() const {
    return this->offsets.size();
  }  // end of PartialQuadratureSpace::getNumberOfElements

  size_type PartialQuadratureSpace::getNumberOfIntegrationPoints() const {
    return this->ng;
  }  // end of PartialQuadratureSpace::getNumberOfIntegrationPoints

  size_type PartialQuadratureSpace::getOffset(const size_type i) const {
    const auto p = this->offsets.find(i);
    if (p == this->offsets.end()) {
      mfem::mfem_error(
          "PartialQuadratureSpace::getOffset: "
          "invalid number of elements");
    }
    return p->second;
  }  // end of PartialQuadratureSpace::getOffset

  PartialQuadratureSpace::~PartialQuadratureSpace() = default;

}  // end of namespace mfem_mgis
