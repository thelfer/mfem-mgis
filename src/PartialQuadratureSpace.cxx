/*!
 * \file   src/PartialQuadratureSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <iterator>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static size_type buildPartialQuadratureSpaceOffsets(
      std::unordered_map<size_type, size_type>& offsets,
      std::unordered_map<size_type, size_type>& number_of_quadrature_points,
      const FiniteElementSpace<parallel>& fespace,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>&
          integration_rule_selector) {
    auto ng = size_type{};
    for (size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(fespace.GetFE(i));
      const auto& tr = *(fespace.GetElementTransformation(i));
      offsets[i] = ng;
      const auto& ir = integration_rule_selector(fe, tr);
      const auto lng = ir.GetNPoints();
      ng += lng;
      number_of_quadrature_points[i] = lng;
    }
    return ng;
  }  // end of buildPartialQuadratureSpaceOffsets

  void PartialQuadratureSpace::treatInvalidElementIndex(const size_type id,
                                                  const size_type i) {
    mgis::raise(
        "PartialQuadratureSpace::getOffset: "
        "invalid element index '" +
        std::to_string(i) + "' for material '" + std::to_string(id) + "'");
  }  // end of treatInvalidElementIndex

  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>& irs)
      : fe_discretization(fed), integration_rule_selector(irs), id(m) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& fespace =
          this->fe_discretization.getFiniteElementSpace<true>();
      this->ng = buildPartialQuadratureSpaceOffsets<true>(
          this->offsets, this->number_of_quadrature_points, fespace, this->id,
          this->integration_rule_selector);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto& fespace =
          this->fe_discretization.getFiniteElementSpace<false>();
      this->ng = buildPartialQuadratureSpaceOffsets<false>(
          this->offsets, this->number_of_quadrature_points, fespace, this->id,
          this->integration_rule_selector);
    }
  }  // end of PartialQuadratureSpace

  const mfem::IntegrationRule& PartialQuadratureSpace::getIntegrationRule(
      const mfem::FiniteElement& e,
      const mfem::ElementTransformation& tr) const {
    return this->integration_rule_selector(e, tr);
  }

  const FiniteElementDiscretization&
  PartialQuadratureSpace::getFiniteElementDiscretization() const {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  PartialQuadratureSpace::~PartialQuadratureSpace() = default;

}  // end of namespace mfem_mgis
