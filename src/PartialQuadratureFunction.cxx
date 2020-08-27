/*!
 * \file   src/PartialQuadratureFunction.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s, const size_type nv)
      : qspace(std::move(s)), data_size(nv) {
    this->values.resize(this->qspace->getNumberOfElements() * this->data_size);
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  mgis::span<real> PartialQuadratureFunction::getIntegrationPointValues(const size_type e,
                                                                        const size_type i) {
    const auto o = this->qspace->getOffset(e) + i * (this->data_size);
    return mgis::span<real>(this->values.data() + o, this->data_size);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  mgis::span<const real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) const {
    const auto o = this->qspace->getOffset(e) + i * (this->data_size);
    return mgis::span<const real>(this->values.data() + o, this->data_size);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

}  // end of namespace mfem_mgis
