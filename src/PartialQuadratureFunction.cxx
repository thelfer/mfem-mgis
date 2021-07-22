/*!
 * \file   src/PartialQuadratureFunction.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s, const size_type nv)
      : qspace(s), data_stride(data_size), data_begin(0), data_size(nv) {
    if (this->data_size <= 0) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "values size");
    }
    this->values_storage.resize(this->qspace->getNumberOfIntegrationPoints() *
                                this->data_size);
    this->values = mgis::span<real>(this->values_storage.data(),
                                    this->values_storage.size());
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s,
      mgis::span<real> v,
      const size_type db,
      const size_type ds)
      : qspace(s), values(v), data_begin(db) {
    if (db < 0) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "start of the data");
    }
    if (ds <= 0) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "data size");
    }
    const auto d = std::div(static_cast<size_type>(v.size()),
                            this->qspace->getNumberOfIntegrationPoints());
    if ((d.rem != 0) || (d.quot <= 0)) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "values size");
    }
    this->data_stride = d.quot;
    if (db > this->data_stride) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "start of the data");
    }
    if (ds == std::numeric_limits<size_type>::max()) {
      this->data_size = this->data_stride;
    } else {
      if (db + ds >= this->data_stride) {
        raise(
            "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
            "data range is outside the stride size");
      }
      this->data_size = ds;
    }
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  real& PartialQuadratureFunction::getIntegrationPointValue(const size_type e,
                                                            const size_type i) {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  const real& PartialQuadratureFunction::getIntegrationPointValue(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  mgis::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  mgis::span<const real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

}  // end of namespace mfem_mgis
