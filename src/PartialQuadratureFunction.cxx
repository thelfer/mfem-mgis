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
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static std::shared_ptr<PartialQuadratureFunction>
  PartialQuadratureFunction_evaluate(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<real(const mfem::FiniteElement&,
                         mfem::ElementTransformation&)> f) {
    const auto& fed = s->getFiniteElementDiscretization();
    const auto& fespace = fed.getFiniteElementSpace<parallel>();
    const auto m = s->getId();
    auto values = std::make_shared<PartialQuadratureFunction>(s, 1);
    for (size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(fespace.GetFE(i));
      auto& tr = *(fespace.GetElementTransformation(i));
      const auto& ir = s->getIntegrationRule(fe, tr);
      for (size_type g = 0; g != ir.GetNPoints(); ++g) {
        // get the gradients of the shape functions
        const auto& ip = ir.IntPoint(g);
        tr.SetIntPoint(&ip);
        values->getIntegrationPointValue(i, g) = f(fe, tr);
      }
    }
    return values;
  }  // end of PartialQuadratureFunction_evaluate

  std::shared_ptr<PartialQuadratureFunction>
  PartialQuadratureFunction::evaluate(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<real(const mfem::FiniteElement&,
                         mfem::ElementTransformation&)> f) {
    const auto& fed = s->getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return PartialQuadratureFunction_evaluate<true>(s, f);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return PartialQuadratureFunction_evaluate<false>(s, f);
  }  // end of evaluate

  std::shared_ptr<PartialQuadratureFunction>
  PartialQuadratureFunction::evaluate(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<real(real, real)> f) {
    const auto n = getSpaceDimension(s->getFiniteElementDiscretization());
    if (n != 2) {
      raise("PartialQuadratureFunction::evaluate: invalid space dimension");
    }
    return PartialQuadratureFunction::evaluate(
        s, [&f](const mfem::FiniteElement&, mfem::ElementTransformation& tr) {
          mfem::Vector p;
          tr.Transform(tr.GetIntPoint(), p);
          return f(p[0], p[1]);
        });
  }  // end of evaluate

  std::shared_ptr<PartialQuadratureFunction>
  PartialQuadratureFunction::evaluate(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<real(real, real, real)> f) {
    const auto n = getSpaceDimension(s->getFiniteElementDiscretization());
    if (n != 3) {
      raise("PartialQuadratureFunction::evaluate: invalid space dimension");
    }
    return PartialQuadratureFunction::evaluate(
        s, [&f](const mfem::FiniteElement&, mfem::ElementTransformation& tr) {
          mfem::Vector p;
          tr.Transform(tr.GetIntPoint(), p);
          return f(p[0], p[1], p[2]);
        });
  }  // end of evaluate

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          std::shared_ptr<const PartialQuadratureSpace> s,
          const size_type nv,
          const size_type db,
          const size_type ds)
      : qspace(s), data_stride(ds), data_begin(db), data_size(nv) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    if (this->data_begin < 0) {
      raise(
          "ImmutablePartialQuadratureFunctionView::"
          "ImmutablePartialQuadratureFunctionView: "
          "invalid start of the data");
    }
    if (this->data_size <= 0) {
      raise(
          "ImmutablePartialQuadratureFunctionView::"
          "ImmutablePartialQuadratureFunctionView: "
          "invalid data size");
    }
    if (this->data_begin + this->data_size > this->data_stride) {
      raise(
          "ImmutablePartialQuadratureFunctionView::"
          "ImmutablePartialQuadratureFunctionView: "
          "invalid data range is outside the stride size");
    }
  }  // end of ImmutablePartialQuadratureFunctionView

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          std::shared_ptr<const PartialQuadratureSpace> s,
          mgis::span<const real> v,
          const size_type db,
          const size_type ds)
      : qspace(s), data_begin(db), data_size(ds) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    if (this->data_begin < 0) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "start of the data");
    }
    if (ds < 0) {
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
    if (this->data_begin >= this->data_stride) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "start of the data");
    }
    if (ds == std::numeric_limits<size_type>::max()) {
      this->data_size = this->data_stride - this->data_begin;
    } else {
      if (this->data_begin + ds > this->data_stride) {
        raise(
            "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
            "data range is outside the stride size");
      }
      this->data_size = ds;
    }
    this->immutable_values = v;
  }  // end of ImmutablePartialQuadratureFunctionView

  const real& ImmutablePartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  mgis::span<const real>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  ImmutablePartialQuadratureFunctionView::
      ~ImmutablePartialQuadratureFunctionView() = default;

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s, const size_type nv)
      : ImmutablePartialQuadratureFunctionView(s, nv, 0, nv) {
    this->local_values_storage.resize(
        this->qspace->getNumberOfIntegrationPoints() * this->data_size);
    this->values = mgis::span<real>(this->local_values_storage);
    this->immutable_values = mgis::span<const real>(this->local_values_storage);
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s,
      mgis::span<real> v,
      const size_type db,
      const size_type ds)
      : ImmutablePartialQuadratureFunctionView(s, v, db, ds),
        values(v) {
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  real& PartialQuadratureFunction::getIntegrationPointValue(const size_type e,
                                                            const size_type i) {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  mgis::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

}  // end of namespace mfem_mgis
