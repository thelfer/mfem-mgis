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
      ImmutablePartialQuadratureFunctionView() = default;

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          std::shared_ptr<const PartialQuadratureSpace> s,
          const size_type nv,
          const size_type db,
          const size_type ds)
      : qspace(s) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    this->data_stride = ds;
    this->data_begin = db;
    this->data_size = nv;
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
          std::span<const real> v,
          const size_type db,
          const size_type ds)
      : qspace(s) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    this->data_begin = db;
    this->data_size = ds;
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

  std::span<const real>
  ImmutablePartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  bool ImmutablePartialQuadratureFunctionView::checkCompatibility(
      const ImmutablePartialQuadratureFunctionView& v) const {
    if (this->getPartialQuadratureSpacePointer() !=
        v.getPartialQuadratureSpacePointer()) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  ImmutablePartialQuadratureFunctionView::
      ~ImmutablePartialQuadratureFunctionView() = default;

  PartialQuadratureFunction::PartialQuadratureFunction(
      PartialQuadratureFunction&& f, const bool local_copy) {
    if (!f.local_values_storage.empty()) {
      // the function holds the memory, just take it from him
      static_cast<PartialQuadratureFunctionDataLayout&>(*this).operator=(f);
      this->qspace = f.qspace;
      this->local_values_storage = std::move(f.local_values_storage);
      this->values = local_values_storage;
      this->immutable_values = local_values_storage;
    } else {
      // the function does not hold the memory
      if (local_copy) {
        this->copy(f);
      } else {
        this->makeView(f);
      }
    }
  }  // end of PartialQuadratureFunction

  PartialQuadratureFunction::PartialQuadratureFunction(
      const PartialQuadratureFunction& v)
      : ImmutablePartialQuadratureFunctionView() {
    this->copy(v);
  }  // end of PartialQuadratureFunction

  PartialQuadratureFunction::PartialQuadratureFunction(
      const ImmutablePartialQuadratureFunctionView& v) {
    this->copy(v);
  }  // end of ImmutablePartialQuadratureFunctionView

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s, const size_type nv)
      : ImmutablePartialQuadratureFunctionView(s, nv, 0, nv) {
    this->local_values_storage.resize(
        this->qspace->getNumberOfIntegrationPoints() * this->data_size);
    this->values = std::span<real>(this->local_values_storage);
    this->immutable_values = std::span<const real>(this->local_values_storage);
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::span<real> v,
      const size_type db,
      const size_type ds)
      : ImmutablePartialQuadratureFunctionView(s, v, db, ds),
        values(v) {
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  void PartialQuadratureFunction::makeView(PartialQuadratureFunction& f) {
    static_cast<PartialQuadratureFunctionDataLayout&>(*this).operator=(f);
    this->qspace = f.qspace;
    this->values = f.values;
    this->immutable_values = f.immutable_values;
  }

  PartialQuadratureFunction& PartialQuadratureFunction::operator=(
      const PartialQuadratureFunction& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  PartialQuadratureFunction& PartialQuadratureFunction::operator=(
      const ImmutablePartialQuadratureFunctionView& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  void PartialQuadratureFunction::copy(
      const ImmutablePartialQuadratureFunctionView& v) {
    this->qspace = v.getPartialQuadratureSpacePointer();
    const auto n = this->qspace->getNumberOfIntegrationPoints();
    this->data_begin = size_type{};
    this->data_size = v.getNumberOfComponents();
    this->data_stride = v.getNumberOfComponents();
    this->local_values_storage.resize(this->data_size * n);
    this->values = local_values_storage;
    this->immutable_values = local_values_storage;
    this->copyValues(v);
  }  // end of copy

  void PartialQuadratureFunction::copyValues(
      const ImmutablePartialQuadratureFunctionView& v) {
    const auto* const v_values = v.getValues().data() + v.getDataOffset();
    const auto vs = v.getDataStride();
    if (vs == v.getNumberOfComponents()) {
      // data are continous in v
      std::copy(v_values, v_values + this->values.size(), this->values.begin());
    } else {
      if (this->data_size == 1) {
        // special case for scalars
        for (size_type i = 0; i != this->values.size(); ++i) {
          this->values[i] = v_values[i * vs];
        }
      } else {
        const auto n =
            this->getPartialQuadratureSpace().getNumberOfIntegrationPoints();
        auto pv = this->values.begin();
        for (size_type i = 0; i != n; ++i) {
          const auto b = v_values + i * vs;
          const auto e = b + this->data_size;
          std::copy(b, e, pv);
          pv += this->data_size;
        }
      }
    }
  }  // end of copy

  real& PartialQuadratureFunction::getIntegrationPointValue(const size_type e,
                                                            const size_type i) {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  std::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

}  // end of namespace mfem_mgis
