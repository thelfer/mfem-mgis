/*!
 * \file   src/PartialQuadratureFunction.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "mfem/mesh/submesh/submesh.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/submesh/psubmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/fem/pgridfunc.hpp"
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
    this->data_begin = db;
    this->data_size = ds;
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_stride = 0;
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
    if (this->data_begin + this->data_size > this->data_stride) {
      raise(
          "PartialQuadratureFunction::PartialQuadratureFunction: invalid "
          "data range is outside the stride size");
    }
    this->immutable_values = v;
  }  // end of ImmutablePartialQuadratureFunctionView

  const real* ImmutablePartialQuadratureFunctionView::data(
      const size_type e, const size_type i) const {
    return this->data(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

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

  real* PartialQuadratureFunction::data(const size_type e, const size_type i) {
    return this->data(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  real& PartialQuadratureFunction::getIntegrationPointValue(const size_type e,
                                                            const size_type i) {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  std::span<real> PartialQuadratureFunction::getIntegrationPointValues(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of PartialQuadratureFunction::getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

  struct PartialQuadratureFunctionsCoefficientBase {
    //
    PartialQuadratureFunctionsCoefficientBase(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
      if (fcts.empty()) {
        raise("no functions defined");
      }
      const auto n = fcts.at(0).getNumberOfComponents();
      for (const auto& f : fcts) {
        const auto mid = f.getPartialQuadratureSpace().getId();
        if (!this->functions.insert({mid, f}).second) {
          raise("multiple functions defined for material '" +
                std::to_string(mid) + "'");
        }
        if (n != f.getNumberOfComponents()) {
          raise("inconsistent number of components");
        }
      }
    }
    //
    PartialQuadratureFunctionsCoefficientBase(
        PartialQuadratureFunctionsCoefficientBase&&) = default;
    PartialQuadratureFunctionsCoefficientBase(
        const PartialQuadratureFunctionsCoefficientBase&) = default;
    PartialQuadratureFunctionsCoefficientBase& operator=(
        PartialQuadratureFunctionsCoefficientBase&&) = default;
    PartialQuadratureFunctionsCoefficientBase& operator=(
        const PartialQuadratureFunctionsCoefficientBase&) = default;
    ~PartialQuadratureFunctionsCoefficientBase() = default;

   protected:
    std::unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        functions;
  };  // end of PartialQuadratureFunctionsCoefficientBase

  struct PartialQuadratureFunctionsScalarCoefficient
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::Coefficient {
    //
    PartialQuadratureFunctionsScalarCoefficient(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(fcts) {
      for (const auto& [mid, f] : functions) {
        static_cast<void>(mid);
        if (f.getNumberOfComponents() != 1) {
          raise("non scalar function given");
        }
      }
    }
    //
    PartialQuadratureFunctionsScalarCoefficient(
        PartialQuadratureFunctionsScalarCoefficient&&) = default;
    PartialQuadratureFunctionsScalarCoefficient(
        const PartialQuadratureFunctionsScalarCoefficient&) = default;
    PartialQuadratureFunctionsScalarCoefficient& operator=(
        PartialQuadratureFunctionsScalarCoefficient&&) = default;
    PartialQuadratureFunctionsScalarCoefficient& operator=(
        const PartialQuadratureFunctionsScalarCoefficient&) = default;
    //
    double Eval(mfem::ElementTransformation& tr,
                const mfem::IntegrationPoint& i) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        return 0.;
      }
      return p->second.getIntegrationPointValue(tr.ElementNo, i.index);
    }  // end of Eval
  };

  struct PartialQuadratureFunctionsVectorCoefficient
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::VectorCoefficient {
    //
    PartialQuadratureFunctionsVectorCoefficient(
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(fcts),
          mfem::VectorCoefficient(fcts.at(0).getNumberOfComponents()) {}
    //
    PartialQuadratureFunctionsVectorCoefficient(
        PartialQuadratureFunctionsVectorCoefficient&&) = default;
    PartialQuadratureFunctionsVectorCoefficient(
        const PartialQuadratureFunctionsVectorCoefficient&) = default;
    PartialQuadratureFunctionsVectorCoefficient& operator=(
        PartialQuadratureFunctionsVectorCoefficient&&) = default;
    PartialQuadratureFunctionsVectorCoefficient& operator=(
        const PartialQuadratureFunctionsVectorCoefficient&) = default;
    //
    void Eval(mfem::Vector& values,
              mfem::ElementTransformation& tr,
              const mfem::IntegrationPoint& ip) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        values = 0.;
      } else {
        const auto rvalues =
            p->second.getIntegrationPointValues(tr.ElementNo, ip.index);
        for (size_type i = 0; i != this->GetVDim(); ++i) {
          values[i] = rvalues[i];
        }
      }
    }  // end of Eval
  };

  template <bool parallel>
  static std::pair<std::unique_ptr<FiniteElementSpace<parallel>>,
                   std::unique_ptr<GridFunction<parallel>>>
  makeGridFunction_impl(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    if (fcts.empty()) {
      raise("no functions defined");
    }
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    auto& fes = fed.getFiniteElementSpace<parallel>();
    auto& mesh = fed.getMesh<parallel>();
    auto fespace = std::make_unique<FiniteElementSpace<parallel>>(
        const_cast<Mesh<parallel>*>(&mesh), fes.FEColl(), n, fes.GetOrdering());
    auto f = std::make_unique<GridFunction<parallel>>(fespace.get());
    return {std::move(fespace), std::move(f)};
  }

  template <>
  std::pair<std::unique_ptr<FiniteElementSpace<true>>,
            std::unique_ptr<GridFunction<true>>>
  makeGridFunction<true>(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    return makeGridFunction_impl<true>(fcts);
  }

  template <>
  std::pair<std::unique_ptr<FiniteElementSpace<false>>,
            std::unique_ptr<GridFunction<false>>>
  makeGridFunction<false>(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    return makeGridFunction_impl<false>(fcts);
  }

  template <bool parallel>
  static std::pair<std::unique_ptr<FiniteElementSpace<parallel>>,
                   std::unique_ptr<GridFunction<parallel>>>
  makeGridFunction_impl(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<parallel>>& mesh) {
    if (fcts.empty()) {
      raise("no functions defined");
    }
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    auto& fes = fed.getFiniteElementSpace<parallel>();
    auto fespace = std::make_unique<FiniteElementSpace<parallel>>(
        mesh.get(), fes.FEColl(), n, fes.GetOrdering());
    auto f = std::make_unique<GridFunction<parallel>>(fespace.get());
    return {std::move(fespace), std::move(f)};
  }

  template <>
  std::pair<std::unique_ptr<FiniteElementSpace<true>>,
            std::unique_ptr<GridFunction<true>>>
  makeGridFunction<true>(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<true>>& mesh) {
    return makeGridFunction_impl<true>(fcts, mesh);
  }

  template <>
  std::pair<std::unique_ptr<FiniteElementSpace<false>>,
            std::unique_ptr<GridFunction<false>>>
  makeGridFunction<false>(
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<false>>& mesh) {
    return makeGridFunction_impl<false>(fcts, mesh);
  }

  template <bool parallel>
  static void updateGridFunction_impl(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    const auto& fes = fed.getFiniteElementSpace<parallel>();
    const auto& fespace = f.FESpace();
    if ((fespace->GetMesh() != &(mesh)) ||  //
        (fespace->GetVDim() != n) ||        //
        (fes.FEColl() != fespace->FEColl()) ||
        (fes.GetOrdering() != fespace->GetOrdering())) {
      raise("inconsistent grid function");
    }
    if (n == 1u) {
      auto c = PartialQuadratureFunctionsScalarCoefficient(fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c = PartialQuadratureFunctionsVectorCoefficient(fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    updateGridFunction_impl<true>(f, fcts);
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    updateGridFunction_impl<false>(f, fcts);
  }

  template <bool parallel>
  static void updateGridFunction_impl(
      GridFunction<parallel> & f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<parallel>>& mesh) {
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& fes = fed.getFiniteElementSpace<parallel>();
    const auto& fespace = f.FESpace();
    if ((fespace->GetMesh() != mesh.get()) ||  //
        (fespace->GetVDim() != n) ||           //
        (fes.FEColl() != fespace->FEColl()) ||
        (fes.GetOrdering() != fespace->GetOrdering())) {
      raise("inconsistent grid function");
    }
    if (n == 1u) {
      auto c = PartialQuadratureFunctionsScalarCoefficient(fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c = PartialQuadratureFunctionsVectorCoefficient(fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<true>>& mesh) {
    updateGridFunction_impl<true>(f, fcts, mesh);
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const std::shared_ptr<SubMesh<false>>& mesh) {
    updateGridFunction_impl<false>(f, fcts, mesh);
  }


}  // end of namespace mfem_mgis
