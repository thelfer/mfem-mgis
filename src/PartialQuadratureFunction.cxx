/*!
 * \file   src/PartialQuadratureFunction.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <set>
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
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    const auto& fespace = s->getFiniteElementSpace<parallel>(ctx) | or_raise;
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

  static void checkPartialQuadratureFunctionConstructorArguments(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::span<const real> v,
      const size_type db,
      const size_type ds) {
    if (s->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      return;
    }
    if (db < 0) {
      raise("invalid start of the data");
    }
    if (ds < 0) {
      raise("invalid data size");
    }
    const auto d = std::div(static_cast<size_type>(v.size()),
                            s->getNumberOfIntegrationPoints());
    if ((d.rem != 0) || (d.quot <= 0)) {
      raise("invalid values size");
    }
    if (db >= d.quot) {
      raise("invalid start of the data");
    }
    if (db + ds > d.quot) {
      raise("data range is outside the stride size");
    }
  }  // end of checkPartialQuadratureFunctionConstructorArguments

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView() = default;

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          std::shared_ptr<const PartialQuadratureSpace> s,
          const size_type nv,
          const size_type db,
          const size_type ds)
      : qspace(s) {
    if (s.get() == nullptr) {
      raise("invalid partial quadrature space pointer");
    }
    this->data_stride = ds;
    this->data_begin = db;
    this->data_size = nv;
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_stride = 0;
      return;
    }
    if (this->data_begin < 0) {
      raise("invalid start of the data");
    }
    if (this->data_size <= 0) {
      raise("invalid data size");
    }
    if (this->data_begin + this->data_size > this->data_stride) {
      raise("invalid data range is outside the stride size");
    }
  }  // end of ImmutablePartialQuadratureFunctionView

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          std::shared_ptr<const PartialQuadratureSpace> s,
          std::span<const real> v,
          const size_type db,
          const size_type ds)
      : qspace(s) {
    if (s.get() == nullptr) {
      raise("invalid partial quadrature space pointer");
    }
    this->data_begin = db;
    this->data_stride = ds;
    this->data_size = ds;
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_stride = 0;
      return;
    }
    checkPartialQuadratureFunctionConstructorArguments(s, v, db, ds);
#pragma message("HERE")
    const auto d = std::div(static_cast<size_type>(v.size()),
                            this->qspace->getNumberOfIntegrationPoints());
    this->data_stride = d.quot;
    this->immutable_values = v;
  }  // end of ImmutablePartialQuadratureFunctionView

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          ImmutablePartialQuadratureFunctionView&&) noexcept = default;

  ImmutablePartialQuadratureFunctionView::
      ImmutablePartialQuadratureFunctionView(
          const ImmutablePartialQuadratureFunctionView&) noexcept = default;

  ImmutablePartialQuadratureFunctionView&
  ImmutablePartialQuadratureFunctionView::operator=(
      ImmutablePartialQuadratureFunctionView&&) noexcept = default;

  ImmutablePartialQuadratureFunctionView&
  ImmutablePartialQuadratureFunctionView::operator=(
      const ImmutablePartialQuadratureFunctionView&) noexcept = default;

  const real* ImmutablePartialQuadratureFunctionView::data(
      const size_type e, const size_type i) const {
    return this->data(this->qspace->getOffset(e) + i);
  }  // end of data

  const real& ImmutablePartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type e, const size_type i) const {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValue

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
      PartialQuadratureFunction&& f) {
    if (!f.local_values_storage.empty()) {
      // the function holds the memory, just take it from him
      static_cast<PartialQuadratureFunctionDataLayout&>(*this).operator=(f);
      this->qspace = f.qspace;
      this->local_values_storage = std::move(f.local_values_storage);
      this->mutable_values = local_values_storage;
      this->immutable_values = local_values_storage;
    } else {
      // the function does not hold the memory
      this->makeView(f);
    }
  }  // end of PartialQuadratureFunction

  [[nodiscard]] std::optional<PartialQuadratureFunction>
  PartialQuadratureFunction::copy(
      Context& ctx, const ImmutablePartialQuadratureFunctionView& v) noexcept {
    try {
      return {PartialQuadratureFunction(v)};
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of ImmutablePartialQuadratureFunctionView

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s, const size_type nv)
      : PartialQuadratureFunctionView(s, nv, 0, nv) {
    if (s.get() == nullptr) {
      raise("invalid partial quadrature space pointer");
    }
    this->local_values_storage.resize(
        this->qspace->getNumberOfIntegrationPoints() * this->data_size);
    this->mutable_values = std::span<real>(this->local_values_storage);
    this->immutable_values = std::span<const real>(this->local_values_storage);
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  std::optional<PartialQuadratureFunction> PartialQuadratureFunction::borrow(
      Context& ctx,
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::span<real> v,
      const size_type db,
      const size_type ds) noexcept {
    try {
      return {PartialQuadratureFunction(s, StorageMode::EXTERNAL_STORAGE, v, db,
                                        ds)};
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of borrow

  PartialQuadratureFunction::PartialQuadratureFunction(
      std::shared_ptr<const PartialQuadratureSpace> s,
      const StorageMode sm,
      std::span<real> v,
      const size_type db,
      const size_type ds) {
    if (s.get() == nullptr) {
      raise("invalid partial quadrature space pointer");
    }
    this->qspace = s;
    //
    this->data_begin = db;
    this->data_size = ds;
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_stride = 0;
      return;
    }
    checkPartialQuadratureFunctionConstructorArguments(s, v, db, ds);
    const auto d = std::div(static_cast<size_type>(v.size()),
                            this->qspace->getNumberOfIntegrationPoints());
    this->data_stride = d.quot;
    //
    if (sm == StorageMode::EXTERNAL_STORAGE) {
      this->mutable_values = v;
      this->immutable_values = std::span<const real>(v);
    } else {
      this->local_values_storage.resize(v.size());
      std::copy(v.begin(), v.end(), this->local_values_storage.begin());
      this->mutable_values = std::span<real>(this->local_values_storage);
      this->immutable_values =
          std::span<const real>(this->local_values_storage);
    }
  }  // end of PartialQuadratureFunction::PartialQuadratureFunction

  void PartialQuadratureFunction::makeView(PartialQuadratureFunction& f) {
    static_cast<PartialQuadratureFunctionDataLayout&>(*this).operator=(f);
    this->qspace = f.qspace;
    this->mutable_values = f.mutable_values;
    this->immutable_values = f.immutable_values;
  }

  PartialQuadratureFunction::PartialQuadratureFunction(
      const ImmutablePartialQuadratureFunctionView& v) {
    this->qspace = v.getPartialQuadratureSpacePointer();
    const auto n = this->qspace->getNumberOfIntegrationPoints();
    this->data_begin = size_type{};
    this->data_size = v.getNumberOfComponents();
    this->data_stride = v.getNumberOfComponents();
    this->local_values_storage.resize(this->data_size * n);
    this->mutable_values = local_values_storage;
    this->immutable_values = local_values_storage;
    this->copyValues(v);
  }  // end of copy

  void PartialQuadratureFunction::copyValues(
      const ImmutablePartialQuadratureFunctionView& v) {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    assign_values(ctx, *this, v) | or_die;
  }  // end of copy

  MFEM_MGIS_EXPORT [[nodiscard]] bool assign_values(
      Context& ctx,
      PartialQuadratureFunctionView f,
      const ImmutablePartialQuadratureFunctionView& v) noexcept {
    //
    auto& qspace = f.getPartialQuadratureSpace();
    //
    if (qspace.getId() != v.getPartialQuadratureSpace().getId()) {
      return ctx.registerErrorMessage("unmatched material");
    }
    //
    const auto n = getSpaceSize(qspace);
    if (n != getSpaceSize(v.getPartialQuadratureSpace())) {
      return ctx.registerErrorMessage("unmatched space size");
    }
    const auto nc = f.getNumberOfComponents();
    if (nc != v.getNumberOfComponents()) {
      return ctx.registerErrorMessage("unmatched number of components");
    }
    //
    if (getSpaceSize(qspace) == 0) {
      return true;
    }
    //
    auto* f_values = f.getValues().data() + f.getDataOffset();
    const auto* const v_values = v.getValues().data() + v.getDataOffset();
    const auto fs = f.getDataStride();
    const auto vs = v.getDataStride();
    const auto f_data_continuous = fs == nc;
    if (f_data_continuous) {
      ctx.assertOrTerminate(f.getDataOffset() == 0,
                            "inconsistent function view, offset shall be null");
      ctx.assertOrTerminate(v.getDataOffset() == 0,
                            "inconsistent function view, offset shall be null");
      if (vs == v.getNumberOfComponents()) {
        // data are also continous in v
        std::copy(v_values, v_values + n, f_values);
      } else {
        if (nc == 1) {
          // special case for scalars
          for (size_type i = 0; i != n; ++i) {
            f_values[i] = v_values[i * vs];
          }
        } else {
          auto pvalues = f_values;
          for (size_type i = 0; i != n; ++i) {
            const auto b = v_values + i * vs;
            const auto e = b + nc;
            std::copy(b, e, pvalues);
            pvalues += nc;
          }
        }
      }
    } else {
      if (nc == 1) {
        for (size_type i = 0; i != n; ++i) {
          f_values[i * fs] = v_values[i * vs];
        }
      } else {
        for (size_type i = 0; i != n; ++i) {
          const auto b = v_values + i * vs;
          const auto e = b + nc;
          std::copy(b, e, f_values + i * fs);
        }
      }
    }
    return true;
  }  // end of assign_values

  real* PartialQuadratureFunctionView::data(const size_type e,
                                            const size_type i) {
    return this->data(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  real& PartialQuadratureFunctionView::getIntegrationPointValue(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValue(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  std::span<real> PartialQuadratureFunctionView::getIntegrationPointValues(
      const size_type e, const size_type i) {
    return this->getIntegrationPointValues(this->qspace->getOffset(e) + i);
  }  // end of getIntegrationPointValues

  PartialQuadratureFunction::~PartialQuadratureFunction() = default;

  /*!
   * \brief base class of the coefficients used to build a grid function from
   * partial quadrature functions.
   *
   * The values of a partial quadrature function are only known at the
   * integration points of its quadrature space. In each element, those values
   * are projected (in the L2 sense) on the shape functions of the element. The
   * coefficients return the value of this local projection, which can be
   * evaluated anywhere in the element and in particular at its nodes.
   *
   * \note The `index` member of the integration point passed to `Eval` is
   * deliberately not used. Its meaning depends on the caller: it is the index
   * of a node when a grid function is projected, the index of a quadrature
   * point in an integrator.
   */
  struct PartialQuadratureFunctionsCoefficientBase {
    /*!
     * \param[in] s: finite element space of the grid function
     * \param[in] fcts: functions
     */
    PartialQuadratureFunctionsCoefficientBase(
        const mfem::FiniteElementSpace& s,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : fespace(&s) {
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
    //
    static void doScalarFunctionsChecks(
        attributes::Throwing,
        const std::unordered_map<size_type,
                                 ImmutablePartialQuadratureFunctionView>&
            fcts) {
      for (const auto& [mid, f] : fcts) {
        static_cast<void>(mid);
        if (f.getNumberOfComponents() != 1) {
          raise("non scalar function given");
        }
      }
    }  // end of checkScalarFunction

    /*!
     * \brief evaluate the local projection of a function
     * \param[out] values: values of the projection at the given point
     * \param[in] f: function
     * \param[in] tr: transformation of the element of the grid function
     * \param[in] ip: point at which the projection is evaluated
     * \param[in] n: number of the element in the mesh of the quadrature space
     */
    void evaluate(mfem::Vector& values,
                  const ImmutablePartialQuadratureFunctionView& f,
                  mfem::ElementTransformation& tr,
                  const mfem::IntegrationPoint& ip,
                  const size_type n) {
      const auto& fe = *(this->fespace->GetFE(tr.ElementNo));
      if (this->element != tr.ElementNo) {
        this->computeNodalValues(f, fe, tr, n);
        this->element = tr.ElementNo;
        tr.SetIntPoint(&ip);
      }
      fe.CalcShape(ip, this->shape);
      this->nodal_values.MultTranspose(this->shape, values);
    }  // end of evaluate

    std::unordered_map<size_type, ImmutablePartialQuadratureFunctionView>
        functions;

   private:
    /*!
     * \brief compute the values at the nodes of an element of the local
     * projection of a function
     * \param[in] f: function
     * \param[in] fe: finite element of the grid function
     * \param[in] tr: transformation of the element of the grid function
     * \param[in] n: number of the element in the mesh of the quadrature space
     */
    void computeNodalValues(const ImmutablePartialQuadratureFunctionView& f,
                            const mfem::FiniteElement& fe,
                            mfem::ElementTransformation& tr,
                            const size_type n) {
      const auto nnodes = fe.GetDof();
      const auto nc = f.getNumberOfComponents();
      this->shape.SetSize(nnodes);
      // mass matrix, integrated exactly
      this->mass_matrix.SetSize(nnodes, nnodes);
      this->mass_matrix = 0.;
      const auto& mir =
          mfem::IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + tr.OrderW());
      for (int i = 0; i != mir.GetNPoints(); ++i) {
        const auto& ip = mir.IntPoint(i);
        tr.SetIntPoint(&ip);
        fe.CalcShape(ip, this->shape);
        mfem::AddMult_a_VVt(ip.weight * tr.Weight(), this->shape,
                            this->mass_matrix);
      }
      // right hand side, integrated with the rule of the quadrature space
      this->rhs.SetSize(nnodes, nc);
      this->rhs = 0.;
      const auto& ir = f.getPartialQuadratureSpace().getIntegrationRule(fe, tr);
      for (int i = 0; i != ir.GetNPoints(); ++i) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        fe.CalcShape(ip, this->shape);
        const auto w = ip.weight * tr.Weight();
        const auto fvalues = f.getIntegrationPointValues(n, i);
        for (int k = 0; k != nnodes; ++k) {
          for (size_type c = 0; c != nc; ++c) {
            this->rhs(k, c) += w * this->shape[k] * fvalues[c];
          }
        }
      }
      //
      this->mass_matrix.Invert();
      this->nodal_values.SetSize(nnodes, nc);
      mfem::Mult(this->mass_matrix, this->rhs, this->nodal_values);
    }  // end of computeNodalValues

    //! \brief finite element space of the grid function
    const mfem::FiniteElementSpace* fespace;
    //! \brief element for which the nodal values have been computed
    int element = -1;
    //! \brief values of the local projection at the nodes of the element
    mfem::DenseMatrix nodal_values;
    //! \brief mass matrix of the element
    mfem::DenseMatrix mass_matrix;
    //! \brief right hand side of the local projection
    mfem::DenseMatrix rhs;
    //! \brief values of the shape functions
    mfem::Vector shape;
  };  // end of PartialQuadratureFunctionsCoefficientBase

  struct PartialQuadratureFunctionsScalarCoefficient final
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::Coefficient {
    //
    PartialQuadratureFunctionsScalarCoefficient(
        const mfem::FiniteElementSpace& s,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(s, fcts), value(1) {
      doScalarFunctionsChecks(throwing, this->functions);
    }
    //
    PartialQuadratureFunctionsScalarCoefficient(
        PartialQuadratureFunctionsScalarCoefficient&&) = delete;
    PartialQuadratureFunctionsScalarCoefficient(
        const PartialQuadratureFunctionsScalarCoefficient&) = default;
    PartialQuadratureFunctionsScalarCoefficient& operator=(
        PartialQuadratureFunctionsScalarCoefficient&&) = delete;
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
      this->evaluate(this->value, p->second, tr, i, tr.ElementNo);
      return this->value[0];
    }  // end of Eval

   private:
    //! \brief value of the local projection
    mfem::Vector value;
  };

  struct PartialQuadratureFunctionsScalarCoefficientII final
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::Coefficient {
    //
    PartialQuadratureFunctionsScalarCoefficientII(
        const mfem::FiniteElementSpace& s,
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(s, fcts),
          elts_mapping(m),
          value(1) {
      doScalarFunctionsChecks(throwing, this->functions);
    }
    //
    PartialQuadratureFunctionsScalarCoefficientII(
        PartialQuadratureFunctionsScalarCoefficientII&&) = delete;
    PartialQuadratureFunctionsScalarCoefficientII(
        const PartialQuadratureFunctionsScalarCoefficientII&) = default;
    PartialQuadratureFunctionsScalarCoefficientII& operator=(
        PartialQuadratureFunctionsScalarCoefficientII&&) = delete;
    PartialQuadratureFunctionsScalarCoefficientII& operator=(
        const PartialQuadratureFunctionsScalarCoefficientII&) = delete;
    //
    double Eval(mfem::ElementTransformation& tr,
                const mfem::IntegrationPoint& i) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        return 0.;
      }
      const auto n = this->elts_mapping[tr.ElementNo];
      this->evaluate(this->value, p->second, tr, i, n);
      return this->value[0];
    }  // end of Eval
   private:
    const mfem::Array<int>& elts_mapping;
    //! \brief value of the local projection
    mfem::Vector value;
  };

  struct PartialQuadratureFunctionsVectorCoefficient final
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::VectorCoefficient {
    //
    PartialQuadratureFunctionsVectorCoefficient(
        const mfem::FiniteElementSpace& s,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(s, fcts),
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
        values.SetSize(this->GetVDim());
        this->evaluate(values, p->second, tr, ip, tr.ElementNo);
      }
    }  // end of Eval
  };

  struct PartialQuadratureFunctionsVectorCoefficientII final
      : public PartialQuadratureFunctionsCoefficientBase,
        public mfem::VectorCoefficient {
    //
    PartialQuadratureFunctionsVectorCoefficientII(
        const mfem::FiniteElementSpace& s,
        const mfem::Array<int>& m,
        const std::vector<ImmutablePartialQuadratureFunctionView>& fcts)
        : PartialQuadratureFunctionsCoefficientBase(s, fcts),
          mfem::VectorCoefficient(fcts.at(0).getNumberOfComponents()),
          elts_mapping(m) {}
    //
    PartialQuadratureFunctionsVectorCoefficientII(
        PartialQuadratureFunctionsVectorCoefficientII&&) = delete;
    PartialQuadratureFunctionsVectorCoefficientII(
        const PartialQuadratureFunctionsVectorCoefficientII&) = default;
    PartialQuadratureFunctionsVectorCoefficientII& operator=(
        PartialQuadratureFunctionsVectorCoefficientII&&) = delete;
    PartialQuadratureFunctionsVectorCoefficientII& operator=(
        const PartialQuadratureFunctionsVectorCoefficientII&) = delete;
    //
    void Eval(mfem::Vector& values,
              mfem::ElementTransformation& tr,
              const mfem::IntegrationPoint& ip) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        values = 0.;
      } else {
        const auto n = this->elts_mapping[tr.ElementNo];
        values.SetSize(this->GetVDim());
        this->evaluate(values, p->second, tr, ip, n);
      }
    }  // end of Eval

   private:
    const mfem::Array<int>& elts_mapping;
  };

  template <bool parallel>
  [[nodiscard]] static bool update_impl(
      Context& ctx,
      PartialQuadratureFunctionView& dest,
      const GridFunction<parallel>& src) noexcept {
    const auto& qspace = dest.getPartialQuadratureSpace();
    const auto& fed = qspace.getFiniteElementDiscretization();
    if constexpr (parallel) {
      if (!fed.describesAParallelComputation()) {
        return ctx.registerErrorMessage(
            "partial quadrature function is not built on a parallel finite "
            "element space");
      }
    } else {
      if (fed.describesAParallelComputation()) {
        return ctx.registerErrorMessage(
            "partial quadrature function is built on a parallel finite "
            "element space");
      }
    }
    const auto ofespace = qspace.template getFiniteElementSpace<parallel>(ctx);
    if (isInvalid(ofespace)) {
      return false;
    }
    const auto& fespace = *ofespace;
    if constexpr (parallel) {
      if (!fed.template isSlibing<parallel>(*(src.ParFESpace()))) {
        return ctx.registerErrorMessage("unmatched finite element space");
      }
    } else {
      if (!fed.template isSlibing<parallel>(*(src.FESpace()))) {
        return ctx.registerErrorMessage("unmatched finite element space");
      }
    }
    if (src.FESpace()->GetMesh() != fespace.GetMesh()) {
      return ctx.registerErrorMessage(
          "the grid function is not defined on the mesh of the partial "
          "quadrature space");
    }
    if (dest.getNumberOfComponents() != src.VectorDim()) {
      return ctx.registerErrorMessage(
          "unmatched number of components (" +
          std::to_string(dest.getNumberOfComponents()) +
          " for the quadrature function and " +
          std::to_string(src.VectorDim()) + " for the grid function)");
    }
    if (dest.getNumberOfComponents() == 1) {
      // scalar case
      for (const auto& [i, offset] : qspace.getOffsets()) {
        static_cast<void>(offset);
        const auto& fe = *(fespace.GetFE(i));
        auto& tr = *(fespace.GetElementTransformation(i));
        const auto& ir = qspace.getIntegrationRule(fe, tr);
        for (size_type g = 0; g != ir.GetNPoints(); ++g) {
          // get the gradients of the shape functions
          const auto& ip = ir.IntPoint(g);
          tr.SetIntPoint(&ip);
          dest.getIntegrationPointValue(i, g) = src.GetValue(tr, ip);
        }
      }
    } else {
      // vectorial case
      const auto nc = dest.getNumberOfComponents();
      for (const auto& [i, offset] : qspace.getOffsets()) {
        static_cast<void>(offset);
        const auto& fe = *(fespace.GetFE(i));
        auto& tr = *(fespace.GetElementTransformation(i));
        const auto& ir = qspace.getIntegrationRule(fe, tr);
        for (size_type g = 0; g != ir.GetNPoints(); ++g) {
          // get the gradients of the shape functions
          const auto& ip = ir.IntPoint(g);
          tr.SetIntPoint(&ip);
          auto v = dest.getIntegrationPointValues(i, g);
          auto tmp = mfem::Vector(v.data(), nc);
          src.GetVectorValue(tr, ip, tmp);
        }
      }
    }
    return true;
  }  // end of update_impl

#ifdef MFEM_USE_MPI
  bool update(Context& ctx,
              PartialQuadratureFunctionView& dest,
              const GridFunction<true>& src) noexcept {
    return update_impl<true>(ctx, dest, src);
  }    // end of update
#endif /* MFEM_USE_MPI */

  bool update(Context& ctx,
              PartialQuadratureFunctionView& dest,
              const GridFunction<false>& src) noexcept {
    return update_impl<false>(ctx, dest, src);
  }  // end of update

  /*!
   * \brief check that the given functions are defined on all the materials or
   * on all the boundaries of the mesh of a grid function, and on those only
   * \param[in, out] ctx: execution context
   * \param[in] fcts: functions
   * \param[in] mesh: mesh of the grid function
   *
   * \note the elements on which no function is defined would contribute to
   * the nodal averages with null values.
   * \note location identifiers are compared, rather than attributes, as the
   * attributes of a submesh defined on boundaries are boundary identifiers.
   * \note a partial quadrature space defined on a boundary is built on the
   * submesh of this boundary, which must be the mesh of the grid function.
   */
  template <bool parallel>
  [[nodiscard]] static bool checkGridFunctionMesh(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<parallel>& mesh) noexcept {
    // locations on which the functions are defined
    auto locations = std::set<LocationIdentifier>{};
    for (const auto& f : fcts) {
      const auto& qspace = f.getPartialQuadratureSpace();
      const auto& fed = qspace.getFiniteElementDiscretization();
      const auto omesh = qspace.template getMesh<parallel>(ctx);
      if (isInvalid(omesh)) {
        return false;
      }
      const auto ol = fed.getLocationIdentifier(ctx, *omesh, qspace.getId());
      if (isInvalid(ol)) {
        return false;
      }
      if (isValid(ol->boundary_identifier) && (&(*omesh) != &mesh)) {
        return ctx.registerErrorMessage(
            "the function defined on '" + qspace.getLocationName() +
            "' can only be projected on the submesh of this boundary");
      }
      if (!locations.empty()) {
        const auto& l = *(locations.begin());
        if (isValid(l.material_identifier) !=
            isValid(ol->material_identifier)) {
          return ctx.registerErrorMessage(
              "functions defined on materials and on boundaries can not be "
              "mixed");
        }
      }
      if (!locations.insert(*ol).second) {
        return ctx.registerErrorMessage("multiple functions defined on '" +
                                        qspace.getLocationName() + "'");
      }
    }
    // locations associated with the mesh of the grid function
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    auto mesh_locations = std::set<LocationIdentifier>{};
    for (const auto a : mesh.attributes) {
      const auto ol = fed.getLocationIdentifier(ctx, mesh, a);
      if (isInvalid(ol)) {
        return false;
      }
      mesh_locations.insert(*ol);
    }
    if (locations != mesh_locations) {
      return ctx.registerErrorMessage(
          "the materials or boundaries on which the functions are defined do "
          "not match those of the mesh of the grid function");
    }
    return true;
  }  // end of checkGridFunctionMesh

  template <bool parallel>
  static std::unique_ptr<GridFunction<parallel>> makeGridFunction_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    if (fcts.empty()) {
      return ctx.registerErrorMessage("no functions defined");
    }
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, fed.getMesh<parallel>())) {
      return {};
    }
    auto m = fed.getFiniteElementSpacesManager();
    auto fespace = m.getFiniteElementSpace<parallel>(ctx, n);
    if (isInvalid(fespace)) {
      return {};
    }
    return std::make_unique<GridFunction<parallel>>(fespace.get());
  }

  template <>
  std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
#ifdef MFEM_USE_MPI
    return makeGridFunction_impl<true>(ctx, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    return makeGridFunction_impl<false>(ctx, fcts);
  }

  template <bool parallel>
  static std::unique_ptr<GridFunction<parallel>> makeGridFunction_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<parallel>& mesh) {
    if (fcts.empty()) {
      return ctx.registerErrorMessage("no functions defined");
    }
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, mesh)) {
      return {};
    }
    auto fespace = fed.getFiniteElementSpacesManager()
                       .template getFiniteElementSpace<parallel>(ctx, mesh, n);
    if (isInvalid(fespace)) {
      return {};
    }
    return std::make_unique<GridFunction<parallel>>(fespace.get());
  }

  template <>
  std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<true>& mesh) {
#ifdef MFEM_USE_MPI
    return makeGridFunction_impl<true>(ctx, fcts, mesh);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<false>& mesh) {
    return makeGridFunction_impl<false>(ctx, fcts, mesh);
  }

  template <bool parallel>
  static std::unique_ptr<GridFunction<parallel>> makeGridFunction_impl(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<parallel>& mesh) {
    if (fcts.empty()) {
      return ctx.registerErrorMessage("no functions defined");
    }
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, mesh)) {
      return {};
    }
    auto fespace = fed.getFiniteElementSpacesManager()
                       .template getFiniteElementSpace<parallel>(ctx, mesh, n);
    if (isInvalid(fespace)) {
      return {};
    }
    return std::make_unique<GridFunction<parallel>>(fespace.get());
  }

  template <>
  std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<true>& mesh) {
#ifdef MFEM_USE_MPI
    return makeGridFunction_impl<true>(ctx, fcts, mesh);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<false>& mesh) {
    return makeGridFunction_impl<false>(ctx, fcts, mesh);
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
    auto ctx = Context{};
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, mesh)) {
      raise(ctx.getErrorMessage());
    }
    if (n == 1u) {
      auto c = PartialQuadratureFunctionsScalarCoefficient(*fespace, fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c = PartialQuadratureFunctionsVectorCoefficient(*fespace, fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
#ifdef MFEM_USE_MPI
    updateGridFunction_impl<true>(f, fcts);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
    updateGridFunction_impl<false>(f, fcts);
  }

  template <bool parallel>
  static void updateGridFunction_impl(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<parallel>& mesh) {
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& fes = fed.getFiniteElementSpace<parallel>();
    const auto& fespace = f.FESpace();
    if ((fespace->GetMesh() != &mesh) ||  //
        (fespace->GetVDim() != n) ||      //
        (fes.FEColl() != fespace->FEColl()) ||
        (fes.GetOrdering() != fespace->GetOrdering())) {
      raise("inconsistent grid function");
    }
    auto ctx = Context{};
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, mesh)) {
      raise(ctx.getErrorMessage());
    }
    if (n == 1u) {
      auto c = PartialQuadratureFunctionsScalarCoefficient(*fespace, fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c = PartialQuadratureFunctionsVectorCoefficient(*fespace, fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<true>& mesh) {
#ifdef MFEM_USE_MPI
    updateGridFunction_impl<true>(f, fcts, mesh);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const Mesh<false>& mesh) {
    updateGridFunction_impl<false>(f, fcts, mesh);
  }

  template <bool parallel>
  static void updateGridFunction_impl(
      GridFunction<parallel>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<parallel>& mesh) {
    const auto n = fcts.at(0).getNumberOfComponents();
    const auto& fed =
        fcts.at(0).getPartialQuadratureSpace().getFiniteElementDiscretization();
    const auto& fes = fed.getFiniteElementSpace<parallel>();
    const auto& fespace = f.FESpace();
    if ((fespace->GetMesh() != &mesh) ||  //
        (fespace->GetVDim() != n) ||      //
        (fes.FEColl() != fespace->FEColl()) ||
        (fes.GetOrdering() != fespace->GetOrdering())) {
      raise("inconsistent grid function");
    }
    auto ctx = Context{};
    if (!checkGridFunctionMesh<parallel>(ctx, fcts, mesh)) {
      raise(ctx.getErrorMessage());
    }
    const auto on_boundaries = fed.isDefinedOnBoundaries(ctx, mesh);
    if (isInvalid(on_boundaries)) {
      raise(ctx.getErrorMessage());
    }
    if (*on_boundaries) {
      // the partial quadrature spaces are built on the submesh
      if (n == 1u) {
        auto c = PartialQuadratureFunctionsScalarCoefficient(*fespace, fcts);
        f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
      } else {
        auto c = PartialQuadratureFunctionsVectorCoefficient(*fespace, fcts);
        f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
      }
      return;
    }
    if (n == 1u) {
      auto c = PartialQuadratureFunctionsScalarCoefficientII(
          *fespace, mesh.GetParentElementIDMap(), fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c = PartialQuadratureFunctionsVectorCoefficientII(
          *fespace, mesh.GetParentElementIDMap(), fcts);
      f.ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<true>(
      GridFunction<true>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<true>& mesh) {
#ifdef MFEM_USE_MPI
    updateGridFunction_impl<true>(f, fcts, mesh);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }

  template <>
  MFEM_MGIS_EXPORT void updateGridFunction<false>(
      GridFunction<false>& f,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const SubMesh<false>& mesh) {
    updateGridFunction_impl<false>(f, fcts, mesh);
  }

}  // end of namespace mfem_mgis

#ifdef MGIS_FUNCTION_SUPPORT

namespace mfem_mgis {

  const PartialQuadratureSpace& getSpace(
      const ImmutablePartialQuadratureFunctionView& f) {
    return f.getPartialQuadratureSpace();
  }

  const PartialQuadratureSpace& getSpace(const PartialQuadratureFunction& f) {
    return f.getPartialQuadratureSpace();
  }

}  // end of namespace mfem_mgis

#endif /* MGIS_FUNCTION_SUPPORT */
