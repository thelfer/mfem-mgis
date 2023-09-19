/*!
 * \file   src/BehaviourIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"

namespace mfem_mgis {

  BehaviourIntegrator::~BehaviourIntegrator() = default;

  template <bool parallel, typename Functor>
  static void performLoopOverIntegrationPoints(Functor& f,
                                               const BehaviourIntegrator& bi) {
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto& fespace = qspace.getFiniteElementDiscretization()
                              .getFiniteElementSpace<parallel>();
    for (const auto [e, eo] : qspace.getOffsets()) {
      const auto& fe = *(fespace.GetFE(e));
      auto& tr = *(fespace.GetElementTransformation(e));
      const auto& ir = bi.getIntegrationRule(fe, tr);
      for (size_type i = 0; i != ir.GetNPoints(); ++i) {
        const auto& ip = ir.IntPoint(i);
        tr.SetIntPoint(&ip);
        f(fe, tr, ip, eo + i);
      }
    }
  }  // end of performLoopOverIntegrationPoints

  template <bool parallel>
  static real BehaviourIntegrator_computeMeasure(
      const BehaviourIntegrator& bi) {
    auto r = real{};
    auto integrate =
        [&r, &bi](const mfem::FiniteElement&, mfem::ElementTransformation& tr,
                  const mfem::IntegrationPoint& ip, const size_type) {
          r += bi.getIntegrationPointWeight(tr, ip);
        };
    performLoopOverIntegrationPoints<parallel>(integrate, bi);
    return r;
  }  // end of BehaviourIntegrator_computeMeasure

  real computeMeasure(const BehaviourIntegrator& bi) {
    const auto& qspace = bi.getPartialQuadratureSpace();
    const auto& fed = qspace.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto lv = BehaviourIntegrator_computeMeasure<true>(bi);
      auto r = real{};
      MPI_Reduce(&lv, &r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return BehaviourIntegrator_computeMeasure<false>(bi);
  }  // end of computeIntegral

  template <bool parallel>
  static real BehaviourIntegrator_computeScalarIntegral(
      const BehaviourIntegrator& bi,
      const ImmutablePartialQuadratureFunctionView& f) {
    const auto* const values = f.getValues().data() + f.getDataOffset();
    auto r = real{};
    if (f.getDataStride() == 1) {
      auto integrate = [&r, &bi, values](const mfem::FiniteElement&,
                                         mfem::ElementTransformation& tr,
                                         const mfem::IntegrationPoint& ip,
                                         const size_type o) {
        const auto v = values[o];
        r += v * bi.getIntegrationPointWeight(tr, ip);
      };
      performLoopOverIntegrationPoints<parallel>(integrate, bi);
    } else {
      const auto s = f.getDataStride();
      auto integrate = [&r, &bi, &s, values](const mfem::FiniteElement&,
                                             mfem::ElementTransformation& tr,
                                             const mfem::IntegrationPoint& ip,
                                             const size_type o) {
        const auto v = values[o * s];
        r += v * bi.getIntegrationPointWeight(tr, ip);
      };
      performLoopOverIntegrationPoints<parallel>(integrate, bi);
    }
    return r;
  }  // end of BehaviourIntegrator_computeScalarIntegral

  template <>
  real computeIntegral(const BehaviourIntegrator& bi,
                       const ImmutablePartialQuadratureFunctionView& f) {
    const auto& qspace = f.getPartialQuadratureSpace();
    if (&qspace != &(bi.getPartialQuadratureSpace())) {
      raise(
          "computeIntegral: the partial quadrature space of the behaviour "
          "integrator don't match the one of the partial quadrature function");
    }
    if (f.getNumberOfComponents() != 1u) {
      raise(
          "computeIntegral: the partial quadrature function is not scalar "
          "valuated");
    }
    const auto& fed = qspace.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto lv = BehaviourIntegrator_computeScalarIntegral<true>(bi, f);
      auto r = real{};
      MPI_Reduce(&lv, &r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return r;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return BehaviourIntegrator_computeScalarIntegral<false>(bi, f);
  }  // end of computeIntegral

}  // end of namespace mfem_mgis
