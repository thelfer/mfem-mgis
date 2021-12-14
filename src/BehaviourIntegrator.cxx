/*!
 * \file   src/BehaviourIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#include <cmath>
#include <cstdlib>
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
    const auto m = qspace.getId();
    for (size_type i = 0; i != fespace.GetNBE(); ++i) {
      if (fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(fespace.GetFE(i));
      auto& tr = *(fespace.GetElementTransformation(i));
      const auto& ir = bi.getIntegrationRule(fe, tr);
      // element offset
      const auto eoffset = qspace.getOffset(i);
      for (size_type j = 0; j != ir.GetNPoints(); ++j) {
        const auto& ip = ir.IntPoint(j);
        tr.SetIntPoint(&ip);
        f(fe, tr, ip, eoffset + i);
      }
    }
  }  // end of performLoopOverIntegrationPoints

  template <bool parallel>
  real BehaviourIntegrator_computeScalarIntegral(
      const BehaviourIntegrator& bi,
      const ImmutablePartialQuadratureFunctionView& f) {
    auto r = real{};
    const auto* const values = f.getValues().data() + f.getInitialDataOffset();
    if (f.getDataStride() == 0) {
      auto integrate = [&bi, &r, values](const mfem::FiniteElement&,
                                         mfem::ElementTransformation& tr,
                                         const mfem::IntegrationPoint& ip,
                                         const size_type o) {
        const auto v = values[o];
        r += v * bi.getIntegrationPointWeight(tr, ip);
      };
      performLoopOverIntegrationPoints<parallel>(integrate, bi);
    } else {
      auto integrate = [&bi, &r, values](const mfem::FiniteElement&,
                                         mfem::ElementTransformation& tr,
                                         const mfem::IntegrationPoint& ip,
                                         const size_type o) {
        const auto v = values[o];
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
