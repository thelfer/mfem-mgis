/*!
 * \file   MFEMMGIS/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/04/2025
 */

#ifndef LIB_MFEM_MGIS_ALGORITHMS_HXX
#define LIB_MFEM_MGIS_ALGORITHMS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionEvaluator.hxx"

namespace mfem_mgis {

  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <size_type N, typename PartialQuadratureFunctionEvaluatorType>
  bool assign(PartialQuadratureFunction&,
              PartialQuadratureFunctionEvaluatorType) requires(N > 0);
  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename PartialQuadratureFunctionEvaluatorType>
  bool assign(PartialQuadratureFunction&,
              PartialQuadratureFunctionEvaluatorType);

  /*!
  template <typename ValueType, typename BinaryOperator>
  ValueType reduce(ImmutablePartialQuadratureFunctionView f,
                   const ValueType init,
                   BinaryOperator op) {
    constexpr bool expects_scalar_function =
        requires(const real v, const ValueType& v2) {
      op(v, v2);
    };
    const auto ne =
        f.getPartialQuadratureSpace().getNumberOfIntegrationPoints();
    auto r = init;
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (expects_scalar_function) {
        const auto v = f.getIntegrationPointValue(i);
        r = op(v, r);
      } else {
        const auto v = f.getIntegrationPointValues(i);
        r = op(v, r);
      }
    }
    return r;
  }
   */

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Algorithms.ixx"

#endif /* LIB_MFEM_MGIS_ALGORITHMS_HXX */
