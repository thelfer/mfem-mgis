/*!
 * \file   MFEMMGIS/Algorithm.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MFEMMGIS_ALGORITHMS_IXX
#define LIB_MFEMMGIS_ALGORITHMS_IXX

namespace mfem_mgis {

  template <size_type N, PartialQuadratureFunctionEvalutorConcept EvaluatorType>
  bool assign(PartialQuadratureFunction& f, EvaluatorType e) requires(N > 0) {
    checkMatchingQuadratureSpaces(f, e);
    raise_if(f.getNumberOfComponents() != N,
             "assign: invalid number of components for the left hand size");
    raise_if(e.getNumberOfComponents() != N,
             "assign: invalid number of components for the right hand size");
    //
    e.check();
    e.allocateWorkspace();
    //
    const auto qspace = f.getPartialQuadratureSpace();
    const auto ne = qspace.getNumberOfIntegrationPoints();
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (N == 1) {
        auto& v = f.getIntegrationPointValue(i);
        v = e(i);
      } else {
        auto lhs_values = f.template getIntegrationPointValues<N>(i);
        const auto& rhs_values = e(i);
        algorithm::copy<N>(rhs_values.begin(), rhs_values.end(),
                           lhs_values.begin());
      }
    }
    return true;
  }  // end of assign

  template <PartialQuadratureFunctionEvalutorConcept EvaluatorType>
  bool assign(PartialQuadratureFunction& f,
              EvaluatorType e) {
    raise_if(&f.getPartialQuadratureSpace() != &e.getPartialQuadratureSpace(),
             "assign: unmatched number of components for the left hand size "
             "and the right hand side");
    raise_if(f.getNumberOfComponents() != e.getNumberOfComponents(),
             "assign: unmatched number of components for the left hand size "
             "and the right hand side");
    //
    e.check();
    e.allocateWorkspace();
    //
    const auto qspace = f.getPartialQuadratureSpace();
    const auto ne = qspace.getNumberOfIntegrationPoints();
    if (f.isScalar()) {
      using result_type = std::invoke_result_t<EvaluatorType, size_type>;
      if constexpr (std::same_as<std::decay_t<result_type>, real>) {
        for (size_type i = 0; i != ne; ++i) {
          auto& lhs_value = f.getIntegrationPointValue(i);
          lhs_value = e(i);
        }
      } else {
        for (size_type i = 0; i != ne; ++i) {
          auto lhs_value = f.getIntegrationPointValue(i);
          lhs_value = *(e(i).begin());
        }
      }
    } else {
      for (size_type i = 0; i != ne; ++i) {
        auto lhs_values = f.getIntegrationPointValues(i);
        const auto& rhs_values = e(i);
        std::copy(rhs_values.begin(), rhs_values.end(), lhs_values.begin());
      }
    }
    return true;
  }  // end of assign

} // end of mfem_mgis

#endif /* LIB_MFEMMGIS_ALGORITHMS_IXX */
