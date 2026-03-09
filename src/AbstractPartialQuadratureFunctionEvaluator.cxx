/*!
 * \file   src/AbstractPartialQuadratureFunctionEvaluator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   07/03/2026
 */

#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/AbstractPartialQuadratureFunctionEvaluator.hxx"

namespace mfem_mgis {

  AbstractPartialQuadratureFunctionEvaluator::
      ~AbstractPartialQuadratureFunctionEvaluator() noexcept = default;

  std::shared_ptr<PartialQuadratureFunction> evaluate(
      Context& ctx,
      const AbstractPartialQuadratureFunctionEvaluator& e) noexcept {
    auto f = make_shared<PartialQuadratureFunction>(
        ctx, e.getPartialQuadratureSpacePointer(), e.getNumberOfComponents());
    if (f.get() == nullptr) {
      return {};
    }
    if (!e.evaluate(ctx, *f)) {
      return {};
    }
    return f;
  };

}  // end of namespace mfem_mgis
