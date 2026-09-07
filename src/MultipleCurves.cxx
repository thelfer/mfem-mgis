/*!
 * \file   MultipleCurves.cxx
 * \brief  This file implements the MultipleCurves class
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#include "MFEMMGIS/PostProcessing/MultipleCurves.hxx"

namespace mfem_mgis {

  MultipleCurves::MultipleCurves() = default;

  MultipleCurves::MultipleCurves(Context &,
                                 PhysicalSystem &,
                                 const Parameters &) {
  }  // end of MultipleCurves

  bool MultipleCurves::addCurve(Context &ctx,
                                std::shared_ptr<const AbstractCurve> c) {
    if (c.get() == nullptr) {
      return ctx.registerErrorMessage("invalid curve");
    }
    this->curves.push_back(std::move(c));
    return true;
  }  // end of addCurve

  std::vector<std::string> MultipleCurves::getDescriptions() const noexcept {
    auto r = std::vector<std::string>{};
    for (const auto &c : this->curves) {
      const auto d = c->getDescriptions();
      r.insert(r.end(), d.begin(), d.end());
    }
    return r;
  }  // end of getDescriptions

  std::optional<std::vector<real>> MultipleCurves::getValues(
      Context &ctx, const TimeStepStage ts) const noexcept {
    auto r = std::vector<real>{};
    for (const auto &c : this->curves) {
      const auto ovalues = c->getValues(ctx, ts);
      if (isInvalid(ovalues)) {
        return {};
      }
      r.insert(r.end(), ovalues->begin(), ovalues->end());
    }
    return r;
  }  // end of getValues

  MultipleCurves::~MultipleCurves() noexcept = default;

}  // end of namespace mfem_mgis
