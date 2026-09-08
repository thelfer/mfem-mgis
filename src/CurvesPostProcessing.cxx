/*!
 * \file   src/CurvesPostProcessing.cxx
 * \brief  This file implements the `CurvesPostProcessing` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#include <locale>
#include <iterator>
#include <algorithm>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/PostProcessing/AbstractCurve.hxx"
#include "MFEMMGIS/PostProcessing/CurvesPostProcessing.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  CurvesPostProcessing::getParametersDescription() noexcept {
    auto d = PostProcessingBase::getParametersDescription();
    const auto cd = CurvesWriter::getParametersDescription();
    d.insert(cd.begin(), cd.end());
    d.insert({{"ExecuteInitialPostProcessing",
               "export values at the beginning of the first time step"}});
    return d;
  }  // end of getParametersDescription

  std::string CurvesPostProcessing::getDescription() noexcept {
    return "a post-processing meant to export values extracted "
           "using instances of the AbstractCurve class to a file. "
           "A list of available curves can be retrieved from the "
           "CurveFactory class";
  }  // end of CurvesPostProcessing

  CurvesPostProcessing::CurvesPostProcessing(Context &ctx,
                                             PhysicalSystem &ps,
                                             const Parameters &params)
      : PostProcessingBase(ps, params, true),
        writer(ctx,
               ps,
               extract(
                   throwing, params, CurvesWriter::getParametersDescription())),
        executeInitialPostProcessing(get_if<bool>(
            throwing, params, "ExecuteInitialPostProcessing", true)) {
    checkParameters(throwing, params,
                    CurvesPostProcessing::getParametersDescription());
  }  // end of CurvesPostProcessing

  std::string CurvesPostProcessing::getName() const noexcept {
    return "Curves";
  }  // end of getName

  bool CurvesPostProcessing::executeInitialPostProcessingTasks(
      Context &ctx, const real t) noexcept {
    if (!this->writer.writeFileHeader(ctx)) {
      return false;
    }
    return this->writer.writeValues(ctx, {.begin = t, .end = t, .dt = 0}, bts);
  }  // end of executeInitialPostProcessingTasks

  bool CurvesPostProcessing::executePostProcessingTasks(
      Context &ctx,
      const TimeStep &ts,
      const bool isPostProcessingRequired) noexcept {
    if ((this->allTimeSteps) || isPostProcessingRequired) {
      return this->writer.writeValues(ctx, ts, ets);
    }
    return true;
  }  // end of executePostProcessingTasks

#pragma message("HERE")
  //   bool CurvesPostProcessing::add(Context &ctx,
  //                                  std::string_view n,
  //                                  const Parameters &params) noexcept {
  //     auto c = this->physicalSystem.addCurve(ctx, n, params);
  //     if (c.get() == nullptr) {
  //       return ctx.registerErrorMessage("invalid curve");
  //     }
  //     this->writer.push_back(c);
  //     return true;
  //   }  // end of add

  bool CurvesPostProcessing::add(Context &ctx,
                                 std::shared_ptr<AbstractCurve> c) noexcept {
    if (c.get() == nullptr) {
      return ctx.registerErrorMessage("invalid curve");
    }
#pragma message("HERE")
    //     if (!this->physicalSystem->addCurve(ctx, c)) {
    //       return false;
    //     }
    return this->writer.addCurve(ctx, c);
  }  // end of add

  CurvesPostProcessing::~CurvesPostProcessing() noexcept = default;

}  // end of namespace mfem_mgis