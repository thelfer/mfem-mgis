/*!
 * \file   src/PointsSetCurvesPostProcessing.cxx
 * \brief  This file implements the `PointsSetCurvesPostProcessing` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#include <locale>
#include <iterator>
#include <algorithm>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/PostProcessing/AbstractCurve.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurvesPostProcessing.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  PointsSetCurvesPostProcessing::getParametersDescription() noexcept {
    auto d = PostProcessingBase::getParametersDescription();
    const auto cd = PointsSetCurvesWriter::getParametersDescription();
    d.insert(cd.begin(), cd.end());
    d.insert({{"ExecuteInitialPostProcessing",
               "export values at the beginning of the first time step"}});
    return d;
  }  // end of getParametersDescription

  std::string PointsSetCurvesPostProcessing::getDescription() noexcept {
    return "a post-processing meant to write values along a points set to a "
           "file";
  }  // end of PointsSetCurvesPostProcessing

  PointsSetCurvesPostProcessing::PointsSetCurvesPostProcessing(
      PhysicalSystem &ps, const Parameters &params)
      : PostProcessingBase(ps, params, true),
        writer(ps,
               extract(throwing,
                       params,
                       PointsSetCurvesWriter::getParametersDescription())),
        executeInitialPostProcessing(get_if<bool>(
            throwing, params, "ExecuteInitialPostProcessing", true)) {
    checkParameters(throwing, params,
                    PointsSetCurvesPostProcessing::getParametersDescription());
  }  // end of PointsSetCurvesPostProcessing

  std::string PointsSetCurvesPostProcessing::getName() const noexcept {
    return "PointsSetCurves";
  }  // end of getName

#ifdef MFEM_USE_MPI

  bool PointsSetCurvesPostProcessing::add(
      Context &ctx, std::string_view n, const GridFunction<true> &f) noexcept {
    return this->writer.add(ctx, n, f);
  }  // end of add

#endif /* MFEM_USE_MPI */

  bool PointsSetCurvesPostProcessing::add(
      Context &ctx, std::string_view n, const GridFunction<false> &f) noexcept {
    return this->writer.add(ctx, n, f);
  }  // end of add

  bool PointsSetCurvesPostProcessing::executeInitialPostProcessingTasks(
      Context &ctx, const real t) noexcept {
    if (!this->writer.writeFileHeader(ctx)) {
      return false;
    }
    return this->writer.writeValues(ctx, {.begin = t, .end = t, .dt = 0}, bts);
  }  // end of executeInitialPostProcessingTasks

  bool PointsSetCurvesPostProcessing::executePostProcessingTasks(
      Context &ctx,
      const TimeStep &ts,
      const bool isPostProcessingRequired) noexcept {
    if ((this->allTimeSteps) || isPostProcessingRequired) {
      return this->writer.writeValues(ctx, ts, ets);
    }
    return true;
  }  // end of executePostProcessingTasks

  PointsSetCurvesPostProcessing::~PointsSetCurvesPostProcessing() noexcept =
      default;

}  // end of namespace mfem_mgis