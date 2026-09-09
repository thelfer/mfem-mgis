/*!
 * \file   src/PointsSetCurves.cxx
 * \brief  This file implements the `PointsSetCurves` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#include <cmath>
#include <locale>
#include <algorithm>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/GridFunctionUtilities.hxx"
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
#include "MFEMMGIS/GridFunctionInterpolator.hxx"
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
#include "MFEMMGIS/PostProcessing/PointsSetCurves.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  PointsSetCurves::getParametersDescription() noexcept {
    auto d = std::map<std::string, std::string>{};
    d.insert({{"PointsSet", "name or definition of the points set"},
              {"ExportCurvilinearAbscissa",
               "boolean stating if the curvilinear abscissa shall be exported"},
              {"ExportCoordinates",
               "boolean stating if the coordinates shall be exported"}});
    return d;
  }  // end of getParametersDescription

  PointsSetCurves::PointsSetCurves(const FiniteElementSpacesManager &manager,
                                   const Parameters &params)
      : fespaces_manager(manager),
        shallExportCurvilinearAbscissa(
            get_if<bool>(throwing, params, "ExportCurvilinearAbscissa", true)),
        shallExportCoordinates(
            get_if<bool>(throwing, params, "ExportCoordinates", true)) {
    checkParameters(throwing, params,
                    PointsSetCurves::getParametersDescription());
    if (!contains(params, "PointsSet")) {
      raise("no points set defined");
    }
    const auto m = this->fespaces_manager.getMeshDiscretization();
    const auto d = ::mfem_mgis::getSpaceDimension(m);
    if (d == 2) {
#if MGIS_HAVE_TFEL
      auto ctx = Context{};
      auto or_raise = ctx.getThrowingFailureHandler();
      this->points =
          makePointsSet<2>(ctx, m, get(throwing, params, "PointsSet")) |
          or_raise;
      if (std::get<std::vector<Point<2>>>(this->points).empty()) {
        raise("points set is empty");
      }
      this->curvilinearAbscissae = computeCurvilinearAbscissae(
          std::get<std::vector<Point<2>>>(this->points));
#endif /* MGIS_HAVE_TFEL */
    } else if (d == 3) {
#if MGIS_HAVE_TFEL
      auto ctx = Context{};
      auto or_raise = ctx.getThrowingFailureHandler();
      this->points =
          makePointsSet<3>(ctx, m, get(throwing, params, "PointsSet")) |
          or_raise;
      if (std::get<std::vector<Point<3>>>(this->points).empty()) {
        raise("points set is empty");
      }
      this->curvilinearAbscissae = computeCurvilinearAbscissae(
          std::get<std::vector<Point<3>>>(this->points));
#endif /* MGIS_HAVE_TFEL */
    } else {
      raise("PointsSetsCurve is only usable in 2D or 3D");
    }
  }  // end of PointsSetCurves

#ifdef MFEM_USE_MPI

  bool PointsSetCurves::add(Context &ctx,
                            std::string_view n,
                            const GridFunction<true> &f) noexcept {
    if (!this->fespaces_manager.manages(*(f.ParFESpace()))) {
      return ctx.registerErrorMessage(
          "the given grid function is defined on a finite element space which "
          "is not managed by the finite element spaces manager of which the "
          "points set curves is built");
    }
    for (const auto &[nf, vf] : this->gridfunctions) {
      static_cast<void>(vf);
      if (nf == n) {
        return ctx.registerErrorMessage("a grid function named '" +
                                        std::string{n} +
                                        "' has already been registred");
      }
    }
    const auto vf =
        std::variant<const GridFunction<true> *, const GridFunction<false> *>{
            &f};
    this->gridfunctions.push_back({std::string{n}, vf});
    return true;
  }  // end of add

#endif /* MFEM_USE_MPI */

  bool PointsSetCurves::add(Context &ctx,
                            std::string_view n,
                            const GridFunction<false> &f) noexcept {
    if (!this->fespaces_manager.manages(*(f.FESpace()))) {
      return ctx.registerErrorMessage(
          "the given grid function is defined on a finite element space which "
          "is not managed by the finite element spaces manager of which the "
          "points set curves is built");
    }
    for (const auto &[nf, vf] : this->gridfunctions) {
      static_cast<void>(vf);
      if (nf == n) {
        return ctx.registerErrorMessage("a grid function named '" +
                                        std::string{n} +
                                        "' has already been registred");
      }
    }
#ifdef MFEM_USE_MPI
    const auto vf =
        std::variant<const GridFunction<true> *, const GridFunction<false> *>{
            &f};
    this->gridfunctions.push_back({std::string{n}, vf});
#else  /* MFEM_USE_MPI */
    this->gridfunctions.push_back({std::string{n}, &f});
#endif /* MFEM_USE_MPI */
    return true;
  }

  bool PointsSetCurves::exportCurvilinearAbscissa() const noexcept {
    return this->shallExportCurvilinearAbscissa;
  }  // end of exportCurvilinearAbscissa

  OptionalReference<const std::vector<real>>
  PointsSetCurves::getCurvilinearAbscissa(Context &ctx) const noexcept {
    if (this->exportCurvilinearAbscissa()) {
      return {&(this->curvilinearAbscissae)};
    }
    return ctx.registerErrorMessage("curvilinear abscissa are not calculated");
  }  // end of getCurvilinearAbscissa

  bool PointsSetCurves::exportCoordinates() const noexcept {
    return this->shallExportCoordinates;
  }  // end of exportCoordinates

  size_type PointsSetCurves::getSpaceDimension() const noexcept {
    if (std::holds_alternative<std::vector<Point<2>>>(points)) {
      return 2;
    }
    return 3;
  }  // end of getSpaceDimension

  std::optional<std::vector<std::vector<real>>> PointsSetCurves::getCoordinates(
      Context &ctx) const noexcept {
    if (!this->exportCoordinates()) {
      return ctx.registerErrorMessage("coordinates are not exported");
    }
    auto r = std::vector<std::vector<real>>{};
#ifdef MGIS_HAVE_TFEL
    if (std::holds_alternative<std::vector<Point<2>>>(points)) {
      const auto &pts = std::get<std::vector<Point<2>>>(points);
      r.resize(2);
      auto &x = r[0];
      auto &y = r[1];
      x.reserve(pts.size());
      y.reserve(pts.size());
      for (std::size_t i = 0; i != pts.size(); ++i) {
        x[i] = pts[i][0];
        y[i] = pts[i][1];
      }
    } else {
      const auto &pts = std::get<std::vector<Point<3>>>(points);
      r.resize(3);
      auto &x = r[0];
      auto &y = r[1];
      auto &z = r[2];
      x.reserve(pts.size());
      y.reserve(pts.size());
      z.reserve(pts.size());
      for (std::size_t i = 0; i != pts.size(); ++i) {
        x[i] = pts[i][0];
        y[i] = pts[i][1];
        z[i] = pts[i][2];
      }
    }
#endif /* MGIS_HAVE_TFEL */
    return r;
  }  // end of getCoordinates

  std::vector<std::string> PointsSetCurves::getValuesDescription()
      const noexcept {
    auto d = std::vector<std::string>{};
    auto add = [&d](const size_type nc, const std::string &n) {
      if (nc == 1) {
        d.push_back("values of '" + n + "'");
      } else {
        for (size_type i = 0; i != nc; ++i) {
          d.push_back("values of the " + std::to_string(i) +
                      "th component of '" + n + "'");
        }
      }
    };
#ifdef MFEM_USE_MPI
    for (const auto &[n, f] : this->gridfunctions) {
      const auto nc = [&f] {
        if (std::holds_alternative<const GridFunction<false> *>(f)) {
          return getNumberOfComponents(
              *(std::get<const GridFunction<false> *>(f)));
        }
        return getNumberOfComponents(
            *(std::get<const GridFunction<true> *>(f)));
      }();
      add(nc, n);
    }
#else  /* MFEM_USE_MPI */
    for (const auto &[n, f] : this->gridfunctions) {
      const auto nc = getNumberOfComponents(*f);
      add(nc, n);
    }
#endif /* MFEM_USE_MPI */
    return d;
  }  // end of getValuesDescription

  std::optional<std::vector<std::vector<real>>> PointsSetCurves::getValues(
      Context &ctx, const TimeStepStage) const noexcept {
    auto r = std::vector<std::vector<real>>{};
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    auto add = [&r](const tfel::math::matrix<real> &values) {
      for (std::size_t nc = 0; nc != values.getNumberOfColumns(); ++nc) {
        auto nrow = std::vector<real>{};
        nrow.resize(values.getNumberOfRows());
        for (std::size_t nr = 0; nr != values.getNumberOfRows(); ++nr) {
          nrow[nr] = values(nr, nc);
        }
        r.push_back(std::move(nrow));
      }
    };
    //
    auto interpolator = GridFunctionInterpolator(this->fespaces_manager);
    if (std::holds_alternative<std::vector<Point<2>>>(points)) {
      const auto &pts = std::get<std::vector<Point<2>>>(points);
      if (!interpolator.addPoints(ctx, pts)) {
        return {};
      }
    } else {
      const auto &pts = std::get<std::vector<Point<3>>>(points);
      if (!interpolator.addPoints(ctx, pts)) {
        return {};
      }
    }
#ifdef MFEM_USE_MPI
    for (const auto &[n, f] : this->gridfunctions) {
      const auto ovalues = [&ctx, &f, &interpolator] {
        if (std::holds_alternative<const GridFunction<false> *>(f)) {
          return interpolator.interpolate(
              ctx, *(std::get<const GridFunction<false> *>(f)));
        }
        return interpolator.interpolate(
            ctx, *(std::get<const GridFunction<true> *>(f)));
      }();
      if (isInvalid(ovalues)) {
        return {};
      }
      add(*ovalues);
    }
#else  /* MFEM_USE_MPI */
    for (const auto &[n, f] : this->gridfunctions) {
      const auto ovalues = interpolator.interpolate(ctx, *f);
      if (isInvalid(ovalues)) {
        return {};
      }
      if (isInvalid(ovalues)) {
        return {};
      }
      add(*ovalues);
    }
#endif /* MFEM_USE_MPI */
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    return r;
  }  // end of getValues

  PointsSetCurves::~PointsSetCurves() noexcept = default;

}  // end of namespace mfem_mgis
