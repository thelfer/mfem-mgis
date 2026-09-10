/*!
 * \file   src/GridFunctionValuesCurve.cxx
 * \brief  This file implements the `GridFunctionValuesCurve` class
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pgridfunc.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
#include "MFEMMGIS/GridFunctionInterpolator.hxx"
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
#include "MFEMMGIS/PostProcessing/GridFunctionValuesCurve.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string> getParametersDescription() noexcept {
    return {{"Points",
             "list of points on which the grid function is interpolated"}};
  }  // end of getParametersDescription()

  std::string GridFunctionValuesCurve::getDescription() noexcept {
    return "curve retrieving the values of the interoplation of a grid "
           "function at a set of points";
  }  // end of getDescription

  GridFunctionValuesCurve::GridFunctionValuesCurve(
      Context &ctx,
      PhysicalSystem &ps,
      const FiniteElementSpacesManager &manager,
      const Parameters &parameters)
      : physicalSystem(ps), fespaces_manager(manager) {
#pragma message("shall check that ps and manager are compatible")
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    auto or_raise = ctx.getThrowingFailureHandler();
    checkParameters(throwing, parameters,
                    GridFunctionValuesCurve::getParametersDescription());
    if (contains(parameters, "Points")) {
      const auto &m = this->physicalSystem.getMeshDiscretization();
      const auto d = getSpaceDimension(m);
      if (d == 2) {
        this->points =
            makePointsSet<2>(ctx, m, get(throwing, parameters, "Points")) |
            or_raise;
      } else if (d == 3) {
        this->points =
            makePointsSet<3>(ctx, m, get(throwing, parameters, "Points")) |
            or_raise;
      } else {
        raise("unsupported space dimension");
      }
    }
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
  }    // end of GridFunctionValuesCurve

  bool GridFunctionValuesCurve::setGridRunction(
      Context &ctx, std::string_view n, const GridFunction<true> &f) noexcept {
    if ((this->parallel_fct != nullptr) || (this->sequential_fct != nullptr)) {
      return ctx.registerErrorMessage("grid function already set");
    }
    this->name = n;
    this->parallel_fct = &f;
    return true;
  }  // end of setGridRunction

  bool GridFunctionValuesCurve::setGridRunction(
      Context &ctx, std::string_view n, const GridFunction<false> &f) noexcept {
    if ((this->parallel_fct != nullptr) || (this->sequential_fct != nullptr)) {
      return ctx.registerErrorMessage("grid function already set");
    }
    this->name = n;
    this->sequential_fct = &f;
    return true;
  }  // end of setGridRunction

  bool GridFunctionValuesCurve::addPoints(Context &ctx,
                                          const Parameter &parameter) noexcept {
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    const auto &m = this->physicalSystem.getMeshDiscretization();
    const auto d = getSpaceDimension(m);
    if (d == 2) {
      const auto opoints = makePointsSet<2>(ctx, parameter);
      if (isInvalid(opoints)) {
        return {};
      }
      this->points = *opoints;
    } else if (d == 3) {
      const auto opoints = makePointsSet<3>(ctx, parameter);
      if (isInvalid(opoints)) {
        return {};
      }
      this->points = *opoints;
    } else {
      return ctx.registerErrorMessage("unsupported space dimension");
    }
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    return true;
  }  // end of addPoints

  std::vector<std::string> GridFunctionValuesCurve::getDescriptions()
      const noexcept {
    if ((this->parallel_fct == nullptr) && (this->sequential_fct == nullptr)) {
      return {};
    }
    if (!this->arePointsDefined()) {
      return {};
    }
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    const auto nc = [this]() {
      if (this->parallel_fct != nullptr) {
        return this->parallel_fct->ParFESpace()->GetVDim();
      }
      return this->sequential_fct->FESpace()->GetVDim();
    }();
    auto d = std::vector<std::string>{};
    if (std::holds_alternative<std::vector<Point<2>>>(this->points)) {
      for (const auto &p : std::get<std::vector<Point<2>>>(this->points)) {
        const auto pt = toString(p);
        if (nc == 1) {
          d.push_back("value of function '" + this->name + "' at point " + pt);
        } else {
          for (size_type i = 0; i != nc; ++i) {
            d.push_back("value of component " + std::to_string(i) +
                        " of function '" + this->name + "' at point " + pt);
          }
        }
      }
    } else {
      for (const auto &p : std::get<std::vector<Point<3>>>(this->points)) {
        const auto pt = toString(p);
        if (nc == 1) {
          d.push_back("value of function '" + this->name + "' at point " + pt);
        } else {
          for (size_type i = 0; i != nc; ++i) {
            d.push_back("value of component " + std::to_string(i) +
                        " of function '" + this->name + "' at point " + pt);
          }
        }
      }
    }
    return d;
#else  /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    return {};
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
  }    // end of getDescriptions

  std::optional<std::vector<real>> GridFunctionValuesCurve::getValues(
      Context &ctx, const TimeStepStage) const noexcept {
    auto oresults = std::optional<std::vector<real>>{};
    if ((this->parallel_fct == nullptr) && (this->sequential_fct == nullptr)) {
      return oresults;
    }
    if (!this->arePointsDefined()) {
      return oresults;
    }
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    auto ointerpolator = [this, &ctx] {
      if (std::holds_alternative<std::vector<Point<2>>>(this->points)) {
        return construct<GridFunctionInterpolator>(
            ctx, this->fespaces_manager,
            std::get<std::vector<Point<2>>>(this->points));
      }
      return construct<GridFunctionInterpolator>(
          ctx, this->fespaces_manager,
          std::get<std::vector<Point<3>>>(this->points));
    }();
    if (isInvalid(ointerpolator)) {
      return oresults;
    }
    const auto ovalues = [this, &ctx, &ointerpolator] {
      if (this->parallel_fct != nullptr) {
        return ointerpolator->interpolate(ctx, *(this->parallel_fct));
      }
      return ointerpolator->interpolate(ctx, *(this->sequential_fct));
    }();
    if (isInvalid(ovalues)) {
      return oresults;
    }
    const auto s = static_cast<std::size_t>(ovalues->getNumberOfRows() *
                                            ovalues->getNumberOfColumns());
    oresults = std::vector<real>{ovalues->data(), ovalues->data() + s};
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    return oresults;
  }  // end of getValues

  bool GridFunctionValuesCurve::arePointsDefined() const noexcept {
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    if (std::holds_alternative<std::vector<Point<2>>>(this->points)) {
      return !std::get<std::vector<Point<2>>>(this->points).empty();
    }
    return !std::get<std::vector<Point<3>>>(this->points).empty();
#else  /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    return false;
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
  }    // end of arePointsDefined

  GridFunctionValuesCurve::~GridFunctionValuesCurve() noexcept = default;

}  // end of namespace mfem_mgis