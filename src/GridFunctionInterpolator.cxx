/*!
 * \file   src/GridFunctionInterpolator.cxx
 * \brief  This file implements the methods of the `GridFunctionInterpolator`
 * class.
 * \author Thomas Helfer
 * \date   04/09/2026
 */

#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pgridfunc.hpp"
#endif /* MFEM_USE_MPI */
#include "mfem/fem/gslib.hpp"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/GridFunctionInterpolator.hxx"

namespace mfem_mgis::internals {

  template <size_type N>
  requires((N == 2) || (N == 3)) static void addPoints_impl(
      std::vector<real>& points, const std::vector<Point<N>>& pts) noexcept {
    auto n = points.size();
    points.resize(points.size() + N * pts.size());
    for (const auto& p : pts) {
      points[n] = p[0];
      points[n + 1] = p[1];
      if constexpr (N == 3) {
        points[n + 2] = p[2];
      }
      n += N;
    }
  }  // end of addPoints_impl

  template <bool parallel, size_type N>
  requires((N == 2) || (N == 3))
      [[nodiscard]] static std::optional<tfel::math::matrix<real>>  //
      interpolate_impl(Context& ctx,
                       const FiniteElementSpacesManager& fespaces_manager,
                       const GridFunction<parallel>& f,
                       std::vector<real>& points) {
    const auto& mesh = fespaces_manager.getMeshDiscretization();
    const auto d = getSpaceDimension(mesh);
    const auto* fespace = [&f] {
      if constexpr (parallel) {
        return f.ParFESpace();
      } else {
        return f.FESpace();
      }
    }();
    if (d != N) {
      return ctx.registerErrorMessage(
          "a grid function defined on a mesh with space dimension " +
          std::to_string(d) +
          " can't be interpolation of points of dimension '" +
          std::to_string(N) + "'");
    }
    //
    if (!fespaces_manager.setNodalFiniteElementSpace(ctx)) {
      return {};
    }
    //
    auto finder = mfem::FindPointsGSLIB{};
    finder.Setup(*(mesh.getMutableMeshPointer<parallel>()));
    finder.SetDefaultInterpolationValue(std::numeric_limits<real>::quiet_NaN());
    //
    auto pts =
        mfem::Vector{points.data(), static_cast<size_type>(points.size())};
    finder.FindPoints(pts, mfem::Ordering::byVDIM);
    for (const auto& c : finder.GetCode()) {
      if (c == 2) {
        return ctx.registerErrorMessage("some points were not found");
      }
    }
    const auto npoints = static_cast<size_type>(points.size() / N);
    auto result = std::optional<tfel::math::matrix<real>>{};
    result = tfel::math::matrix<real>(npoints, fespace->GetVDim());
    auto tmp = std::vector<real>{};
    auto data = static_cast<real*>(nullptr);
    if (fespace->GetOrdering() == mfem::Ordering::byVDIM) {
      data = result->data();
    } else {
      tmp.resize(static_cast<std::vector<real>::size_type>(npoints *
                                                           fespace->GetVDim()));
      data = tmp.data();
    }
    auto interpolated_values = mfem::Vector(
        data, static_cast<size_type>(npoints * fespace->GetVDim()));
    finder.Interpolate(f, interpolated_values);
    if (fespace->GetOrdering() == mfem::Ordering::byNODES) {
      for (size_type i = 0; i != npoints; ++i) {
        for (size_type j = 0; j != fespace->GetVDim(); ++j) {
          const auto idx = static_cast<size_type>(j * npoints + i);
          result->operator()(i, j) = tmp[idx];
        }
      }
    }
    return result;
  }  // end of interpolate_impl

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  GridFunctionInterpolator::GridFunctionInterpolator(
      const FiniteElementSpacesManager& m) noexcept
      : fespaces_manager(m) {}

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx,
      const FiniteElementSpacesManager& m,
      const std::vector<Point<2>>& pts)
      : GridFunctionInterpolator(m) {
    auto or_raise = ctx.getThrowingFailureHandler();
    this->addPoints(ctx, pts) | or_raise;
  }  // end of GridFunctionInterpolator

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx,
      const FiniteElementSpacesManager& m,
      const std::vector<Point<3>>& pts)
      : GridFunctionInterpolator(m) {
    auto or_raise = ctx.getThrowingFailureHandler();
    this->addPoints(ctx, pts) | or_raise;
  }  // end of GridFunctionInterpolator

  GridFunctionInterpolator::GridFunctionInterpolator(
      const FiniteElementDiscretization& fed) noexcept
      : fespaces_manager(fed.getFiniteElementSpacesManager()) {}

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const std::vector<Point<2>>& pts)
      : GridFunctionInterpolator(
            ctx, fed.getFiniteElementSpacesManager(), pts) {
  }  // end of GridFunctionInterpolator

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const std::vector<Point<3>>& pts)
      : GridFunctionInterpolator(
            ctx, fed.getFiniteElementSpacesManager(), pts) {
  }  // end of GridFunctionInterpolator

  bool GridFunctionInterpolator::addPoints(
      Context& ctx, const std::vector<Point<2>>& pts) noexcept {
    const auto d =
        getSpaceDimension(this->fespaces_manager.getMeshDiscretization());
    if (d != 2) {
      return ctx.registerErrorMessage("can't add 2D points to a " +
                                      std::to_string(d) + "D mesh");
    }
    internals::addPoints_impl<2>(this->points, pts);
    return true;
  }  // end of GridFunctionInterpolator::addPoints

  bool GridFunctionInterpolator::addPoints(
      Context& ctx, const std::vector<Point<3>>& pts) noexcept {
    const auto d =
        getSpaceDimension(this->fespaces_manager.getMeshDiscretization());
    if (d != 3) {
      return ctx.registerErrorMessage("can't add 3D points to a " +
                                      std::to_string(d) + "D mesh");
    }
    internals::addPoints_impl<3>(this->points, pts);
    return true;
  }  // end of GridFunctionInterpolator::addPoints

#ifdef MFEM_USE_MPI
  std::optional<tfel::math::matrix<real>> GridFunctionInterpolator::interpolate(
      Context& ctx, const GridFunction<true>& f) noexcept {
    if (!this->fespaces_manager.manages(*(f.ParFESpace()))) {
      return ctx.registerErrorMessage(
          "the given grid function is defined on a finite element space which "
          "is not managed by the finite element spaces manager of which the "
          "interpolator is built");
    }
    const auto d =
        getSpaceDimension(this->fespaces_manager.getMeshDiscretization());
    if (d == 2) {
      return internals::interpolate_impl<true, 2>(ctx, this->fespaces_manager,
                                                  f, this->points);
    }
    return internals::interpolate_impl<true, 3>(ctx, this->fespaces_manager, f,
                                                this->points);
  }    // end of interpolate
#endif /* MFEM_USE_MPI */

  std::optional<tfel::math::matrix<real>> GridFunctionInterpolator::interpolate(
      Context& ctx, const GridFunction<false>& f) noexcept {
    if (!this->fespaces_manager.manages(*(f.FESpace()))) {
      return ctx.registerErrorMessage(
          "the given grid function is defined on a finite element space which "
          "is not managed by the finite element spaces manager of which the "
          "interpolator is built");
    }
    const auto d =
        getSpaceDimension(this->fespaces_manager.getMeshDiscretization());
    if (d == 2) {
      return internals::interpolate_impl<false, 2>(ctx, this->fespaces_manager,
                                                   f, this->points);
    }
    return internals::interpolate_impl<false, 3>(ctx, this->fespaces_manager, f,
                                                 this->points);
  }  // end of interpolate

  GridFunctionInterpolator::~GridFunctionInterpolator() = default;

}  // end of namespace mfem_mgis
