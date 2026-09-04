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
#include "MFEMMGIS/GridFunctionInterpolator.hxx"

namespace mfem_mgis::internals {

  template <size_type N>
  requires((N == 2) || (N == 3)) [[nodiscard]] static bool setSpaceDimension(
      Context& ctx, std::optional<size_type>& space_dimension) noexcept {
    if (!space_dimension.has_value()) {
      space_dimension = N;
      return true;
    }
    if (*space_dimension != N) {
      return ctx.registerErrorMessage(
          "the space dimension is already can't be set to " +
          std::to_string(N) + "as it is already set to " +
          std::to_string(*space_dimension));
    }
    return true;
  }  // end of setSpaceDimension

  template <size_type N>
  requires((N == 2) || (N == 3)) [[nodiscard]] static bool addPoints_impl(
      Context& ctx,
      std::optional<size_type>& space_dimension,
      std::vector<real>& points,
      const std::vector<Point<N>>& pts) noexcept {
    if (!setSpaceDimension<N>(ctx, space_dimension)) {
      return ctx.registerErrorMessage(
          "GridFunctionInterpolator::addPoints failed");
    }
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
    return true;
  }  // end of addPoints_impl

  template <bool parallel, size_type N>
  requires((N == 2) || (N == 3))
      [[nodiscard]] static std::optional<tfel::math::matrix<real>>  //
      interpolate_impl(Context& ctx,
                       GridFunction<parallel>& f,
                       std::vector<real>& points) {
    auto* const fespace = [&f] {
      if constexpr (parallel) {
        return f.ParFESpace();
      } else {
        return f.FESpace();
      }
    }();
    auto* const mesh = [&fespace] {
      if constexpr (parallel) {
        return fespace->GetParMesh();
      } else {
        return fespace->GetMesh();
      }
    }();
    if (mesh->SpaceDimension() != N) {
      return ctx.registerErrorMessage(
          "a grid function defined on a mesh with space dimension " +
          std::to_string(mesh->SpaceDimension()) +
          " can't be interpolation of points of dimension '" +
          std::to_string(N) + "'");
    }
    //
    // set the node fespace if required
    const auto shall_set_nodal_fespace = [mesh, fespace] {
      // nodes is a pointer to a grid function, even in parallel
      const auto* const nodes = mesh->GetNodes();
      if (nodes == nullptr) {
        return true;
      }
      return nodes->FESpace() != fespace;
    }();
    if (shall_set_nodal_fespace) {
      mesh->SetNodalFESpace(const_cast<FiniteElementSpace<parallel>*>(fespace));
    }
    //
    auto finder = mfem::FindPointsGSLIB{};
    finder.Setup(*mesh);
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

  GridFunctionInterpolator::GridFunctionInterpolator() = default;

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx, const std::vector<Point<2>>& pts) {
    auto or_raise = ctx.getThrowingFailureHandler();
    this->addPoints(ctx, pts) | or_raise;
  }  // end of GridFunctionInterpolator

  GridFunctionInterpolator::GridFunctionInterpolator(
      Context& ctx, const std::vector<Point<3>>& pts) {
    auto or_raise = ctx.getThrowingFailureHandler();
    this->addPoints(ctx, pts) | or_raise;
  }  // end of GridFunctionInterpolator

  bool GridFunctionInterpolator::addPoints(
      Context& ctx, const std::vector<Point<2>>& pts) noexcept {
    return internals::addPoints_impl<2>(ctx, this->space_dimension,
                                        this->points, pts);
  }  // end of GridFunctionInterpolator::addPoints

  bool GridFunctionInterpolator::addPoints(
      Context& ctx, const std::vector<Point<3>>& pts) noexcept {
    return internals::addPoints_impl<3>(ctx, this->space_dimension,
                                        this->points, pts);
  }  // end of GridFunctionInterpolator::addPoints

#ifdef MFEM_USE_MPI
  std::optional<tfel::math::matrix<real>> GridFunctionInterpolator::interpolate(
      Context& ctx, GridFunction<true>& f) noexcept {
    if (!this->space_dimension.has_value()) {
      return ctx.registerErrorMessage("no points defined");
    }
    if (*(this->space_dimension) == 2) {
      return internals::interpolate_impl<true, 2>(ctx, f, this->points);
    }
    return internals::interpolate_impl<true, 3>(ctx, f, this->points);
  }    // end of interpolate
#endif /* MFEM_USE_MPI */

  std::optional<tfel::math::matrix<real>> GridFunctionInterpolator::interpolate(
      Context& ctx, GridFunction<false>& f) noexcept {
    if (!this->space_dimension.has_value()) {
      return ctx.registerErrorMessage("no points defined");
    }
    if (*(this->space_dimension) == 2) {
      return internals::interpolate_impl<false, 2>(ctx, f, this->points);
    }
    return internals::interpolate_impl<false, 3>(ctx, f, this->points);
  }  // end of interpolate

  GridFunctionInterpolator::~GridFunctionInterpolator() = default;

}  // end of namespace mfem_mgis
