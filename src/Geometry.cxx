/*!
 * \file   src/Geometry.cxx
 * \brief  This file implements the functions declared in
 * `MFEMMGIS/Geometry.hxx`
 * \author Thomas Helfer
 * \date   04/09/2026
 */

#include "TFEL/Math/Discretization1D.hxx"
#include "MFEMMGIS/Geometry.hxx"

namespace mfem_mgis::internals {

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] static std::optional<Point<N>> makePoint_impl(
          Context& ctx, const Parameter& p) noexcept {
    const auto opt = convert<std::vector<real>>(ctx, p);
    if (isInvalid(opt)) {
      return ctx.registerErrorMessage("can't extract point from parameter");
    }
    if (static_cast<size_type>(opt->size()) != N) {
      return ctx.registerErrorMessage(
          "can't extract point from parameter, invalid number of values");
    }
    if constexpr (N == 2) {
      return Point<2>{opt->at(0), opt->at(1)};
    } else {
      return Point<3>{opt->at(0), opt->at(1), opt->at(2)};
    }
  }  // end of makePoint_impl

  static std::optional<std::vector<real>> getUniformDensityWeights(
      Context& ctx, const Parameters& p) noexcept {
    if (!checkParameters(ctx, p,
                         std::map<std::string, std::string>{
                             {"NumberOfPoints", "number of points"}})) {
      return {};
    }
    if (!contains(p, "NumberOfPoints")) {
      return ctx.registerErrorMessage(
          "number of points undefined (no parameter 'NumberOfPoints')");
    }
    if (!is<int>(throwing, p, "NumberOfPoints")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'NumberOfPoints', expected an integer");
    }
    const auto n = get<int>(throwing, p, "NumberOfPoints");
    if (n < 2) {
      return ctx.registerErrorMessage("invalid number of point (" +
                                      std::to_string(n) + ")");
    }
    auto weights = std::vector<real>(static_cast<std::size_t>(n));
    const auto dx = static_cast<real>(1) / (n - 1);
    for (std::size_t i = 0; i != n - 1; ++i) {
      weights[i] = dx * i;
    }
    weights.back() = 1;
    return weights;
  }  // end of getUniformDensityWeights

  static std::optional<std::vector<real>>
  getDensityWeightsFollowingAGeometricProgression(
      Context& ctx, const Parameters& p) noexcept {
    if (!checkParameters(
            ctx, p,
            std::map<std::string, std::string>{
                {"InitialDensity", "targeted density at the inital point"},
                {"FinalDensity", "targeted density at the final point"},
                {"NumberOfPoints", "number of points"}})) {
      return {};
    }
    if (!contains(p, "InitialDensity")) {
      return ctx.registerErrorMessage(
          "intial density undefined (no parameter 'InitialDensity')");
    }
    if (!contains(p, "FinalDensity")) {
      return ctx.registerErrorMessage(
          "final density undefined (no parameter 'FinalDensity')");
    }
    if (!contains(p, "NumberOfPoints")) {
      return ctx.registerErrorMessage(
          "number of points undefined (no parameter 'NumberOfPoints')");
    }
    //
    if (!is<int>(throwing, p, "NumberOfPoints")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'NumberOfPoints', expected an integer");
    }
    const auto n = get<int>(throwing, p, "NumberOfPoints");
    if (n < 2) {
      return ctx.registerErrorMessage("invalid number of point (" +
                                      std::to_string(n) + ")");
    }
    if (!is<double>(throwing, p, "InitialDensity")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'InitialDensity', expected a floating "
          "point number");
    }
    if (!is<double>(throwing, p, "FinalDensity")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'FinalDensity', expected a floating "
          "point number");
    }
    const auto di = get<double>(throwing, p, "InitialDensity");
    const auto de = get<double>(throwing, p, "FinalDensity");
    if (!(di > 0)) {
      return ctx.registerErrorMessage(
          "invalid initial density, expected a positive number");
    }
    if (!(de > 0)) {
      return ctx.registerErrorMessage(
          "invalid final density, expected a positive number");
    }
    auto weights = std::vector<real>();
    const auto ok = MGIS_INVOKE(
        ctx, ::tfel::math::geometricDiscretization<std::vector<real>>, weights,
        0, 1, di, de, n);
    if (!ok) {
      return {};
    }
    weights.back() = 1;
    return weights;
  }  // end of getDensityWeightsFollowingAGeometricProgression

  static std::optional<std::vector<real>> getDiscretizationWeights(
      Context& ctx, const Parameters& p) noexcept {
    const auto odiscretization = extractFactoryArgument(ctx, p);
    if (isInvalid(odiscretization)) {
      return {};
    }
    const auto [n, params] = *odiscretization;
    if (n == "Uniform") {
      return getUniformDensityWeights(ctx, params);
    } else if (n == "GeometricProgression") {
      return getDensityWeightsFollowingAGeometricProgression(ctx, params);
    }
    return ctx.registerErrorMessage(
        "unknown curve discretization generator '" + n +
        "'. Currently the only supported generators are 'Uniform' and "
        "'GeometricProgression'");
  }  // end of getDiscretizationWeights

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] static std::optional<std::vector<Point<N>>> makeLine_impl(
          Context& ctx, const Parameters& p) noexcept {
    if (!checkParameters(
            ctx, p,
            std::map<std::string, std::string>{
                {"InitialPoint", "initial point"},
                {"FinalPoint", "final point"},
                {"Discretization",
                 "description of the discretization of the line"}})) {
      return {};
    }
    if (!contains(p, "InitialPoint")) {
      return ctx.registerErrorMessage(
          "intial point undefined (no parameter named 'InitialPoint')");
    }
    if (!contains(p, "FinalPoint")) {
      return ctx.registerErrorMessage(
          "final point undefined (no parameter named 'FinalPoint')");
    }
    if (!contains(p, "Discretization")) {
      return ctx.registerErrorMessage(
          "discretization undefined (no parameter named 'Discretization')");
    }
    if (!is<Parameters>(throwing, p, "Discretization")) {
      return ctx.registerErrorMessage(
          "invalid type for the parameter 'Discretization'");
    }
    const auto op0 = makePoint_impl<N>(ctx, get(throwing, p, "InitialPoint"));
    const auto op1 = makePoint_impl<N>(ctx, get(throwing, p, "FinalPoint"));
    const auto oweights = getDiscretizationWeights(
        ctx, get<Parameters>(throwing, p, "Discretization"));
    if (!areValid(op0, op1, oweights)) {
      return {};
    }
    ctx.assertOrTerminate(oweights->size() >= 2,
                          "internal error: invalid discretization weights");
    const auto dp = *op1 - *op0;
    auto pts = std::vector<Point<N>>(oweights->size());
    for (std::size_t i = 0; i != oweights->size(); ++i) {
      pts[i] = *op0 + oweights->at(i) * dp;
    }
    pts.back() = *op1;
    return pts;
  }  // end of makeLine_impl

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] static std::optional<
          std::vector<Point<N>>> makePointsOnCurve_impl(Context& ctx,
                                                        const Parameters&
                                                            p) noexcept {
    const auto ostrategy = extractFactoryArgument(ctx, p);
    if (isInvalid(ostrategy)) {
      return {};
    }
    const auto [n, params] = *ostrategy;
    if (n == "Line") {
      return makeLine_impl<N>(ctx, params);
    }
    return ctx.registerErrorMessage(
        "unknown curve generator '" + n +
        "'. Currently the only supported generator is 'Line'");
  }  // end of makePointsOnCurve_impl

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  template <>
  std::optional<Point<2>> makePoint<2>(Context& ctx,
                                       const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<2>(ctx, p);
  }  // end of makePoint

  template <>
  std::optional<Point<3>> makePoint<3>(Context& ctx,
                                       const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<3>(ctx, p);
  }  // end of makePoint

  template <>
  std::optional<std::vector<Point<2>>> makePointsOnCurve<2>(
      Context& ctx, const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<2>(ctx, p);
  }  // end of makePointsOnCurve

  template <>
  std::optional<std::vector<Point<3>>> makePointsOnCurve<3>(
      Context& ctx, const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<3>(ctx, p);
  }  // end of makePointsOnCurve

}  // end of namespace mfem_mgis
