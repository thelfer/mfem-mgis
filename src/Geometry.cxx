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
      static std::optional<Point<N>> makePoint_impl(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept {
    if (is<std::string>(p)) {
      const auto n = get<std::string>(throwing, p);
      if (pts.contains(n)) {
        return pts.at(n);
      }
      return ctx.registerErrorMessage("no point named '" + n + "' registred");
    }
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
    if (!checkParameters(ctx, p,
                         std::map<std::string, std::string>{
                             {"NormalizedInitialDensity",
                              "targeted density at the inital point"},
                             {"NormalizedFinalDensity",
                              "targeted density at the final point"},
                             {"NumberOfPoints", "number of points"}})) {
      return {};
    }
    if (!contains(p, "NormalizedInitialDensity")) {
      return ctx.registerErrorMessage(
          "intial density undefined (no parameter 'NormalizedInitialDensity')");
    }
    if (!contains(p, "NormalizedFinalDensity")) {
      return ctx.registerErrorMessage(
          "final density undefined (no parameter 'NormalizedFinalDensity')");
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
    if (!is<double>(throwing, p, "NormalizedInitialDensity")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'NormalizedInitialDensity', expected a "
          "floating "
          "point number");
    }
    if (!is<double>(throwing, p, "NormalizedFinalDensity")) {
      return ctx.registerErrorMessage(
          "invalid type for parameter 'NormalizedFinalDensity', expected a "
          "floating "
          "point number");
    }
    const auto di = get<double>(throwing, p, "NormalizedInitialDensity");
    const auto de = get<double>(throwing, p, "NormalizedFinalDensity");
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
      static std::optional<std::vector<Point<N>>> makeLine_impl(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameters& p) noexcept {
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
    const auto op0 =
        makePoint_impl<N>(ctx, pts, get(throwing, p, "InitialPoint"));
    const auto op1 =
        makePoint_impl<N>(ctx, pts, get(throwing, p, "FinalPoint"));
    const auto oweights = getDiscretizationWeights(
        ctx, get<Parameters>(throwing, p, "Discretization"));
    if (!areValid(op0, op1, oweights)) {
      return {};
    }
    ctx.assertOrTerminate(oweights->size() >= 2,
                          "internal error: invalid discretization weights");
    const auto dp = *op1 - *op0;
    auto npts = std::vector<Point<N>>(oweights->size());
    for (std::size_t i = 0; i != oweights->size(); ++i) {
      npts[i] = *op0 + oweights->at(i) * dp;
    }
    npts.back() = *op1;
    return npts;
  }  // end of makeLine_impl

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      static std::optional<std::vector<Point<N>>> makePointsOnCurve_impl(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameters& p) noexcept {
    const auto ostrategy = extractFactoryArgument(ctx, p);
    if (isInvalid(ostrategy)) {
      return {};
    }
    const auto [n, params] = *ostrategy;
    if (n == "Line") {
      return makeLine_impl<N>(ctx, pts, params);
    }
    return ctx.registerErrorMessage(
        "unknown curve generator '" + n +
        "'. Currently the only supported generator is 'Line'");
  }  // end of makePointsOnCurve_impl

  template <size_type N>
  requires((N == 2) || (N == 3))                   //
      static std::optional<std::vector<Point<N>>>  //
      makePointsSet_impl(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept {
    if (is<Parameters>(p)) {
      return makePointsOnCurve_impl<N>(ctx, pts, get<Parameters>(throwing, p));
    }
    if (!is<std::vector<Parameter>>(p)) {
      return ctx.registerErrorMessage(
          "can't extract points set from parameter (invalid parameter type)");
    }
    const auto& parameters = get<std::vector<Parameter>>(throwing, p);
    auto points = std::vector<Point<N>>{};
    points.reserve(parameters.size());
    for (const auto& parameter : parameters) {
      const auto opt = makePoint<N>(ctx, pts, parameter);
      if (isInvalid(opt)) {
        return {};
      }
      points.push_back(*opt);
    }
    return points;
  }  // end of makePointsSet_impl

}  // end of namespace mfem_mgis::internals

namespace mfem_mgis {

  template <>
  std::optional<Point<2>> makePoint<2>(Context& ctx,
                                       const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<2>(ctx, {}, p);
  }  // end of makePoint

  template <>
  std::optional<Point<3>> makePoint<3>(Context& ctx,
                                       const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<3>(ctx, {}, p);
  }  // end of makePoint

  template <>
  std::optional<std::vector<Point<2>>> makePointsSet<2>(
      Context& ctx, const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePointsSet_impl<2>(ctx, {}, p);
  }  // end of makePointsSet

  template <>
  std::optional<std::vector<Point<3>>> makePointsSet<3>(
      Context& ctx, const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePointsSet_impl<3>(ctx, {}, p);
  }  // end of makePointsSet

  template <>
  std::optional<std::vector<Point<2>>> makePointsOnCurve<2>(
      Context& ctx, const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<2>(ctx, {}, p);
  }  // end of makePointsOnCurve

  template <>
  std::optional<std::vector<Point<3>>> makePointsOnCurve<3>(
      Context& ctx, const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<3>(ctx, {}, p);
  }  // end of makePointsOnCurve

  template <>
  std::string toString<2>(const Point<2>& pt) noexcept {
    return '(' + std::to_string(pt[0]) + ", " + std::to_string(pt[1]) + ')';
  }  // end of toString

  template <>
  std::string toString<3>(const Point<3>& pt) noexcept {
    return '(' + std::to_string(pt[0]) + ", " +  //
           std::to_string(pt[1]) + ", " +        //
           std::to_string(pt[2]) + ')';
  }  // end of toString

  template <>
  std::optional<Point<2>> makePoint<2>(
      Context& ctx,
      const std::map<std::string, Point<2>, std::less<>>& pts,
      const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<2>(ctx, pts, p);
  }  // end of makePoint<2>

  template <>
  std::optional<std::vector<Point<2>>> makePointsSet<2>(
      Context& ctx,
      const std::map<std::string, Point<2>, std::less<>>& pts,
      const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePointsSet_impl<2>(ctx, pts, p);
  }  // end of makePointsSet<2>

  template <>
  std::optional<std::vector<Point<2>>> makePointsOnCurve<2>(
      Context& ctx,
      const std::map<std::string, Point<2>, std::less<>>& pts,
      const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<2>(ctx, pts, p);
  }  // end of makePointsOnCurve<2>

  template <>
  std::optional<Point<3>> makePoint<3>(
      Context& ctx,
      const std::map<std::string, Point<3>, std::less<>>& pts,
      const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePoint_impl<3>(ctx, pts, p);
  }  // end of makePoint<3>

  template <>
  std::optional<std::vector<Point<3>>> makePointsSet<3>(
      Context& ctx,
      const std::map<std::string, Point<3>, std::less<>>& pts,
      const Parameter& p) noexcept {
    return ::mfem_mgis::internals::makePointsSet_impl<3>(ctx, pts, p);
  }  // end of makePointsSet<3>

  template <>
  std::optional<std::vector<Point<3>>> makePointsOnCurve<3>(
      Context& ctx,
      const std::map<std::string, Point<3>, std::less<>>& pts,
      const Parameters& p) noexcept {
    return ::mfem_mgis::internals::makePointsOnCurve_impl<3>(ctx, pts, p);
  }  // end of makePointsOnCurve<3>

}  // end of namespace mfem_mgis
