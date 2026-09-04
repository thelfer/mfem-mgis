/*!
 * \file   GeometryTest.cxx
 * \brief  Tests of the `GeometryTest` class
 * \author Thomas Helfer
 * \date   02/09/2026
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Math/Vector/tvectorIO.hxx"
#include "MFEMMGIS/Geometry.hxx"

struct {
  const char* mesh_file = nullptr;
  int order = 1;
  bool parallel = false;
} parameters;

struct GeometryTest final : public tfel::tests::TestCase {
  GeometryTest()
      : tfel::tests::TestCase("MFEMMGIS", "GeometryTest") {
  }  // end of GeometryTest

  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    return this->result;
  }  // end of execute

 private:
  // create a line in 2D using a uniform discretization
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    const auto params = Parameters{
        {"Line",
         Parameters{
             {"InitialPoint", std::vector<Parameter>{-2, 0.5}},  //
             {"FinalPoint", std::vector<Parameter>{6, 8}},
             {"Discretization",
              Parameters{{"Uniform", Parameters{{"NumberOfPoints", 10}}}}}}}};
    const auto opts = makePointsOnCurve<2>(ctx, params);
    TFEL_TESTS_ASSERT(isValid(opts));
    TFEL_TESTS_ASSERT(opts->size() == 10);
    const auto p0 = Point<2>{-2, 0.5};
    const auto p1 = Point<2>{6, 8};
    const auto v = (p1 - p0) / 9;
    for (std::size_t i = 0; i != 10; ++i) {
      const auto pa = opts->at(i);
      const auto pb = p0 + i * v;
      TFEL_TESTS_CHECK(norm(pb - pa) < 1e-14);
    }
  }  // end of test1
  // create a line in 3D
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    const auto params = Parameters{
        {"Line",
         Parameters{
             {"InitialPoint", std::vector<Parameter>{-2, 0.5, -4}},  //
             {"FinalPoint", std::vector<Parameter>{6, 8, 3.}},
             {"Discretization",
              Parameters{{"Uniform", Parameters{{"NumberOfPoints", 10}}}}}}}};
    const auto opts = makePointsOnCurve<3>(ctx, params);
    TFEL_TESTS_ASSERT(isValid(opts));
    TFEL_TESTS_ASSERT(opts->size() == 10);
    const auto p0 = Point<3>{-2, 0.5, -4};
    const auto p1 = Point<3>{6, 8, 3.};
    const auto v = (p1 - p0) / 9;
    for (std::size_t i = 0; i != 10; ++i) {
      const auto pa = opts->at(i);
      const auto pb = p0 + i * v;
      TFEL_TESTS_CHECK(norm(pb - pa) < 1e-14);
    }
  }  // end of test2
  // check that parsing issues are correctly detected
  void test3() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    // invalid generator
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Lne",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // invalid type of parameter for the Line generator
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(ctx, {{"Line", 1}})));
  }  // end of test3
  // check that parsing issues are correctly detected for the line generator
  void test4() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    // invalid type to define the initial point
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", 1},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // inconsistent space dimension
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5, -4}},
                          {"FinalPoint", std::vector<Parameter>{6, 8, 3.}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // unexpected parameter
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPont", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // inconsistent space dimension for the initial point
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5, 2}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // inconsistent space dimension for the final point
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8, 2}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // missing initial point
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // missing final point
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              10}}}}}}}})));
    // missing discretization
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx,
        {{"Line", Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                             {"FinalPoint", std::vector<Parameter>{6, 8}}}}})));
    // invalid discretization
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Unform", Parameters{{"NumberOfPoints",
                                                             10}}}}}}}})));
    // invalid parameter for the discretization
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOPoints",
                                                              10}}}}}}}})));
    // invalid number of points
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              1}}}}}}}})));
    // invalid number of points (2)
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization",
                           Parameters{{"Uniform", Parameters{{"NumberOfPoints",
                                                              -2}}}}}}}})));
    // invalid type of parameter for the discretization
    TFEL_TESTS_CHECK(isInvalid(makePointsOnCurve<2>(
        ctx, {{"Line",
               Parameters{{"InitialPoint", std::vector<Parameter>{-2, 0.5}},
                          {"FinalPoint", std::vector<Parameter>{6, 8}},
                          {"Discretization", Parameters{{"Uniform", 1}}}}}})));
  }  // end of test4
  // create a line in 2D using a discretization based on a geometric progression
  void test5() {
    using namespace mfem_mgis;
    static constexpr auto eps = real{1.e-8};
    static constexpr auto references = std::array<real, 11>{
        2.,         2.2320018,  2.55317371, 2.99778834, 3.613291, 4.46536267,
        5.64492887, 7.27786241, 9.53841545, 12.6678141, 17};
    auto ctx = Context{};
    auto exe = [this, &ctx](const real x0, const auto y0, const real x1,
                            const auto y1) {
      const auto params = Parameters{
          {"Line",
           Parameters{
               {"InitialPoint", std::vector<Parameter>{x0, y0}},  //
               {"FinalPoint", std::vector<Parameter>{x1, y1}},
               {"Discretization",
                Parameters{{"GeometricProgression",
                            Parameters{{"NumberOfPoints", 10},
                                       {"InitialDensity", real{1} / 150},
                                       {"FinalDensity", real{1} / 3}}}}}}}};
      const auto opts = makePointsOnCurve<2>(ctx, params);
      TFEL_TESTS_ASSERT(isValid(opts));
      TFEL_TESTS_ASSERT(opts->size() == 11u);
      return *opts;
    };
    const auto pts1 = exe(2, 0, 17, 0);
    const auto pts2 = exe(0, 2, 0, 17);
    for (std::size_t i = 0; i != 11; ++i) {
      TFEL_TESTS_CHECK(std::abs(pts1.at(i)[0] - references.at(i)) < eps);
      TFEL_TESTS_CHECK(std::abs(pts1.at(i)[1]) < eps);
      TFEL_TESTS_CHECK(std::abs(pts2.at(i)[1] - references.at(i)) < eps);
      TFEL_TESTS_CHECK(std::abs(pts2.at(i)[0]) < eps);
    }
  }  // end of test5
};

TFEL_TESTS_GENERATE_PROXY(GeometryTest, "GeometryTest");

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("GeometryTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
