/*!
 * \file   GridFunctionInterpolatorTest.cxx
 * \brief  Tests of the `GridFunctionInterpolatorTest` class
 * \author Thomas Helfer
 * \date   02/09/2026
 */

#include <cmath>
#include <limits>
#include <vector>
#include <variant>
#include <cstdlib>
#include <optional>
#include "mfem/general/optparser.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/gslib.hpp"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/GridFunctionUtilities.hxx"
#include "MFEMMGIS/GridFunctionInterpolator.hxx"

struct {
  const char* mesh_file = nullptr;
  int order = 1;
  bool parallel = false;
} parameters;

struct GridFunctionInterpolatorTest final : public tfel::tests::TestCase {
  GridFunctionInterpolatorTest()
      : tfel::tests::TestCase("MFEMMGIS", "GridFunctionInterpolatorTest") {
  }  // end of GridFunctionInterpolatorTest
  tfel::tests::TestResult execute() override {
    if (parameters.parallel) {
#ifdef MFEM_USE_MPI
      this->test1<true>();
      this->test2<true>();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
      this->test2<false>();
    }
    return this->result;
  }

 private:
  //
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto ofed = construct<FiniteElementDiscretization>(
        ctx,
        Parameters{{"MeshFileName", parameters.mesh_file},
                   {"FiniteElementFamily", "H1"},
                   {"FiniteElementOrder", parameters.order},
                   {"UnknownsSize", 2},
                   {"NumberOfUniformRefinements", 0},  // faster for testing
                   {"Parallel", parameters.parallel}});
    TFEL_TESTS_ASSERT(isValid(ofed));
    auto& fespace = ofed->template getFiniteElementSpace<parallel>();
    auto& m = ofed->template getMesh<parallel>();
    m.SetNodalFESpace(&fespace);
    //
    auto ox = makeGridFunction<parallel>(ctx, fespace, 2);
    TFEL_TESTS_ASSERT(isValid(ox));
    m.SetNodalGridFunction(ox->f.get());
    //
    mfem::FindPointsGSLIB finder;
    finder.Setup(m);
    finder.SetDefaultInterpolationValue(std::numeric_limits<real>::quiet_NaN());
    //
    auto pts = mfem::Vector(6);
    pts[0] = 0.5830000000000004;
    pts[1] = 0.2;
    pts[2] = 0.5830000000000004;
    pts[3] = 0.1;
    pts[4] = 0.5830000000000004;
    pts[5] = 0.4;
    finder.FindPoints(pts, mfem::Ordering::byVDIM);
    //
    const auto codes = finder.GetCode();
    TFEL_TESTS_ASSERT(codes.Size() == 3);
    TFEL_TESTS_CHECK(codes[0] != 2);
    TFEL_TESTS_CHECK(codes[1] != 2);
    TFEL_TESTS_CHECK(codes[2] == 2);
    //
    auto x_values = mfem::Vector(4);
    finder.Interpolate(*(ox->f), x_values);
    if (fespace.GetOrdering() == mfem::Ordering::byNODES) {
      for (size_type i = 0; i != 2; ++i) {
        TFEL_TESTS_CHECK(std::abs(x_values[i] - pts[2 * i]) < 1e-14);
        TFEL_TESTS_CHECK(std::abs(x_values[3 + i] - pts[2 * i + 1]) < 1e-14);
      }
      TFEL_TESTS_CHECK(std::isnan(x_values[2]));
      TFEL_TESTS_CHECK(std::isnan(x_values[5]));
    } else {
      for (size_type i = 0; i != 2; ++i) {
        TFEL_TESTS_CHECK(std::abs(x_values[2 * i] - pts[2 * i]) < 1e-14);
        TFEL_TESTS_CHECK(std::abs(x_values[2 * i + 1] - pts[2 * i + 1]) <
                         1e-14);
      }
      TFEL_TESTS_CHECK(std::isnan(x_values[4]));
      TFEL_TESTS_CHECK(std::isnan(x_values[5]));
    }
  }  // end of test1
     //
  template <bool parallel>
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto ofed = construct<FiniteElementDiscretization>(
        ctx,
        Parameters{{"MeshFileName", parameters.mesh_file},
                   {"FiniteElementFamily", "H1"},
                   {"FiniteElementOrder", parameters.order},
                   {"UnknownsSize", 2},
                   {"NumberOfUniformRefinements", 0},  // faster for testing
                   {"Parallel", parameters.parallel}});
    TFEL_TESTS_ASSERT(isValid(ofed));
    auto& fespace = ofed->template getFiniteElementSpace<parallel>();
    auto& m = ofed->template getMesh<parallel>();
    //
    auto ocoords = makeGridFunction<parallel>(ctx, fespace, 2);
    TFEL_TESTS_ASSERT(isValid(ocoords));
    m.SetNodalGridFunction(ocoords->f.get());
    //
    const auto pts = std::vector<GridFunctionInterpolator::Point<2>>{
        {0.583, 0.2}, {0.583, 0.1}};
    auto ointerpolator = construct<GridFunctionInterpolator>(ctx, pts);
    TFEL_TESTS_ASSERT(isValid(ointerpolator));
    const auto ovalues = ointerpolator->interpolate(ctx, *(ocoords->f));
    TFEL_TESTS_ASSERT(isValid(ovalues));
    TFEL_TESTS_CHECK_EQUAL(ovalues->getNumberOfRows(), 2);
    TFEL_TESTS_CHECK_EQUAL(ovalues->getNumberOfColumns(), 2);
    for (size_type i = 0; i != 2; ++i) {
      const auto x = ovalues->operator()(i, 0);
      const auto y = ovalues->operator()(i, 1);
      TFEL_TESTS_CHECK(std::abs(x - pts[i][0]) < 1e-14);
      TFEL_TESTS_CHECK(std::abs(y - pts[i][1]) < 1e-14);
    }
  }  // end of test2
};

TFEL_TESTS_GENERATE_PROXY(GridFunctionInterpolatorTest,
                          "GridFunctionInterpolatorTest");

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  //
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parameters.parallel, "-p", "--parallel", "", "--serial",
                 "enable or disable parallel calulation", false);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  if (parameters.mesh_file == nullptr) {
    mfem_mgis::getOutputStream() << "no mesh file specified\n";
    return EXIT_FAILURE;
  }
  if (parameters.order <= 0) {
    mfem_mgis::getOutputStream() << "invalid finite element order\n";
    return EXIT_FAILURE;
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (::mfem_mgis::getMPIrank() == 0) {
    m.addXMLTestOutput(parameters.parallel
                           ? "ParallelGridFunctionInterpolatorTest.xml"
                           : "GridFunctionInterpolatorTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
