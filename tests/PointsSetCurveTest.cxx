/*!
 * \file   tests/PointsSetCurvesWriter.cxx
 * \brief
 * \author Thomas Helfer
 * \date   09/09/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/GridFunctionUtilities.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx"
#include "UnitTestingUtilities.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int linearsolver = 0;
  int order = 1;
  int parallel = 0;
};  // end of struct TestParameters

auto parameters = TestParameters{};

static void parseCommandLineOptions(int argc, char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.linearsolver, "-ls", "--linearsolver",
                 "identifier of the linear solver: 0 -> CG, 1 -> GMRES, 2 -> "
                 "UMFPack (serial), 3-> MUMPS(serial), 2 -> HypreFGMRES "
                 "(//), 3 -> HyprePCG (//), 4 -> HypreGMRES (//)");
  args.AddOption(&parameters.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parameters.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (parameters.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

struct PointsSetCurveTest final : public tfel::tests::TestCase {
  PointsSetCurveTest()
      : tfel::tests::TestCase("MFEMMGIS", "PointsSetCurveTest") {
  }  // end of PointsSetCurveTest
  tfel::tests::TestResult execute() override {
    if (parameters.parallel) {
#ifdef MFEM_USE_MPI
      this->test1<true>();
#endif /* MFEM_USE_MPI */
    } else {
      this->test1<false>();
    }
    return this->result;
  }

 private:
  //
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto lOA =
        dict{{"Line", dict{{"InitialPoint", "O"},
                           {"FinalPoint", "A"},
                           {"Discretization",
                            dict{{"Uniform", dict{{"NumberOfPoints", 20}}}}}}}};
    auto lCB =
        dict{{"Line", dict{{"InitialPoint", "C"},
                           {"FinalPoint", "B"},
                           {"Discretization",
                            dict{{"Uniform", dict{{"NumberOfPoints", 20}}}}}}}};
    auto ofed = construct<FiniteElementDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", parameters.order},
                  {"UnknownsSize", 2},
                  {"NumberOfUniformRefinements", 0},
                  {"Parallel", parameters.parallel == 1},
                  {"Points", dict{{"O", list{0, 0}},
                                  {"A", list{1, 0}},
                                  {"B", list{1, 0.2}},
                                  {"C", list{0, 0.2}}}},
                  {"PointsSets", dict{{"lOA", lOA},  //
                                      {"lCB", lCB}}}});
    TFEL_TESTS_ASSERT(isValid(ofed));
    //
    auto ospace = make_shared<PartialQuadratureSpace>(
        ctx, *ofed, 5,
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation& tr) noexcept
        -> const mfem::IntegrationRule& {
          const auto order = 2 * tr.OrderGrad(&e);
          return mfem::IntRules.Get(e.GetGeomType(), order);
        });
    TFEL_TESTS_ASSERT(isValid(ospace));
    auto fct = PartialQuadratureFunction::evaluate(
        ospace, [](const real x, const real y) noexcept { return cos(x) * y; });
    auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
        ctx, ofed->getFiniteElementSpace<parallel>(), parameters);
    TFEL_TESTS_ASSERT(isValid(s));
    const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
    TFEL_TESTS_ASSERT(isValid(oresult));
    //
    auto owriter = construct<PointsSetCurvesWriter>(
        ctx, *ofed,
        dict{{"File", parameters.parallel == 1
                          ? "ParallelPointsSetCurveTestResults.txt"
                          : "PointsSetCurveTestResults.txt"},
             {"PointsSet", "lCB"},
             {"Precision", 14},
             {"ExportCurvilinearAbscissa", true},
             {"ExportCoordinates", true}});
    TFEL_TESTS_ASSERT(isValid(owriter));
    TFEL_TESTS_ASSERT(owriter->add(ctx, "L2Projection", *(oresult->result)));
    TFEL_TESTS_ASSERT(owriter->writeFileHeader(ctx));
    TFEL_TESTS_ASSERT(
        owriter->writeValues(ctx, {.begin = 0, .end = 0, .dt = 0}, ets));
  }  // end of test1
};

TFEL_TESTS_GENERATE_PROXY(PointsSetCurveTest, "PointsSetCurveTest");

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // options treatment
  initialize(argc, argv);
  parseCommandLineOptions(argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (::mfem_mgis::getMPIrank() == 0) {
    m.addXMLTestOutput(parameters.parallel ? "ParallelPointsSetCurveTest.xml"
                                           : "PointsSetCurveTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
