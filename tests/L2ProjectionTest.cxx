/*!
 * \file   tests/L2ProjectionTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   20/01/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "UnitTestingUtilities.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int linearsolver = 0;
  int order = 1;
  int parallel = 0;
};  // end of struct TestParameters

static void parseCommandLineOptions(TestParameters& params,
                                    int argc,
                                    char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&params.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&params.linearsolver, "-ls", "--linearsolver",
                 "identifier of the linear solver: 0 -> CG, 1 -> GMRES, 2 -> "
                 "UMFPack (serial), 3-> MUMPS(serial), 2 -> HypreFGMRES "
                 "(//), 3 -> HyprePCG (//), 4 -> HypreGMRES (//)");
  args.AddOption(&params.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&params.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (params.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

// single component test
template <bool parallel>
bool test(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, 5,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto fct = PartialQuadratureFunction::evaluate(
      space, [](const real x, const real y) noexcept { return cos(x) * y; });
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("Result-test1");
  if constexpr (parallel) {
#ifdef MFEM_USE_MPI
    exporter.SetMesh(&(fed.getMesh<true>()));
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  } else {
    exporter.SetMesh(&(fed.getMesh<false>()));
  }
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("Result", oresult->result.get());
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
  return true;
}  // end of test

// multi-component test
template <bool parallel>
bool test2(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  auto space = std::make_shared<PartialQuadratureSpace>(
      fed, 5,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto f1 = PartialQuadratureFunction::evaluate(
      space, [](const real x, const real y) noexcept { return cos(x) * y; });
  auto f2 = PartialQuadratureFunction::evaluate(
      space,
      [](const real x, const real y) noexcept { return exp(-x) * y * y; });
  auto fct = std::make_shared<PartialQuadratureFunction>(space, 2);
  auto f1_values = f1->getValues();
  auto f2_values = f2->getValues();
  auto fct_values = fct->getValues();
  for (size_type idx = 0; idx != space->getNumberOfIntegrationPoints(); ++idx) {
    fct_values[2 * idx] = f1_values[idx];
    fct_values[2 * idx + 1] = f2_values[idx];
  }
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {*fct});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("Result-test2");
  if constexpr (parallel) {
#ifdef MFEM_USE_MPI
    exporter.SetMesh(&(fed.getMesh<true>()));
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  } else {
    exporter.SetMesh(&(fed.getMesh<false>()));
  }
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("Result", oresult->result.get());
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
  return true;
}  // end of test2

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  // options treatment
  auto params = TestParameters{};
  initialize(argc, argv);
  parseCommandLineOptions(params, argc, argv);
  const auto success = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test<false>(ctx, params);
  }();
  const auto success2 = [&ctx, &params] {
    if (params.parallel) {
#ifdef MFEM_USE_MPI
      return test2<true>(ctx, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return test2<false>(ctx, params);
  }();
  if (!(success && success2)) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
