/*!
 * \file   GridFunctionTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   26/01/2026
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
  // creation of the grid function
  auto f = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
  auto c = mfem::FunctionCoefficient(
      [](const mfem::Vector& x) { return cos(x[0]) * exp(-2 * x[1]); });
  f.ProjectCoefficient(c);
  //
  auto qspace = std::make_shared<PartialQuadratureSpace>(
      fed, 5,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto qf = PartialQuadratureFunction(qspace);
  if (!update(ctx, qf, f)) {
    return false;
  }
  // some test
  auto qf2 = PartialQuadratureFunction::evaluate(
      qspace,
      [](const real x, const real y) noexcept { return cos(x) * exp(-2 * y); });
  for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
    const auto v1 = qf.getIntegrationPointValue(i);
    const auto v2 = qf2->getIntegrationPointValue(i);
    if (std::abs(v1 - v2) > 2e-2) {
      ctx.registerErrorMessage("invalid value at integration point '" +
                               std::to_string(i) + "' (computed " +
                               std::to_string(v1) + ", expected " +
                               std::to_string(v2) + ", error " +
                               std::to_string(std::abs(v1 - v2)) + ")");
      return false;
    }
  }
  // back on nodes
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {qf});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("GridFunctionResult-test1");
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
  //
  return true;
}

// multi-component test
template <bool parallel>
bool test2(mfem_mgis::Context& ctx, const TestParameters& params) {
  using namespace mfem_mgis;
  auto fed = FiniteElementDiscretization{
      {{"MeshFileName", params.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", params.order},
       {"UnknownsSize", 2},
       {"NumberOfUniformRefinements", parallel ? 1 : 0},
       {"Parallel", parallel}}};
  // creation of the grid function
  auto f = GridFunction<parallel>{&(fed.getFiniteElementSpace<parallel>())};
  auto c = mfem::VectorFunctionCoefficient(
      2, [](const mfem::Vector& x, mfem::Vector& v) {
        v.SetSize(2);
        v[0] = std::cos(x[0]) * std::exp(-2 * x[1]);
        v[1] = std::exp(x[0] * x[0]) * std::exp(x[1]);
      });
  f.ProjectCoefficient(c);
  //
  auto qspace = std::make_shared<PartialQuadratureSpace>(
      fed, 5,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
      -> const mfem::IntegrationRule& {
        const auto order = 2 * tr.OrderGrad(&e);
        return mfem::IntRules.Get(e.GetGeomType(), order);
      });
  auto qf = PartialQuadratureFunction(qspace, 2);
  if (!update(ctx, qf, f)) {
    return false;
  }
  // some test
  auto qf2 = PartialQuadratureFunction::evaluate(
      qspace, [](const real x, const real y) noexcept {
        return std::cos(x) * std::exp(-2 * y);
      });
  auto qf3 = PartialQuadratureFunction::evaluate(
      qspace, [](const real x, const real y) noexcept {
        return std::exp(x * x) * std::exp(y);
      });
  for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
    const auto v1 = qf.getIntegrationPointValues(i);
    const auto v2 = qf2->getIntegrationPointValue(i);
    const auto v3 = qf3->getIntegrationPointValue(i);
    if (std::max(std::abs(v1[0] - v2), std::abs(v1[1] - v3)) > 2e-2) {
      ctx.registerErrorMessage(
          "invalid value at integration point '" + std::to_string(i) +
          "' (computed (" + std::to_string(v1[0]) + ", " +
          std::to_string(v1[1]) + "), expected (" + std::to_string(v2) + +", " +
          std::to_string(v3) + "), error " +
          std::to_string(std::max(std::abs(v1[0] - v2), std::abs(v1[1] - v3))) +
          ")");
      return false;
    }
  }
  // back on nodes
  auto s = mfem_mgis::unit_tests::getLinearSolver<parallel>(
      ctx, fed.getFiniteElementSpace<parallel>(), params);
  if (isInvalid(s)) {
    return false;
  }
  const auto oresult = computeL2Projection<parallel>(ctx, s, {qf});
  if (isInvalid(oresult)) {
    return false;
  }
  mfem::ParaViewDataCollection exporter("GridFunctionResult-test2");
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
}

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
