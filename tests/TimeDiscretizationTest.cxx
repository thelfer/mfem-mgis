/*!
 * \file   tests/TimeDiscretizationTests.cxx
 * \author Thomas Helfer
 * \date   06/04/2023
 */

#include <algorithm>
#include "MFEMMGIS/Simulation.hxx"

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

struct TimeDiscretizationTest final : public tfel::tests::TestCase {
  TimeDiscretizationTest()
      : tfel::tests::TestCase("MFEMMGIS", "TimeDiscretizationTest") {
  }  // end of TimeDiscretizationTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    return this->result;
  }

 private:
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto td2 = construct<Simulation::TimesDescription>(ctx, 0, -1);
    TFEL_TESTS_CHECK(isInvalid(td2));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("invalid time sequence"));
    auto td3 = construct<Simulation::TimesDescription>(ctx, 0, 0);
    TFEL_TESTS_CHECK(isInvalid(td3));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("time increment too small"));
    auto td4 = construct<Simulation::TimesDescription>(ctx, 1, 2);
    TFEL_TESTS_ASSERT(isValid(td4));
    const auto times = std::vector<real>(td4->cbegin(), td4->cend());
    TFEL_TESTS_ASSERT(times.size() == 2);
    TFEL_TESTS_CHECK(std::abs(times[0] - 1) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[1] - 2) < 1e-14);
  }
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto td1 = construct<Simulation::TimesDescription>(ctx, 0, 1, 0);
    TFEL_TESTS_CHECK(isInvalid(td1));
    TFEL_TESTS_CHECK(ctx.getErrorMessage().starts_with(
        "invalid number of temporal sequences"));
    auto td2 = construct<Simulation::TimesDescription>(ctx, 0, -1, 1);
    TFEL_TESTS_CHECK(isInvalid(td2));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("invalid time sequence"));
    auto td3 = construct<Simulation::TimesDescription>(ctx, 0, 0, 10);
    TFEL_TESTS_CHECK(isInvalid(td3));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("time increment too small"));
    auto td4 = construct<Simulation::TimesDescription>(ctx, 0, 1, 3);
    TFEL_TESTS_ASSERT(isValid(td4));
    const auto times = std::vector<real>(td4->cbegin(), td4->cend());
    TFEL_TESTS_ASSERT(times.size() == 4);
    TFEL_TESTS_CHECK(std::abs(times[0] - 0) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[1] - real{1} / 3) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[2] - real{2} / 3) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[3] - 1) < 1e-14);
  }
  void test3() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto td1 = construct<Simulation::TimesDescription>(ctx, std::vector<real>{0});
    TFEL_TESTS_CHECK(isInvalid(td1));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("invalid number of times"));
    auto td2 =
        construct<Simulation::TimesDescription>(ctx, std::vector<real>{0, -1, 2});
    TFEL_TESTS_CHECK(isInvalid(td2));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("negative time increment"));
    auto td3 =
        construct<Simulation::TimesDescription>(ctx, std::vector<real>{0, 1, 1});
    TFEL_TESTS_CHECK(isInvalid(td3));
    TFEL_TESTS_CHECK(
        ctx.getErrorMessage().starts_with("invalid time increment"));
    auto td4 =
        construct<Simulation::TimesDescription>(ctx, std::vector<real>{0, 1, 2, 5});
    TFEL_TESTS_ASSERT(isValid(td4));
    const auto times = std::vector<real>(td4->cbegin(), td4->cend());
    TFEL_TESTS_ASSERT(times.size() == 4);
    TFEL_TESTS_CHECK(std::abs(times[0] < 0) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[1] < 1) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[2] < 2) < 1e-14);
    TFEL_TESTS_CHECK(std::abs(times[3] < 5) < 1e-14);
  }
};

TFEL_TESTS_GENERATE_PROXY(TimeDiscretizationTest, "TimeDiscretizationTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("TimeDiscretizationTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
