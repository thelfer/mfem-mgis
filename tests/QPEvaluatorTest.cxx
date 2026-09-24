/*!
 * \file   QPEvaluatorTest.cxx
 * \brief
 * \author th202608
 * \date   22/09/2026
 */

#include <array>
#include <numeric>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/QPEvaluator/QPEvaluator.hxx"

struct QPEvaluatorTest final : public tfel::tests::TestCase {
  QPEvaluatorTest()
      : tfel::tests::TestCase("MFEMMGIS", "QPEvaluatorTest") {
  }  // end of QPEvaluatorTest
  tfel::tests::TestResult execute() override {
    this->test1<1>();
    this->test1<2>();
    this->test1<3>();
    this->test1<4>();
    this->test1<5>();
    this->test1<6>();
    this->test1<7>();
    this->test1<8>();
    this->test1<9>();
    this->test1<10>();
    this->test1<11>();
    this->test1<12>();
    return this->result;
  }

 private:
  template <mfem_mgis::size_type N>
  void test1() {
    using namespace mfem_mgis;
    auto in = std::array<size_type, N>{};
    auto out = std::array<size_type, N>{};
    std::iota(in.begin(), in.end(), N);
    algorithm::copy<N>(in.begin(), in.end(), out.begin());
    TFEL_TESTS_CHECK(in == out);
  }
};

TFEL_TESTS_GENERATE_PROXY(QPEvaluatorTest, "QPEvaluatorTest");

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("QPEvaluatorTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
