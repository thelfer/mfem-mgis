/*!
 * \file   tests/ParametersValidatorTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   19/09/2026
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Parameters.hxx"

struct ParametersValidatorTest final : public tfel::tests::TestCase {
  ParametersValidatorTest()
      : tfel::tests::TestCase("MFEMMGIS", "ParametersValidatorTest") {
  }  // end of ParametersValidatorTest

  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    return this->result;
  }  // end of execute

 private:
  //
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto validator = ParametersValidator{}.add<std::string>("a");
    //
    auto d1 = Parameters{};
    // a is not required, so this is ok
    TFEL_TESTS_CHECK(validator.validate(ctx, d1));
    // a has not the good type
    auto d2 = Parameters{{"a", 12}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d2));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "parameter 'a' does not hold the expected type");
    // b is not authorized, only a is
    auto d3 = Parameters{{"b", 12}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d3));  // b is not authorized
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "invalid parameter 'b'. Valid parameters are:\n- 'a' (undocumented)");
  }
  void test2() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto validator =
        ParametersValidator{}.addKeysIncompatibilityCheck({"a", "b"});
    //
    // a and b are not required, so this is ok
    auto d1 = Parameters{};
    TFEL_TESTS_CHECK(validator.validate(ctx, d1));
    // a is authorized, so this is ok
    auto d2 = Parameters{{"a", 12}};
    TFEL_TESTS_CHECK(validator.validate(ctx, d2));
    // b is authorized, so this is ok
    auto d3 = Parameters{{"b", 12}};
    TFEL_TESTS_CHECK(validator.validate(ctx, d3));
    //
    auto d4 = Parameters{{"a", "t"}, {"b", 12}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d4));  // b is not authorized
    const auto e = ctx.getRawErrorMessage();
    TFEL_TESTS_CHECK(
        (e ==
         "parameters 'a' and 'b' are exclusive: only one shall be defined") ||
        (e ==
         "parameters 'b' and 'a' are exclusive: only one shall be defined"));
  }
  void test3() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto validator = ParametersValidator{}.addStrictlyPositiveIntegerCheck(
        "a", {.required = true});
    // a is required
    auto d1 = Parameters{};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d1));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "required parameter 'a' is missing");
    // a is not an integer
    auto d2 = Parameters{{"a", "b"}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d2));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "parameter 'a' is not an integer");
    // a is a strictly positive integer
    auto d3 = Parameters{{"a", 12}};
    TFEL_TESTS_CHECK(validator.validate(ctx, d3));
    // a is a negative integer
    auto d4 = Parameters{{"a", -12}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d4));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "parameter 'a' is not strictly postive");
    // a is a null
    auto d5 = Parameters{{"a", 0}};
    TFEL_TESTS_CHECK(!validator.validate(ctx, d5));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "parameter 'a' is not strictly postive");
  }
};

TFEL_TESTS_GENERATE_PROXY(ParametersValidatorTest, "ParametersValidatorTest");

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ParametersValidatorTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
