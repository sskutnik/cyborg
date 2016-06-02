#ifndef REACTOR_TESTS_H_
#define REACTOR_TESTS_H_

#include <gtest/gtest.h>
#include "reactor.h"
#include "context.h"
#include "test_context.h"

namespace reactor {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class reactorTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc_;
  cyborg::reactor* src_facility_;

  virtual void SetUp();
  virtual void TearDown();
  void InitParameters();
  void SetUpReactor();
  void TestInitState(cyborg::reactor* fac);

  std::string in_r1, in_c1, out_c1;
  double power_cap, fuel_capacity, cap_factor, enrichment, mod_density;
  int cycle_length, reactor_lifetime;
  
};
} // namespace reactor
#endif // REACTOR_TESTS_H_
