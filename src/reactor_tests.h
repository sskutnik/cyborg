#ifndef REACTOR_TESTS_H_
#define REACTOR_TESTS_H_

#include <gtest/gtest.h>
#include "reactor.h"
#include "context.h"
#include "test_context.h"

namespace cyborg {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class ReactorTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc_;
  cyborg::reactor* src_facility_;

  virtual void SetUp();
  virtual void TearDown();
  void InitParameters();
  void SetUpReactor();

  void set_cycle_time(const int);
  void TestInitState(cyborg::reactor* fac);

  std::string in_r1, in_c1, out_c1;
  double power_cap, core_capacity, enrichment, mod_density;
  int cycle_length, reactor_lifetime;
  int n_assem_core, n_assem_spent, n_assem_batch; 
  double assem_mass;
};
} // namespace reactor
#endif // REACTOR_TESTS_H_
