#ifndef REACTOR_TESTS_H_
#define REACTOR_TESTS_H_

#include <gtest/gtest.h>
#include "reactor.h"
#include "context.h"
#include "test_context.h"

namespace cyborg {
namespace ReactorTests {

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
  void set_cycle_step(const int);
  void TestInitState(cyborg::reactor* fac);

  std::vector<std::string> in_r1, in_c1;
  std::string out_c1;
  double power_cap, core_capacity, enrichment, mod_density;
  int cycle_time, reactor_lifetime;
  int n_assem_core, n_assem_spent, n_assem_batch; 
  double assem_mass;
  bool refresh_recipe;
};

} //namespace ReactorTests
} // namespace cyborg

#endif // REACTOR_TESTS_H_
