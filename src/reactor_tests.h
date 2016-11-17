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
  cyborg::Reactor* src_facility_;

  virtual void SetUp();
  virtual void TearDown();
  void InitParameters();
  void SetUpReactor();

  void set_cycle_time(const int);
  void set_cycle_step(const int);
  void set_discharged(const bool);
  int get_cycle_time();
  void TestInitState(cyborg::Reactor* fac);

  std::vector<std::string> in_r1, in_c1;
  std::string out_c1;
  std::string fuel_type;
  double power_cap, core_capacity, enrichment, mod_density;
  int cycle_time, refuel_time;
  int n_assem_core, n_assem_spent, n_assem_batch; 
  double assem_mass;
};

} //namespace ReactorTests
} // namespace cyborg

#endif // REACTOR_TESTS_H_
