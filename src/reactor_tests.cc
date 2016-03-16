#include <gtest/gtest.h>

#include "reactor_tests.h"  

using cyborg::reactor;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace reactor {

void reactorTest::SetUp() {
  src_facility_ = new cyborg::reactor(tc_.get());
  // InitParameters();
  // SetUpReactor();
}

void reactorTest::TearDown() {
  delete src_facility_;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility_->str());
  // Test reactor specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, Tick) {
  ASSERT_NO_THROW(src_facility_->Tick());
  // Test reactor specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, Tock) {
  EXPECT_NO_THROW(src_facility_->Tock());
  // Test reactor specific behaviors of the Tock function here
}
} // namespace reactor

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* reactorConstructor(cyclus::Context* ctx) {
  return new cyborg::reactor(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(reactor, FacilityTests,
                        ::testing::Values(&reactorConstructor));
INSTANTIATE_TEST_CASE_P(reactor, AgentTests,
                        ::testing::Values(&reactorConstructor));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
