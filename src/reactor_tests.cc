#include <gtest/gtest.h>

#include "reactor_tests.h"  

using cyborg::reactor;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace reactor {

void reactorTest::SetUp() {
  src_facility_ = new cyborg::reactor(tc_.get());
  InitParameters();
  SetUpReactor();
}

void reactorTest::TearDown() {
  delete src_facility_;
}

void reactorTest::InitParameters(){
  in_r1 = "in_r1";
  in_c1 = "in_c1";
  out_c1 = "out_c1";
  power_cap = 100.0;
  fuel_capacity = 20.0;
  cycle_length = 1;
  cap_factor = 0.9;
  reactor_lifetime = 10;
  enrichment = 4.0;

  cyclus::CompMap v;
  v[922350000] = 0.04;
  v[922380000] = 0.96;
  cyclus::Composition::Ptr recipe = cyclus::Composition::CreateFromAtom(v);
  tc_.get()->AddRecipe(in_r1, recipe);
}

void reactorTest::SetUpReactor(){
  src_facility_->fuel_recipe = in_r1;
  src_facility_->fresh_fuel = in_c1;
  src_facility_->spent_fuel = out_c1;
  src_facility_->power_cap = power_cap;
  src_facility_->fuel_capacity = fuel_capacity;
  src_facility_->cycle_length = cycle_length;
  src_facility_->cap_factor = cap_factor;
  src_facility_->reactor_lifetime = reactor_lifetime;
  src_facility_->enrichment = enrichment;
}

void reactorTest::TestInitState(cyborg::reactor* fac){
  EXPECT_EQ(in_r1, fac->fuel_recipe);
  EXPECT_EQ(in_c1, fac->fresh_fuel);
  EXPECT_EQ(out_c1, fac->spent_fuel);
  EXPECT_EQ(power_cap, fac->power_cap);
  EXPECT_EQ(fuel_capacity, fac->fuel_capacity);
  EXPECT_EQ(cycle_length, fac->cycle_length);
  EXPECT_EQ(cap_factor, fac->cap_factor);
  EXPECT_EQ(reactor_lifetime, fac->reactor_lifetime);
  EXPECT_EQ(enrichment, fac->enrichment);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, InitialState) {
  // Test things about the initial state of the facility here
  TestInitState(src_facility_);
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
