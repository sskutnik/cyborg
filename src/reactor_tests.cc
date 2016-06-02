#include "reactor_tests.h"  
#include "facility_tests.h"
#include "agent_tests.h"

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
  power_cap = 100.0; // MWt
  fuel_capacity = 20.0E3; // kg
  cycle_length = 1;
  cap_factor = 0.9;
  reactor_lifetime = 10;
  enrichment = 4.0;
  mod_density = 0.0; // Setting to 0 to test auto-interpolation of density

  cyclus::CompMap v;
  v[922350000] = (enrichment/100.0);
  v[922380000] = (100.0 - enrichment)/100.0;
  v[80160000]  = 2.0;
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
  src_facility_->mod_density = mod_density;

  src_facility_->fuel.capacity(src_facility_->fuel_capacity);
   
  // Create an input material buffer of fresh fuel (100 MTU)
  cyclus::Composition::Ptr rec = tc_.get()->GetRecipe(in_r1);
  cyclus::Material::Ptr recmat = cyclus::Material::CreateUntracked(src_facility_->fuel.space()*5.0, rec);
  src_facility_->fresh_inventory.Push(recmat); 
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
  EXPECT_EQ(mod_density, fac->mod_density);
  std::cerr << "TEST: fac->mod_density = " << fac->mod_density << std::endl;
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

//\TODO Add a TickDecom test

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(reactorTest, Tock) {
  src_facility_->fuel.capacity(src_facility_->fuel_capacity*1000);
  src_facility_->assembly_type = "w17x17";

  //std::cerr << "Reactor assembly type: " << src_facility_->assembly_type << std::endl;
  //std::cerr << "Tock: moderator density = " << src_facility_->mod_density << std::endl;
  std::cerr << "lib data path: " << src_facility_->lib_path << std::endl; 
  
  try{ src_facility_->Tock(); }
  catch( std::exception& ex) {
    std::cerr << "Exception thrown: " << ex.what() << std::endl;
  }
  catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }  

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
