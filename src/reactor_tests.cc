#include "reactor_tests.h"  
#include "facility_tests.h"
#include "agent_tests.h"

using cyborg::reactor;
using pyne::nucname::id;
using cyclus::QueryResult;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace cyborg {

Composition::Ptr c_uox() {
  cyclus::CompMap m;
  m[id("u235")] = 0.04;
  m[id("u238")] = 0.96;
  return Composition::CreateFromMass(m);
};


void ReactorTest::SetUp() {
  src_facility_ = new cyborg::reactor(tc_.get());
  InitParameters();
  SetUpReactor();
}

void ReactorTest::TearDown() {
  delete src_facility_;
}

void ReactorTest::InitParameters(){
  in_r1 = std::vector<std::string>(1,"in_r1");
  in_c1 = std::vector<std::string>(1,"in_c1");
  out_c1 = "out_c1";
  power_cap = 800.0; // MWt
  cycle_time = 12; // months
  reactor_lifetime = 480;
  enrichment = 4.0;
  //mod_density = 0.0; // Setting to 0 to test auto-interpolation of density
  n_assem_core = 75;
  n_assem_spent = 0;
  n_assem_batch = 25;
  assem_mass = 325.0; // kg 
  refresh_recipe = true;
 
  cyclus::CompMap v;
  v[922350000] = (enrichment/100.0);
  v[922380000] = (100.0 - enrichment)/100.0;
  //v[80160000]  = 2.0;
  cyclus::Composition::Ptr recipe = cyclus::Composition::CreateFromMass(v);
  tc_.get()->AddRecipe(in_r1[0], recipe);
}

void ReactorTest::SetUpReactor(){
  src_facility_->fuel_recipes = in_r1;
  src_facility_->fuel_incommods = in_c1;
  src_facility_->spent_fuel = out_c1;
  src_facility_->power_cap = power_cap;
  src_facility_->cycle_time = cycle_time;
  src_facility_->reactor_lifetime = reactor_lifetime;
  //src_facility_->enrichment = enrichment;
  //src_facility_->mod_density = mod_density;
  src_facility_->assem_size = assem_mass;
  src_facility_->core.capacity(assem_mass*n_assem_core);
  src_facility_->spent.capacity(3.0*assem_mass*n_assem_core);
  src_facility_->fresh.capacity(2.0*assem_mass*n_assem_core);
  src_facility_->assembly_type = "w17x17";
  src_facility_->n_assem_batch = n_assem_batch;
  src_facility_->n_assem_core = n_assem_core;
  src_facility_->refreshSpentRecipe = refresh_recipe;
  //src_facility_->fuel.capacity(src_facility_->fuel_capacity);
   
  // Create an input material buffer of fresh fuel (255 MTU), i.e., 5 * core capacity
  cyclus::Composition::Ptr rec = tc_.get()->GetRecipe(in_r1[0]);
  cyclus::Material::Ptr recmat = cyclus::Material::CreateUntracked(src_facility_->core.space()*5.0, rec); 

  for(size_t i=0; i < n_assem_core*2; ++i) { 
     if(src_facility_->fresh.space() > assem_mass) src_facility_->fresh.Push(recmat->ExtractQty(assem_mass));  
  }
}

void ReactorTest::set_cycle_time(const int t) { src_facility_->cycle_time = t; }
void ReactorTest::set_cycle_step(const int t) { src_facility_->cycle_step = t; }

void ReactorTest::TestInitState(cyborg::reactor* fac){
  EXPECT_EQ(in_r1[0], fac->fuel_recipes[0]);
  EXPECT_EQ(in_c1[0], fac->fuel_incommods[0]);
  EXPECT_EQ(out_c1, fac->spent_fuel);
  EXPECT_EQ(power_cap, fac->power_cap);
  //EXPECT_EQ(fuel_capacity, fac->fuel_capacity);
  EXPECT_EQ(cycle_time, fac->cycle_time);
  //EXPECT_EQ(cap_factor, fac->cap_factor);
  EXPECT_EQ(reactor_lifetime, fac->reactor_lifetime);
  //EXPECT_EQ(enrichment, fac->enrichment);
  //EXPECT_EQ(mod_density, fac->mod_density);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, InitialState) {
  std::string config =
     "  <fuel_incommods> <val>LEU</val>  </fuel_incommods>  "
     "  <fuel_recipes>  <val>uox</val>  </fuel_recipes>  "
     "  <fuel_type>UOX</fuel_type> "
     "  <spent_fuel> used fuel </spent_fuel>  "
     ""
     "  <power_cap>400.0 </power_cap> "
     "  <power_name> Electric LOOOOOOVE </power_name> " 
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>166.7</assem_size>  "
     "  <n_assem_fresh>6</n_assem_fresh>"
     "  <n_assem_core>6</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  "
     "  <tags> "
     "    <item> <tag>Moderator Density</tag> <value>0.723</value> </item>"
     "  </tags>";
 
  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cyborg:reactor"), config, simdur);
  sim.AddSource("LEU").Finalize();
  sim.AddRecipe("uox", c_uox());

  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  // 7 for initial core, 3 per time step for each new batch for remainder
  EXPECT_EQ(6+3*(simdur-1), qr.rows.size());
  //

  // Test things about the initial state of the facility here
  TestInitState(src_facility_);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility_->str());
  // Test reactor specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Tick) {
  // Need to manually load first core
  src_facility_->Load_();

  // Loop through a full reactor cycle
  for(size_t i=0; i < src_facility_->get_cycle_time() + 2; ++i) {
    //std::cerr << "Setting cycle time to: " << i << "  limit = " << src_facility_->get_cycle_time() + 2 << std::endl; 
    set_cycle_step(i); 
    //EXPECT_NO_THROW(src_facility_->Tick()); 
    src_facility_->Tick(); 
  }
  // Test reactor-specific behaviors of the Tick function here

}

//\TODO Add a TickDecom test

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Tock) {
  //std::cerr << "Reactor assembly type: " << src_facility_->assembly_type << std::endl;
  //std::cerr << "Tock: moderator density = " << src_facility_->mod_density << std::endl;
  //std::cerr << "lib data path: " << src_facility_->lib_path << std::endl;  

 /*
  src_facility_->reactor_time = src_facility_->cycle_time - 1; 
  try{ 
     src_facility_->Tock(); // no depletion
     src_facility_->Tock(); // depletion
     src_facility_->Tock(); // no depletion

  }
  catch( std::exception& ex) {
    std::cerr << "Exception thrown: " << ex.what() << std::endl;
  }
  catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }  
*/

  // Test for one step before, during, and after a depletion step
  set_cycle_time(src_facility_->get_cycle_time() - 1); 
  EXPECT_NO_THROW(src_facility_->Tock()); // no depletion
  EXPECT_NO_THROW(src_facility_->Tock()); // depletion
  EXPECT_NO_THROW(src_facility_->Tock()); // no depletion

  // Test reactor specific behaviors of the Tock function here
}
} // namespace cyborg

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
