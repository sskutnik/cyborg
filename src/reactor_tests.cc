#include "reactor_tests.h"  
#include "facility_tests.h"
#include "agent_tests.h"

using cyborg::reactor;
using pyne::nucname::id;
using cyclus::Cond;
using cyclus::QueryResult;
using cyclus::toolkit::MatQuery;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace cyborg {
namespace ReactorTests {

Composition::Ptr c_uox() {
  cyclus::CompMap m;
  m[id("u235")] = 0.04;
  m[id("u238")] = 0.96;
  return Composition::CreateFromMass(m);
};


void ReactorTest::SetUp() {
  src_facility_ = new cyborg::reactor(tc_.get());
  InitParameters();
  //SetUpReactor();
}

void ReactorTest::TearDown() {
  delete src_facility_;
}

void ReactorTest::InitParameters(){
  in_r1 = std::vector<std::string>(1,"in_r1");
  in_c1 = std::vector<std::string>(1,"in_c1");
  out_c1 = "out_c1";
  fuel_type = "UOX";
  power_cap = 800.0; // MWt
  cycle_time = 12; // months
  refuel_time = 1; // months
  reactor_lifetime = 480;
  enrichment = 4.0;
  //mod_density = 0.0; // Setting to 0 to test auto-interpolation of density
  n_assem_core = 75;
  n_assem_spent = 0;
  n_assem_batch = 25;
  assem_mass = 325.0; // kg 
 
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
  src_facility_->refuel_time = refuel_time;
  src_facility_->reactor_lifetime = reactor_lifetime;
  src_facility_->assem_size = assem_mass;
  src_facility_->core.capacity(assem_mass*n_assem_core);
  src_facility_->spent.capacity(3.0*assem_mass*n_assem_core);
  src_facility_->fresh.capacity(2.0*assem_mass*n_assem_core);
  src_facility_->assembly_type = "w17x17";
  src_facility_->fuel_type = "UOX";
  src_facility_->n_assem_batch = n_assem_batch;
  src_facility_->n_assem_core = n_assem_core;
  src_facility_->core_power_frac = std::vector<double>(n_assem_core/n_assem_batch, 
                                      static_cast<double>(n_assem_batch) / static_cast<double>(n_assem_core));
   
  // Create an input material buffer of fresh fuel (255 MTU), i.e., 5 * core capacity
  cyclus::Composition::Ptr rec = tc_.get()->GetRecipe(in_r1[0]);
  cyclus::Material::Ptr recmat = cyclus::Material::CreateUntracked(src_facility_->core.space()*5.0, rec); 

  for(size_t i=0; i < n_assem_core*2; ++i) { 
     if(src_facility_->fresh.space() > assem_mass) src_facility_->fresh.Push(recmat->ExtractQty(assem_mass));  
  }
}

void ReactorTest::set_cycle_time(const int t) { src_facility_->cycle_time = t; }
void ReactorTest::set_cycle_step(const int t) { src_facility_->cycle_step = t; }
void ReactorTest::set_discharged(const bool d) { src_facility_->discharged = d; }
int  ReactorTest::get_cycle_time() { return src_facility_->cycle_time; }

void ReactorTest::TestInitState(cyborg::reactor* fac){
  EXPECT_EQ(in_r1[0], fac->fuel_recipes[0]) << "Reactor fuel input recipe not set correclty!\n";
  EXPECT_EQ(in_c1[0], fac->fuel_incommods[0]) << "Reactor input fuel commodity not set correctly!\n";
  EXPECT_EQ(out_c1, fac->spent_fuel) << "Reactor spent fuel commodity not set correctly!\n";
  EXPECT_EQ(power_cap, fac->power_cap) << "Reactor capacity not set correctly!\n";
  EXPECT_EQ(cycle_time, fac->cycle_time) << "Facility lifetime not set correctly!\n";
  EXPECT_EQ(reactor_lifetime, fac->reactor_lifetime) << "Reactor lifetime not set correctly!\n";
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility_->str());
  // Test reactor specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Tick) {
  // Directly initialize reactor
  SetUpReactor();
  TestInitState(src_facility_);
 
  // Need to manually load first core
  src_facility_->Load_();
  set_discharged(false);

  // Loop through a full reactor cycle
  for(size_t i=0; i < get_cycle_time() + 2; ++i) {
    set_cycle_step(i); 
    //EXPECT_NO_THROW(src_facility_->Tick()); 
    src_facility_->Tick(); 
  }
  // Test reactor-specific behaviors of the Tick function here...

}

//\TODO Add a TickDecom test

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorTest, Tock) {
  // Directly initialize reactor
  SetUpReactor();
  TestInitState(src_facility_);

  // Need to manually load first core
  src_facility_->Load_();
  set_discharged(false);

  // Loop through a full reactor cycle
  for(size_t i=0; i < get_cycle_time() + 2; ++i) {
    set_cycle_step(i); 
    EXPECT_NO_THROW(src_facility_->Tock()); 
  }
  // Test reactor-specific behaviors of the Tick function here...
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Reproduces cycamore::reactor test to verify correct number of assemblies 
// are popped from core each cycle
TEST(ReactorXMLTests, BatchSizes) {
  std::string config =
     "  <fuel_incommods> <val>LEU</val>  </fuel_incommods>  "
     "  <fuel_recipes>  <val>uox</val>  </fuel_recipes>  "
     "  <fuel_type>UOX</fuel_type> "
     "  <spent_fuel> used fuel </spent_fuel>  "
     ""
     "  <power_cap>400.0 </power_cap> "
     "  <power_name> Electric LOOOOOOVE </power_name> " 
     "  <cycle_time>2</cycle_time>  "
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
  // 6 for initial core, 3 per every other time step for each new batch for remainder
  EXPECT_EQ(6+3*(simdur/2+1), qr.rows.size());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test for proper retirement of cyborg::reactor at end of lifetime
// Based on cycamore::reactor test with minor modifications
TEST(ReactorXMLTests, Retire) {

   int simDur = 50;
   int rxLife = 36;
   int n_assem_core = 3;
   int n_assem_batch = 1;
   int cycleTime = 7;
   int refuelTime = 0;
   double assemSize = 300; // kg

   std::stringstream simInput;
   simInput << "  <fuel_recipes> <val>fresh_uox</val> </fuel_recipes>  "
            << "  <fuel_incommods> <val>LEU</val> </fuel_incommods>   "
            << "  <spent_fuel>UNF</spent_fuel>  "
            << "  <fuel_type>UOX</fuel_type> "
            << "  <cycle_time>" << cycleTime << "</cycle_time>  "
            << "  <refuel_time>" << refuelTime << "</refuel_time>  "
            << "  <assem_size>" << assemSize << "</assem_size>  "
            << "  <n_assem_fresh>" << n_assem_batch << "</n_assem_fresh>  "
            << "  <n_assem_core>" << n_assem_core << "</n_assem_core>  "
            << "  <n_assem_batch>" << n_assem_batch << "</n_assem_batch>  "
            << "  <power_cap>30</power_cap>  ";
 
   cyclus::MockSim sim(cyclus::AgentSpec(":cyborg:reactor"), simInput.str(), simDur, rxLife);
   sim.AddSource("LEU").Finalize();
   sim.AddSink("UNF").Finalize();
   sim.AddRecipe("fresh_uox", c_uox());
   int id = sim.Run();
   
   // Ensure that reactor stops requesting fresh fuel as it approaches retirement
   int num_assem_recv = 
      static_cast<int>(ceil(static_cast<double>(rxLife) / static_cast<double>(cycleTime)) 
      + (n_assem_core - n_assem_batch));

   std::vector<Cond> conds;
   conds.push_back(Cond("ReceiverId", "==", id));
   QueryResult qr = sim.db().Query("Transactions", &conds);
   EXPECT_EQ(num_assem_recv, qr.rows.size()) << "Failed to stop ordering near retirement.";
   
   // reactor should discharge all fuel before/by retirement  
   // Due to discrete material handling, assemblies in == assemblies out, 
   // (i.e., sell transactions == buy transactions)
   conds.clear();
   conds.push_back(Cond("SenderId", "==", id));
   qr = sim.db().Query("Transactions", &conds);
   EXPECT_EQ(num_assem_recv, qr.rows.size()) 
      << "Failed to discharge all material by retirement time";

   qr = sim.db().Query("Transactions", NULL);
/*
   for(int i=0; i < qr.rows.size(); ++i) {
     std::cerr << qr.GetVal<int>("SenderId", i) << "->" << qr.GetVal<int>("ReceiverId", i) << std::endl;
   }
*/
  // Check that the reactor records the power entry on the time step it retires if operating
   int time_online = rxLife / (cycleTime + refuelTime) * cycleTime 
                       + std::min(rxLife % (cycleTime + refuelTime), cycleTime);
   conds.clear();
   conds.push_back(Cond("AgentId", "==", id));
   conds.push_back(Cond("Value", ">", 0));
   qr = sim.db().Query("TimeSeriesPower", &conds);
   EXPECT_EQ(time_online, qr.rows.size())
       << "failed to generate power for the correct number of time steps";
}
} // namespace ReactorTests
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
