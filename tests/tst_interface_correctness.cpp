/*!
 *  \file - tst_interface_correctness.cpp
 *  \author - Nicholas C. Sly, Steven E. Skutnik
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include <gtest/gtest.h>

#include "cyclus_origen_interface.h"
#include "orglib_default_location.h"
#include "error.h"
#include "Origen/Core/dc/ConcentrationConverter.h"
#include "Origen/Core/dc/ConcentrationUnit.h"
#include "Origen/Core/dc/FakeFactory.h"
#include "Origen/Core/dc/Library.h"
#include "Origen/Core/dc/TagManager.h"

using namespace OrigenInterface;

// iso: isotope (generally zzzaaa) and val: value in GATOMS.
double unitConv(int iso, double val)
{
  auto conv = Origen::ConcentrationConverter();
  conv.convert_to(Origen::ConcentrationUnit::KILOGRAMS,iso,Origen::ConcentrationUnit::GATOMS,val);
  return val;
}

class OrigenResultTester : public ::testing::Test {
protected:

  OrigenResultTester() {
  } 

  void SetUp() {
      // Setting the path to the libraries
      lib_path = ORIGEN_LIBS_DEFAULT;

      // Setting ID tags for w17x17 assembly
      id_tags["Assembly Type"] = "w17x17";

      // Setting Interp tags for enrichment
      interp_tags["Enrichment"] = 2.7;

      // Setting material IDS
      ids.push_back(922340);
      ids.push_back(922350);
      ids.push_back(922360);
      ids.push_back(922380);

      // Setting material initial concentrations in kg
      concs.push_back(0.534); // for u-234
      concs.push_back(27.); // for u-235
      concs.push_back(0.276); // for u-236
      concs.push_back(972.19); // for u-238

      // Setting powers and breaks in Watts.
      powers.insert(powers.end(), 10, 20.0e6);
      powers.insert(powers.end(), 9, 0.0);
      powers.insert(powers.end(), 10, 19.0e6);
      powers.insert(powers.end(), 9, 0.0);
      powers.insert(powers.end(), 10, 18.0e6);
      powers.insert(powers.end(), 9, 0.0);

      //  Setting times (must be cumulative).. in days      
      // Irradiation cycle #1
      for( unsigned int i=0; i < 11; ++i) { times.push_back(540.0/10 * i); }
      // Decay #1
      for( unsigned int i=9; i > 0; --i) { times.push_back(540 + 40.0/pow(3,(i-1))); }

      // Irradiation #2
      for( unsigned int i=1; i < 11; ++i) { times.push_back(580 + 540.0/10 * i); }
      // Decay #2 (40 days)
      for( unsigned int i=9; i > 0; --i) { times.push_back(1120.0 + 40.0/pow(3,(i-1))); }

      // Irradiation #3
      for( unsigned int i=1; i < 11; ++i) { times.push_back(1160 + 540.0/10 * i); }
      // Decay #3 (final)
      for( unsigned int i=9; i > 0; --i) { times.push_back(1700.0 + 1825/pow(3,(i-1))); }
      
       // std::cerr << "Pushed back: " << 1700.0 + 1825.0/pow(3,(i-1)) << std::endl; }

      std::cout << "powers.size() " << powers.size() << "  times.size() = " << times.size() << std::endl;
  }

  void TearDown() {    
  }

  std::string lib_path;
  std::map<std::string,std::string> id_tags;
  std::map<std::string,double> interp_tags;
  std::vector<int> ids;
  std::vector<double> powers, times, concs;
  cyclus2origen tester;
};


TEST_F(OrigenResultTester,CorrectnessTest){
  tester.set_lib_path(lib_path); 
  tester.set_id_tags(id_tags);
  tester.set_parameters(interp_tags);
  EXPECT_NO_THROW(tester.interpolate());
  tester.set_powers(powers);
  tester.set_time_steps(times);
  EXPECT_NO_THROW(tester.set_materials(ids,concs));
  EXPECT_NO_THROW(tester.solve());

  std::map<int,double> masses_out;

  tester.get_masses_at_easy(0,masses_out);

  double test234 = 0.534;

  auto tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Initial mass for u-234 does not match previous external calculation results.";

  // Value 445.7490 grams comes from an externally executed combination of obiwan and origen depletion.
  test234 = 0.445749;

  masses_out.clear();
  tester.get_masses_at_easy(10,masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Mass of u-234 after initial burn does not match previous external calculation results.";

  test234 = 0.3710746;

  masses_out.clear();
  tester.get_masses_at_easy(29,masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Mass of u-234 after second burn does not match previous external calculation results.";

  test234 = 0.3711255;

  masses_out.clear();
  tester.get_masses_at_easy(38,masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Mass of u-234 after second decay does not match previous external calculation results.";

  test234 = 0.3085457;

  masses_out.clear();
  tester.get_masses_at_easy(48,masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Mass of u-234 after final burn does not match previous external calculation results.";

  test234 = 0.3146113;

  masses_out.clear();
  tester.get_masses_at_easy(56,masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Mass of u-234 after final burn does not match previous external calculation results.";

  test234 = 0.3146113;

  masses_out.clear();
  tester.get_masses_final_easy(masses_out);

  tmp234 = masses_out[922340];

  EXPECT_NEAR(tmp234,test234,0.001) << "Resulting mass for u-234 does not match previous external calculation results.";
}
