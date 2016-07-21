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
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(20.e6);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(19.e6);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(18.e6);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);
      powers.push_back(0.);

      //  Setting times (must be cumulative).. in days
      times.push_back(0.);
      times.push_back(54.);
      times.push_back(108.);
      times.push_back(162.);
      times.push_back(216.);
      times.push_back(270.);
      times.push_back(324.);
      times.push_back(378.);
      times.push_back(432.);
      times.push_back(486.);
      times.push_back(540.);
      times.push_back(540.01);
      times.push_back(540.03);
      times.push_back(540.1);
      times.push_back(540.3);
      times.push_back(541.);
      times.push_back(543.);
      times.push_back(550.);
      times.push_back(570.);
      times.push_back(580.);
      times.push_back(634.);
      times.push_back(688.);
      times.push_back(742.);
      times.push_back(796.);
      times.push_back(850.);
      times.push_back(904.);
      times.push_back(958.);
      times.push_back(1012.);
      times.push_back(1066.);
      times.push_back(1120.);
      times.push_back(1120.01);
      times.push_back(1120.03);
      times.push_back(1120.1);
      times.push_back(1120.3);
      times.push_back(1121.);
      times.push_back(1123.);
      times.push_back(1130.);
      times.push_back(1150.);
      times.push_back(1160.);
      times.push_back(1214.);
      times.push_back(1268.);
      times.push_back(1322.);
      times.push_back(1376.);
      times.push_back(1430.);
      times.push_back(1484.);
      times.push_back(1538.);
      times.push_back(1592.);
      times.push_back(1646.);
      times.push_back(1700.);
      times.push_back(1701.);
      times.push_back(1703.);
      times.push_back(1710.);
      times.push_back(1730.);
      times.push_back(1800.);
      times.push_back(2000.);
      times.push_back(2700.);
      times.push_back(3525.);
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
