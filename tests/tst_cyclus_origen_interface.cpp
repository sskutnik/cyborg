/*!
 *  \file - tst_cyclus_origen_interface.cpp
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

using namespace OrigenInterface;

class OrigenInterfaceTester : public ::testing::Test {
protected:

  OrigenInterfaceTester() {
  } 

  void SetUp() {
     fluxes.push_back(1.e14);
     fluxes.push_back(2.e14);
     fluxes.push_back(3.e14);

     // Powers in watts.
     powers.push_back(1.e7);
     powers.push_back(2.e7);
     powers.push_back(3.e7);

     // Time assignments in days.
     times.push_back(0.);
     times.push_back(100.);
     times.push_back(200.);
     times.push_back(300.); 

     // Initial concentrations; default units of kg
     concs.push_back(4.0);
     concs.push_back(96.0);

     // IDs can be either zzzaaai or pizzzaaa/sizzzaaa.
     ids.push_back(922350);
     ids.push_back(922380);

     // Can specify libraries directly, or use library path.
     lib_names.push_back("ce14_e20.arplib");
     lib_names.push_back("ce14_e30.arplib");

     lib_path = ORIGEN_LIBS_DEFAULT;

     id_tags["Assembly Type"] = "ce14x14";
     id_tags["Fuel Type"] = "Uranium";
     id_tags["Something"] = "Else";

     params["Enrichment"] = concs[0]/(concs[0] + concs[1]);
     params["Moderator Density"] = 0.45;
     params["Fuel Temperature"] = 9000;
  }

  void TearDown() {    
  }

  std::vector<std::string> lib_names;
  std::string lib_path;
  std::map<std::string,std::string> id_tags;
  std::map<std::string,double> params;
  std::vector<int> ids;
  std::vector<int> iblank;
  std::vector<double> fluxes, powers, times, concs;
  std::vector<double> dblank;
  cyclus2origen tester;
};

TEST_F(OrigenInterfaceTester,libManipulation)
{
  tester.set_lib_names(lib_names);

  std::vector<std::string> get_names;
  tester.get_lib_names(get_names);
  EXPECT_EQ(get_names,lib_names) << "Setting or getting library names failed.";

  tester.add_lib_names({"ce14_e40.arplib"});

  lib_names.push_back("ce14_e40.arplib");

  get_names.clear();
  tester.get_lib_names(get_names);
  EXPECT_EQ(get_names,lib_names) << "Adding a new library name failed.";

  tester.remove_lib_names({"ce14_e40.arplib"});
  lib_names.pop_back();

  get_names.clear();
  tester.get_lib_names(get_names);
  EXPECT_EQ(get_names,lib_names) << "Removing a library name failed.";
}

TEST_F(OrigenInterfaceTester,idTagManipulation){

  EXPECT_THROW(tester.remove_id_tag(id_tags.begin()->first),cyclus::StateError);

  tester.set_id_tags(id_tags);

  std::vector<std::string> names;
  std::vector<std::string> values;

  tester.get_id_tags(names,values);

  for(size_t i = 0; i < names.size(); i++){
    EXPECT_EQ(id_tags[names[i]],values[i]) << "ID Tags not properly emplaced on interface object with map.";
  }

  tester.remove_id_tag("Something");
  tester.remove_id_tag("Fuel Type");

  names.clear();values.clear();

  tester.get_id_tags(names,values);

  EXPECT_EQ(names.size(),1) << "ID Tag removal failed.";
  EXPECT_EQ(id_tags[names[0]],values[0]) << "ID Tag removal resulted in an incorrect list of remaining tags.";

  tester.set_id_tag("Fuel Type","Uranium");

  names.clear(); values.clear();

  tester.get_id_tags(names,values);
  EXPECT_EQ(names.size(),2) << "ID Tag addition failed.";
  for(size_t i = 0; i < names.size(); i++){
    EXPECT_EQ(id_tags[names[i]],values[i]) << "ID Tag addition did not result in a correct list of ID Tags.";
  }
}

TEST_F(OrigenInterfaceTester,parameterManipulation){

  EXPECT_THROW(tester.remove_parameter(params.begin()->first),cyclus::StateError);

  tester.set_parameters(params);

  std::vector<std::string> names;
  std::vector<double> values;

  tester.get_parameters(names,values);
  EXPECT_EQ(names.size(),params.size()) << "Parameter setting or getting failed to return the correct number of parameters.";
  for(size_t i = 0; i < names.size(); i++){
    EXPECT_EQ(params[names[i]],values[i]) << "Parameter setting and getting did not return the correct set of parameters.";
  }

  tester.remove_parameter("Fuel Temperature");
  tester.remove_parameter("Moderator Density");
  names.clear(); values.clear();

  tester.get_parameters(names,values);
  EXPECT_EQ(names.size(),(params.size()-2)) << "Parameter removal did not result in the correct number of remaining parameters.";

  for(size_t i = 0; i < names.size(); i++){
    EXPECT_EQ(params[names[i]],values[i]) << "Parameter removal did not result in a correct set of remaining parameters.";
  }

  tester.add_parameter("Moderator Density",0.45);
  names.clear(); 
  values.clear();

  tester.get_parameters(names,values);
  EXPECT_EQ(names.size(),(params.size()-1)) << "Parameter addition did not result in the correct number of remaining parameters.";

  for(size_t i = 0; i < names.size(); i++){
    EXPECT_EQ(params[names[i]],values[i]) << "Parameter addition did not result in the correct set of remaining parameters.";
  }
}

TEST_F(OrigenInterfaceTester,interpolationTest){
  // Tests for failure, not correctness. Test for correctness in tst_interface_correctness.cpp.

  EXPECT_THROW(tester.interpolate(),cyclus::StateError);

  id_tags.erase("Something");
  tester.set_id_tags(id_tags);

  params.erase("Moderator Density");
  params.erase("Fuel Temperature");
  tester.set_parameters(params);

  EXPECT_THROW(tester.interpolate(),cyclus::ValueError);

  tester.set_lib_path(lib_path);
  std::cout << "Expect next line to be warning about unspecified tag 'Moderator Density'.\n";
  EXPECT_NO_THROW(tester.interpolate());
}

TEST_F(OrigenInterfaceTester,materialTest){

  id_tags.erase("Something");
  params.erase("Fuel Temperature");
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);

  EXPECT_THROW(tester.set_materials(ids,concs),cyclus::StateError);

  EXPECT_NO_THROW(tester.interpolate());

  EXPECT_THROW(tester.set_materials(iblank,concs),cyclus::StateError);
  EXPECT_THROW(tester.set_materials(ids,dblank),cyclus::StateError);

  iblank.push_back(1);
  EXPECT_THROW(tester.set_materials(iblank,concs),cyclus::ValueError);

  tester.set_powers(powers);
  tester.set_time_steps(times);

  EXPECT_NO_THROW(tester.set_materials(ids,concs));
  /* Commented out until such time as fluxes are properly implemented.
   ** Currently hindered by the need to translate to flux in order to
   ** interpolate the Origen library over burnup to get transition
   ** structure data at the appropriate burnups, which is the only way
   ** to get a relevant conversion from flux to power/burnup.
   ** 
   ** Followup: Probably not that hard as it may be built into one of
   ** The ORIGEN objects.  Currently neglected out of ambivalance.
   TEST_F(OrigenInterfaceTester,fluxTest){
   EXPECT_THROW(tester.set_fluxes(dblank),cyclus::StateError);
    tester.set_fluxes(fluxes);
  */
}

TEST_F(OrigenInterfaceTester,powerTest){
  EXPECT_THROW(tester.set_powers(dblank),cyclus::StateError);
  EXPECT_NO_THROW(tester.set_powers(powers));

  std::vector<double> out_power = tester.get_powers();

  EXPECT_TRUE(out_power.size()>0);
  EXPECT_TRUE(out_power.size()==powers.size());

  for(size_t i = 0; i < out_power.size(); ++i) EXPECT_FLOAT_EQ(out_power[i],powers[i]) << "Powers fetched from interface object do not match those put in.";

  const double scaling_factor = 0.75;

  tester.set_power_scaling_factor(scaling_factor);

  out_power.clear();

  out_power = tester.get_powers();

  for(size_t i = 0; i < out_power.size(); ++i) EXPECT_FLOAT_EQ(out_power[i],powers[i]*scaling_factor) << "Powers fetched after scaling do not match those put in with scaling factor.";
}

TEST_F(OrigenInterfaceTester,solveTest){

  EXPECT_THROW(tester.solve(),cyclus::StateError);
  id_tags.erase("Something");
  params.erase("Fuel Temperature");
  params["Moderator Density"] = 0.7332;
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);
  tester.interpolate();
  EXPECT_THROW(tester.set_materials(ids,concs),cyclus::ValueError);

  tester.set_powers(powers);
  EXPECT_THROW(tester.set_materials(ids,concs),cyclus::ValueError);
  tester.delete_powers();

  dblank.push_back(1.);
  tester.set_time_steps(times);
  tester.set_powers(dblank);
  EXPECT_THROW(tester.set_materials(ids,concs),cyclus::ValueError);
  tester.delete_powers();

  tester.set_powers(powers);
  tester.set_materials(ids,concs);
  ASSERT_NO_THROW(tester.solve());

  // Check burnups for each step
  EXPECT_FLOAT_EQ(tester.burnup_at(0),0.);
  EXPECT_FLOAT_EQ(tester.burnup_at(1),10000.);
  EXPECT_FLOAT_EQ(tester.burnup_at(2),30000.);
  EXPECT_FLOAT_EQ(tester.burnup_at(3),60000.);
  EXPECT_FLOAT_EQ(tester.burnup_last(),60000.);

  // TODO: Add tests for other end-of-cycle parameters
  
  tester.reset_material();
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);
  tester.set_powers(powers);
  tester.add_time_step(500);
  tester.add_power(5.e6);
  tester.set_lib_path("/");
  EXPECT_NO_THROW(tester.interpolate());
  tester.set_materials(ids,concs);
  ASSERT_NO_THROW(tester.solve());

  EXPECT_FLOAT_EQ(tester.burnup_at(4),70000.);
  EXPECT_FLOAT_EQ(tester.burnup_last(),70000.);

}

TEST_F(OrigenInterfaceTester,resultTest){

  const size_t ORIGEN_LIB_SIZE = 2237;

  id_tags.erase("Something");
  params.erase("Fuel Temperature");
  params["Moderator Density"] = 0.7332;
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);
  tester.interpolate();

  tester.set_powers(powers);
  tester.set_time_steps(times); 
  ASSERT_TRUE(powers.size() == times.size()-1);
  tester.set_materials(ids,concs);

  ASSERT_NO_THROW(tester.solve());

  std::vector<int> ids_out;
  tester.get_ids(ids_out);
  EXPECT_EQ(ORIGEN_LIB_SIZE,ids_out.size()) << "Resulting ID vector is not of the correct size.";


  std::vector<std::vector<double> > masses_out;
  tester.get_masses(masses_out);
   
   
  EXPECT_EQ(times.size(),masses_out.size()) << "get_masses() returned an unexpected number of concentration vectors!";
  for(size_t i = 0; i < times.size(); i++){
    EXPECT_EQ(ORIGEN_LIB_SIZE,masses_out[i].size()) << "Masses vector #" << i 
              << " for time " << times[i] << " is of the incorrect size.";

    std::vector<double> mass_out;
    tester.get_masses_at(i,mass_out);
    for(size_t j = 0; j < masses_out[i].size(); j++){
      EXPECT_EQ(mass_out[j],masses_out[i][j]) 
         << "Disagreement between return of get_masses() and get_masses_at() for " 
         << ids_out [j] << " at time " << times[i] << ".";
    }
  }

  std::vector<double> mass_out;
  tester.get_masses_final(mass_out);
  size_t numTimes = times.size();

  EXPECT_EQ(mass_out.size(),masses_out[numTimes-1].size()) << "Size mismatch between final mass vector size.";
  for(size_t i = 0; i < masses_out[numTimes-1].size(); ++i){
    EXPECT_EQ(mass_out[i],masses_out[numTimes-1][i]) << "Disagreement between return of get_masses() and get_masses_final().";
  }

  std::string tm_test = tester.get_tag_manager_string();
  std::stringstream tm_ref;
  tm_ref << "{\n   \"Assembly Type\" : \"ce14x14\",\n" << \
               "   \"Enrichment\" : 0.040,\n" << \
               "   \"Fuel Type\" : \"Uranium\",\n" << \
               "   \"Moderator Density\" : 0.73320,\n" << \
               "   \"Powers (W)\" : \"1e+07,2e+07,3e+07\",\n" << \
               "   \"Times (d)\" : \"0,100,200,300\"\n}\n";
  EXPECT_EQ(tm_ref.str(),tm_test) << "Didn't get the tag manager string comparison right...";
}
