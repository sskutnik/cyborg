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

     powers.push_back(1.e7);
     powers.push_back(2.e7);
     powers.push_back(3.e7);

     times.push_back(0.);
     times.push_back(100.);
     times.push_back(200.);
     times.push_back(300.);  

     // Initial concentrations; default units of kg
     concs.push_back(4.0);
     concs.push_back(96.0);

     ids.push_back(922350);
     ids.push_back(922380);

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
  // Tests for failure, not correctness.
  // No methods currently in place to test for correctness.

  EXPECT_THROW(tester.interpolate(),cyclus::StateError);

  id_tags.erase("Something");
  tester.set_id_tags(id_tags);

  params.erase("Moderator Density");
  params.erase("Fuel Temperature");
  tester.set_parameters(params);

  EXPECT_THROW(tester.interpolate(),cyclus::ValueError);

  tester.set_lib_path(lib_path);
  std::cout << "Expect next line to be warning about unspecified tag 'Moderator Density'." << std::endl;
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

  EXPECT_NO_THROW(tester.set_materials(ids,concs));
}
/* Commented out until such time as fluxes are properly implemented.
** Currently hindered by the need to translate to flux in order to
** interpolate the Origen library over burnup to get transition
** structure data at the appropriate burnups, which is the only way
** to get a relevant conversion from flux to power/burnup.
TEST_F(OrigenInterfaceTester,fluxTest){
  EXPECT_THROW(tester.set_fluxes(dblank),cyclus::StateError);
  tester.set_fluxes(fluxes);
}
*/
TEST_F(OrigenInterfaceTester,powerTest){
  EXPECT_THROW(tester.set_powers(dblank),cyclus::StateError);
  EXPECT_NO_THROW(tester.set_powers(powers));
  // Test that powers were set correctly
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
  tester.set_materials(ids,concs);
  EXPECT_THROW(tester.solve(),cyclus::StateError);

  tester.set_powers(powers);
  EXPECT_THROW(tester.solve(),cyclus::ValueError);
  tester.delete_powers();

  dblank.push_back(1.);
  tester.set_time_steps(times);
  tester.set_powers(dblank);
  EXPECT_THROW(tester.solve(),cyclus::ValueError);
  tester.delete_powers();

  tester.set_powers(powers);
  ASSERT_NO_THROW(tester.solve());

  // Check burnups for each step
  EXPECT_FLOAT_EQ(tester.burnup_at(0),10000.);
  EXPECT_FLOAT_EQ(tester.burnup_at(1),30000.);
  EXPECT_FLOAT_EQ(tester.burnup_at(2),60000.);
  EXPECT_FLOAT_EQ(tester.burnup_last(),60000.);

  // TODO: Add tests for other end-of-cycle parameters
  
  tester.reset_material();
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);
  tester.set_powers(powers);
  tester.add_time_step(500);
  tester.add_power(5.e6);
  tester.interpolate();
  tester.set_materials(ids,concs);
  ASSERT_NO_THROW(tester.solve());

  EXPECT_FLOAT_EQ(tester.burnup_at(3),70000.);
  EXPECT_FLOAT_EQ(tester.burnup_last(),70000.);

}

TEST_F(OrigenInterfaceTester,resultTest){
  id_tags.erase("Something");
  params.erase("Fuel Temperature");
  params["Moderator Density"] = 0.7332;
  tester.set_id_tags(id_tags);
  tester.set_parameters(params);
  tester.set_lib_path(lib_path);
  tester.interpolate();

  tester.set_materials(ids,concs);
  tester.set_powers(powers);
  tester.set_time_steps(times); 
  ASSERT_TRUE(powers.size() == times.size()-1);

  ASSERT_NO_THROW(tester.solve());

  std::vector<int> ids_out;
  tester.get_ids(ids_out);
  EXPECT_EQ(1946,ids_out.size()) << "Resulting ID vector is not of the correct size.";


  std::vector<std::vector<double> > masses_out;
  tester.get_masses(masses_out);
   
   
  EXPECT_EQ(times.size(),masses_out.size()) << "get_masses() returned an unexpected number of concentration vectors!";
  for(size_t i = 0; i < times.size(); i++){
    EXPECT_EQ(1946,masses_out[i].size()) << "Masses vector #" << i << " for time " << times[i] << " is of the incorrect size.";

    std::vector<double> mass_out;
    tester.get_masses_at(i,mass_out);
    for(size_t j = 0; j < masses_out[i].size(); j++){
      EXPECT_EQ(mass_out[j],masses_out[i][j]) << "Disagreement between return of get_masses() and get_masses_at() for " << ids_out [j] << " at time " << times[i] << ".";
    }
  }

  std::vector<double> mass_out;
  tester.get_masses_final(mass_out);
  size_t numTimes = times.size();

  EXPECT_EQ(mass_out.size(),masses_out[numTimes-1].size()) << "Size mismatch between final mass vector size.";
  for(size_t i = 0; i < masses_out[numTimes-1].size(); ++i){
    EXPECT_EQ(mass_out[i],masses_out[numTimes-1][i]) << "Disagreement between return of get_masses() and get_masses_final().";
  }
}
/*

int main(int argc, char ** argv){
  OrigenInterface::cyclus2origen tester;
// Library names should be specified as string literals.  Names can
// be given as absolute or relative pathnames if the library is not
// in the current directory.

  tester.set_lib_names({"ce14_e20.arplib","ce14_e30.arplib","ce14_e40.arplib"});

  tester.list_lib_names();
  std::cout << std::endl;

  tester.remove_lib_names({"ce14_e40.arplib","ce14_e30.arplib"});

  tester.list_lib_names();
  std::cout << std::endl;

  tester.add_lib_names({"ce14_e30.arplib"});

  tester.list_lib_names();
  std::cout << std::endl;

  tester.remove_lib_names({"ce14_e30.arplib","ce14_e20.arplib"});

  tester.set_lib_path("/home/nsly/scale_dev_data/arplibs");
  tester.list_lib_names();
  std::cout << std::endl;

// Setting the ID tags are done one-by-one giving the name and
// value of the ID tags as string literals.  These will be used
// to collect the appropriate libraries from the ones provided,
// allowing the list of libraries to actually include more than
// those actually used for this calculation.  This does require
// that the libraries have the tags set.  This can be checked
// using the obiwan tool included with Origen.
//
// usage -- obiwan tag [library name]

  tester.set_id_tag("Assembly Type","ce14x14");
  tester.set_id_tag("Fuel Type","Uranium");
  tester.set_id_tag("Not real","Fake");

  tester.list_id_tags();
  std::cout << std::endl;

  tester.remove_id_tag("Not real");
  tester.remove_id_tag("Fuel Type");

  tester.list_id_tags();
  std::cout << std::endl;

  tester.set_id_tag("Fuel Type","Uranium");

  tester.list_id_tags();
  std::cout << std::endl;

  tester.remove_id_tag("Assembly Type");
  tester.remove_id_tag("Fuel Type");

  std::map<std::string,std::string> id_tags;

  id_tags["Assembly Type"] = "ce14x14";
  id_tags["Fuel Type"] = "Uranium";

  tester.set_id_tags(id_tags);

  tester.list_id_tags();
  std::cout << std::endl;

// Adding the parameters adds the interpolable tags to the Tag
// Manager on the cyclus2origen object (here: tester).  These
// will be used by the interpolate function to interpolate
// between the libraries set in set_lib_names that match the
// ID tags set in set_id_tags to generate a single library that
// will be used in the set_materials function.  This allows for
// the use of problem-specific libraries that may not precisely
// match those already on-disk.

  tester.add_parameter("Enrichment",2.5);
  tester.add_parameter("Fuel Temp",9000);

  tester.list_parameters();
  std::cout << std::endl;

  tester.remove_parameter("Fuel Temp");

  tester.list_parameters();
  std::cout << std::endl;

  tester.remove_parameter("Enrichment");

  std::map<std::string,double> interp_tags;

  interp_tags["Enrichment"] = 2.5;
  interp_tags["Fuel Temp"] = 9000.;

  tester.set_parameters(interp_tags);

  tester.list_parameters();
  std::cout << std::endl;

  tester.remove_parameter("Fuel Temp");

  tester.list_parameters();
  std::cout << std::endl;

// Interpolation using the aforementioned libraries and tags.
// This will generate a single library at the user-defined
// point in the parameter space.  This library's cross sections
// will be used as the basis for the depletion calculation.

  std::cout << "Expect warning about undefined interp value:" << std::endl;
  tester.interpolate();

// Setting the volume of the material in the problem.

//  tester.set_volume(3);

// Setting the vector of nuclide IDs using the sizzzaaa format.
// Will allow for alternative formats in the future.

  std::vector<int> ids;
  ids.push_back(20092235);
  ids.push_back(20092238);

// Setting the number densities of the isotopes specified in
// the ids vector.  The position of the densities in this
// vector should directly correspond with the isotope vector.

  std::vector<double> numden;
  numden.push_back(5);
  numden.push_back(95);

// Uses the library generated by interpolate and the volume
// along with the IDs and number densities provided in the
// input parameters to create the Material object used to
// simulate the depletion calculation.

  tester.set_materials(ids,numden);

// Specifying the fluxes to be used in the depletion
// calculation.  There should be a flux for each burn step.

//  std::vector<double> fluxes(3,1.0e14);
//  tester.set_flux(fluxes);

  std::vector<double> powers(4,200);
  tester.set_powers(powers);

// Specifying the times for which isotopics will be kept.
// Should have one more element than the fluxes, as the
// burn steps will occur between the times specified here.
// Units are specified in their own call using strings
// to spell out the units (i.e. "seconds", "days","years").
// This unit will apply to the entire times vector.

  std::vector<double> times(4,0.0);
  times[1] = 200;
  times[2] = 500;
  times[3] = 5000;
  tester.set_time_steps(times);
  tester.set_time_units("seconds");

// Calling the solver to calculate the isotopics after
// every depletion step.

  tester.solve();

// Getting the concentration vectors and printing the number
// of vectors received.  Should be equal to the size of the times vector.

  std::vector<std::vector<double>> concs_vecs;
  tester.get_concentrations(concs_vecs);
  std::cout << "Total concentration vectors calculated: " << concs_vecs.size() << " == " << times.size() << " times." << std::endl;

// Getting the concentrations after the first burn step
// and printing the size of the vector (number of concentrations.)

  std::vector<double> concs;
  tester.get_concentrations_at(1,concs);
  std::cout << "Concentration set 2 has size: " << concs.size() << "." << std::endl;

// Getting the concentrations after the first burn step
// and printing the size of the vector (number of concentrations.)

  std::vector<double> concs_final;
  tester.get_concentrations_final(concs_final);
  std::cout << "Last concentration set has size: " << concs_final.size() << "." << std::endl;

// Getting the nuclide ids (sizzzaaa format) for the entire
// problem and printing the size of the vector.  Should equal
// the number of concentrations at any step.

  std::vector<int> id_out;
  tester.get_ids(id_out);
  std::cout << "ID vector has size: " << id_out.size() << "." << std::endl;

  return 0;
}
*/
