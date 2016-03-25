/*!
 *  \file - tst_cyclus_origen_interface.cpp
 *  \author - Nicholas C. Sly
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include <gtest/gtest.h>

#include "cyclus_origen_interface.h"
#include "error.h"

using namespace OrigenInterface;

class myError : public cyclus::Error{
  public:
    myError(std::string msg) : Error(msg) {}
};

class OrigenInterfaceTester : public ::testing::Test {
protected:

//  void SetUp()
//  {
//// define setup that requires initial values with 'new' calls here
//  }

// define variables that don't require new calls here.
  cyclus2origen tester;
};

TEST_F(OrigenInterfaceTester,libManipulation)
{
  std::vector<std::string> lib_names;
  lib_names.push_back("ce14_e20.arplib");
  lib_names.push_back("ce14_e30.arplib");
  lib_names.push_back("ce14_e40.arplib");
  tester.set_lib_names(lib_names);

  std::vector<std::string> get_names;
  tester.get_lib_names(get_names);
  if(get_names!=lib_names){
    myError except=myError("Libraries not emplaced onto interface object properly.");
    except.what();
  }
}

//TEST_F(OrigenInterfaceTester,
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
