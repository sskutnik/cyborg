#include <string>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include "Origen/Core/Definitions.h"
#include "Origen/Core/dc/ConcentrationUnit.h"
#include "Origen/Core/dc/Library.h"
#include "Origen/Core/dc/Material.h"
#include "Origen/Core/dc/Power.h"
#include "Origen/Core/dc/Time.h"
//#include "Origen/Core/xf/Solver.h"
//#include "Origen/Solver/cram/Solver_cram.h"
#include "Origen/Solver/SolverSelector.h"

/*!
 * \note - 
 *   The internal units for Origen::Material are:
 *   - amount of nuclides \f$ {cm^2}/{barn} \f$
 *   - time in \f$ s \f$ 
 *   - flux in \f$ {\#}/{cm^2 s} \f$
 *   - power in \f$ W \f$
 */

/*!
 * \note - Envisioned usage is:
 *         set_lib_name->set_id_tag(somewhat optional)->
 *         add_parameter(optional)->interpolate->
 *         set_volume->set_materials->set_flux->
 *         set_time_steps->
 *         set_**_units(as necessary)->
 *         set_solver(if made available)->solve.
 *         
 *         After that, concentrations at various
 *         time steps and corresponding ids should
 *         be available through the get_* calls.
 *
 *         Streamlined version would be:
 *         set_lib_path->set_id_tag->add_parameter->interpolate
 *           --Gets you a problem specific library at default burnups.
 *
 *         set_materials_with_masses->solve(times,fluxes)->
 *         get_ids->get_masses.
 *           --Gets you material concentrations at specified burnups.
 *
 *         set_materials->solve(times,fluxes)->get_ids->get_concentrations
 *           --After initial solve, materials can be set based on output
 *             of previous solve to solve for second cycle.
 */

namespace OrigenInterface {

class cyclus2origen {

public:

    /*!
     ** \brief - Function to set the library names including paths
     **          that are either relative or absolute paths.
     ** \param - String representing a reactor assembly type for
     **          which libraries exist on disk.
     */
    void set_lib_names(const std::vector<std::string>&);

    /*!
     ** \brief - Function to set the path to a directory of libraries.
     ** \param - String prepresenting a path to tagged ORIGEN
     **          libraries for use in interpolation and depletion.
     */
    void set_lib_path(const std::string);

    /*!
     ** \brief - Function to add library names to the already
     **          existing list.
     ** \param - Vector of library names as strings.  Can be in
     **          absolute or relative path formats.
     */
    void add_lib_names(const std::vector<std::string>&);

    /*!
     ** \brief - Function to remove library names from internally
     **          stored list.
     ** \param - Vector of library names as strings to be removed.
     **          Must be exactly the same.
     */
    void remove_lib_names(const std::vector<std::string>&);

    /*!
     ** \brief - Function to list library names currently stored
     **          to std::cout.
     */
    void list_lib_names() const;

    /*!
     ** \brief - Function to get a copy of library names in memory.
     ** \param - Empty vector of strings to be populated with library names.
     */
    void get_lib_names(std::vector<std::string>&) const;

    /*!
     ** \brief - Function to set the id_tags used by the TagManager
     **          during the set_materials function to determine which
     **          libraries to use during interpolation and depletion.
     ** \param - String representing an id_tag name on TagManager.
     ** \param - String representing an id_tag value on TagManager.
     */
    void set_id_tag(const std::string, const std::string);

    /*!
     ** \brief - Function to set multiple id_tags at a time using
     **          a standard map container holding all of the
     **          name-value pairs for the tags.
     ** \param - Map containing the name-value pairs for the tags.
     */
    void set_id_tags(const std::map<std::string,std::string>&);

    /*!
     ** \brief - Delete an ID tag.
     ** \param - String name of ID tag to be removed.
     */
    void remove_id_tag(const std::string);

    /*!
     ** \brief - List ID tags currently stored and their values to cout.
     */
    void list_id_tags() const;

    /*!
     ** \brief - Get a copy of ID tags currently stored.
     ** \param - Empty vector of strings to be populated with ID tag names.
     ** \param - Empty vector of strings to be populated with ID tag values.
     */
    void get_id_tags(std::vector<std::string>&, std::vector<std::string>&) const;

    /*!
     ** \brief - Function to set the initial materials to be depleted
     **          that accepts masses as the quantity input.  This function
     **          converts the masses into concentration units (cm_2_barn)
     **          and then uses that to call the set_materials function.
     ** \param - Vector of ids.  As this function calls set_materials,
     **          the types of IDs accepted are the same.
     ** \param - Vector of masses corresponding to the vector of IDs.
     **          Units are currently hard-coded to kilograms.
     */
    void set_materials_with_masses(std::vector<int>&,const std::vector<double>&);

    /*!
    ** \brief - Function to set the initial materials to be depleted.
    ** \param - Vector of ids.  Works with ids formatted in either
    **          zzzaaai or pizzzaaa.
    ** \param - Vector of concentrations.  Should directly correspond
    **              to the vector of ids.
    */
    void set_materials(std::vector<int>&, const std::vector<double>&);

    /*!
     ** \brief - Function to set the units used in the concentrations vectors.
     ** \param - String declaring the units used for the materials concentrations.
     */
    void set_mat_units(const std::string);

    /*!
     ** \brief - Function to set the volume of material being depleted.
     ** \param - Volume of the material being processed.  in g/cc?
     */
//    void set_volume(const double vol){b_vol=vol;};

    /*!
     ** \brief - Function to set the times at which a burn step will end,
     **          and excepting the final value, begin.
     ** \param - Vector of values for the time the materials will be
     **          exposed to the specified fluxes/burned/decayed.
     */
    void set_time_steps(const std::vector<double>& times){b_times=times;};

    /*!
     ** \brief - Function to specify the time unit used in the vector
     **          in set_time_steps.  
     ** \param - String declaring the units for the time step values.
     **          Accepts full words in all lower case (i.e. - 'seconds',
     **          'days','years').
     */
    void set_time_units(const char* time_units);

    /*!
     ** \brief - Function to add a time step on to the end of an existing
     **          set of times and fluxes.
     ** \param - New time to target for calculation.  Must be greater
     **          than the greatest current time.  If we want to add a
     **          new time step in the middle, the calculation will have
     **          to be redone.
     ** \param - New flux to correspond with the new time.
     */
//    void add_time_step(double, double, double);

    /*!
     ** \brief - Function to input the fluxes for each burn step.
     ** \param - Vector of fluxes for each time step.
     **          Units of n/cm2*s.
     */
    void set_flux(const std::vector<double>&);

    /*!
     ** \brief - Function to input the powers for each burn step.
     ** \param - Vector of powers for each time step.
     **          Units of watts.
     */
    void set_powers(const std::vector<double>&);

    /*!
     ** \brief - Setting a new tag or changing the value
     **          of an already existing tag.
     ** \param - String for interp_tag name on TagManager.
     ** \param - Double for interp_tag value on TagManager.
     */
    void add_parameter(const std::string, const double);

    /*!
     ** \brief - Method for adding parameters all at once
     **          by using a std::map.
     ** \param - std::map of parameter/interp_tag names and
     **          values, std::string and double, respectively.
     */
    void add_parameters(const std::map<std::string,double>&);

    /*!
     ** \brief - Function to set multiple parameters at once.
     ** \param - Map of name-value pairs for the parameter tags.
     */
    void set_parameters(const std::map<std::string,double>&);

    /*!
     ** \brief - Remove an interpolable tag.
     ** \param - String name of tag to be removed.
     */
    void remove_parameter(const std::string);

    /*!
     ** \brief - Print out list of paramters (interp_tags).
     */
    void list_parameters() const;

    /*!
     ** \brief - Get copies of parameter names and values.
     ** \param - Empty vector of strings to be populated with
     **          currently stored parameter names.
     ** \param - Emtpy vector of doubles to be populated with
     **          currently stored parameter values.
     */
    void get_parameters(std::vector<std::string>&,std::vector<double>&) const;

// choose solver?  allow for choice at all?
    // 0: cram
    // 1: ??
    // 2: ??
//    void set_solver(const int svlr=0);

    /*!
     ** \brief - Function to interpolate between the specified
     **          libraries in set_lib_names that have the IDs
     **          specified in set_id_tag to the parameter values
     **          specified in add_parameter to generate a new
     **          library to be used in the depletion calculation.
     */
    void interpolate();

    /*!
     ** \brief - Function to call the solver.  Assumes time and
     **          flux/power vectors have been populated by their
     **          respective functions.
     */
    void solve();
    
    /*!
     ** \brief - Function to call the solver given a vector of
     **          times and fluxes.
     ** \param - Vector of times at which concentrations are
     **          to be calculated.
     ** \param - Vector of fluxes for the burn steps between
     **          each of the time points in the previous parameter.
     **          As such, this vector should have one less element
     **          than the times vector.  If empty, the powers
     **          parameter will be used.
     ** \param - Vector of powers for the burn steps between
     **          each of the time points in the times parameter.
     **          As such, this vector should have one less element
     **          than the times vector.  If empty, the fluxes
     **          parameter will be used.
     ** \note - Only the fluxes or powers parameters should be
     **         populated vectors.  The other should have size==0.
     */
    
    void solve(std::vector<double>& times,std::vector<double>& fluxes,std::vector<double>& powers);

    /*!
     ** \brief - Function to retrieve all of the concentrations at
     **          all of the times from the calculation.
     ** \param - Vector of vectors to be filled with all final 
     **          concentration values.
     */
    void get_concentrations(std::vector<std::vector<double>>&) const;

    /*!
     ** \brief - Function to retrieve the concentrations at a
     **          particular time from the times vector.
     ** \param - Integer indicating the burn step for which 
     **          concentrations are desired.
     ** \param - Vector to be filled with the final concentration
     **          values for a particular burn step.
     */
    void get_concentrations_at(int, std::vector<double>&) const;

    /*!
     ** \brief - Function to retrieve the concentrations at the end
     **          of the final burn step.
     ** \param - Vector to be filled with resulting concentrations
     **          of final burn step.
     */
    void get_concentrations_final(std::vector<double>&) const;

    void get_masses(std::vector<std::vector<double>>&) const;

    void get_masses_at(int,std::vector<double>&) const;

    void get_masses_final(std::vector<double>&) const;

    /*!
     ** \brief - Function to return the IDs that correspond to the
     **          concentrations fetched by the previous functions.
     ** \param - Vector of vectors to be filled with final ids,
     **          which are applicable to every burn step. Given in
     **          pizzzaaa format.
     */
    void get_ids(std::vector<int>&) const;

    /*!
     ** \brief - Function to return the ID vector using the zzzaaai
     **          format available in SCALE and in Cyclus.
     ** \param - Vector of ints into which the IDs will be placed.
     */

    void get_ids_zzzaaai(std::vector<int>&) const;

protected:

    /*!
     ** \brief - Function to take the library resulting from interpolate()
     **          and interpolate over the burnup dimension to the burnups
     **          calculated using the times and fluxes/powers.
     ** \param - Library which will be interpolated.
     ** \param - Vector of times for new burnup calculation.  Assumed to
     **          be in units of seconds unless b_time_units is set to a
     **          different unit.
     ** \param - Vector of fluxes that can be used to determine the burnups.
     **          Should have one less element than times vector.
     ** \param - Vector of powers that can be used to determine the burnups.
     **          Units are in watts.
     ** \note - It is assumed that either fluxes or powers will be an empty
     **         vector.  Anything else will throw an error.
     */
    void prob_spec_lib(Origen::SP_Library,std::vector<double> &times,std::vector<double> &fluxes,std::vector<double> &powers);

    Origen::SP_Library b_lib;
    Origen::SP_TagManager b_tm;
    Origen::SP_Material b_mat;
    Origen::SP_NuclideSet b_nucset;
//    Solver_cram b_slv;
    std::vector<std::string> b_lib_names;
    std::string b_lib_path;
    std::string b_interp_name;
    const double b_vol=1.;
    std::vector<int> b_ids;
    std::vector<double> b_concs;
    std::vector<double> b_burnups;
    std::vector<double> b_fluxes;
    std::vector<double> b_powers;
    std::vector<double> b_times;
    Origen::ConcentrationUnit b_concs_units;
    Origen::Time::UNITS b_time_units; // Default is seconds. Also accepts minutes,hours,days, and years.
    Origen::Power::UNITS b_power_units; // Default is watts.
};
} // end namspace

