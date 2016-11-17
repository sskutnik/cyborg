#ifndef CYCLUS_CYBORG_REACTOR_H_
#define CYCLUS_CYBORG_REACTOR_H_

#include <string>
#include "cyclus.h"
//#include "cyclus_origen_interface.h"

/// @class Reactor
///
/// cyborg::Reactor is an ORIGEN-based reactor archetype capable of dynamically
/// calculating spent fuel recipes based on reactor physics and depletion 
/// calculations using ORIGEN
///
/// @section intro Introduction
///
/// The CyBORG Reactor is largely derived from the cycamore::Reactor class,
/// save for the fact that it uses ORIGEN for depletion to calculate discharge
/// fuel recipes as-needed. As a result, CyBORG infers reactor data interpolation
/// parameters from the input fuel recipe and fuel type (UOX/MOX/other), as well
/// as user-provided interpolation tags.
///
/// @section agentparams Agent Parameters
/// Place a description of the required input parameters which define the
/// agent implementation.
///
/// @section optionalparams Optional Parameters
/// Place a description of the optional input parameters to define the
/// agent implementation.
///
/// @section detailed Detailed Behavior
/// Place a description of the detailed behavior of the agent. Consider
/// describing the behavior at the tick and tock as well as the behavior
/// upon sending and receiving materials and messages.

//  Forward declaration for Cyclus-ORIGEN interface layer
namespace OrigenInterface {
   class cyclus2origen;
}

namespace cyborg {

//  Forward declaration used for test harness
namespace ReactorTests {
   class ReactorTest;
}

class Reactor : public cyclus::Facility,
    public cyclus::toolkit::CommodityProducer {
 public:
  /// Constructor for Reactor Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit Reactor(cyclus::Context* ctx);

 
  #pragma cyclus note {"doc": "CyBORG is an ORIGEN-based Reactor archetype.",\
                       "niche": "Reactor"}

  /// A verbose printer for the Reactor
  virtual std::string str();

  /// Sets up the Reactor's trade requests
  virtual void EnterNotify();

  /// The handleTick function specific to the Reactor.
  /// @param time the time of the tick
  virtual void Tick();

  /// The handleTick function specific to the Reactor.
  /// @param time the time of the tock
  virtual void Tock();

  virtual void Load_();

  /// Transmute the batch that is about to be discharged from the core to its
  /// fully burnt state 
  void Transmute_();
  
  /// Transmute the specified number of assemblies in the core to their
  /// fully burnt state 
  void Transmute_(int n_assem, int n_cycles=-1, double last_cycle=1.0);
 

  /// Records a Reactor event to the output db with the given name and note val.
  void Record(std::string name, std::string val);

  /// Discharge a batch from the core 
  bool Discharge_();

  /// Discharge the specified number of assemblies from the core 
  bool Discharge_(const int);

  ///  Deplete the material for this cycle
  ///  @param mat   The material to be depleted
  ///  @param n_assem Number of assemblies to be depleted (using assem_size => total HM mass) 
  ///  @param n_cycles Number of cycles to perform depletion (<= # of batches; used for final cycle)
  ///  @param cycle_length Fractional length (0,1] of last cycle (used for final / decom cycle)
  cyclus::Composition::Ptr Deplete_(cyclus::Material::Ptr, const int, const int, const double);

  double fuel_capacity() { return ( this->fresh.space() + this->core.space() ); }

  /// Store fuel info index for the given resource received on incommod.
  void index_res(cyclus::Resource::Ptr m, std::string incommod);

  /// Material buy policy procedures

  virtual void AcceptMatlTrades(const std::vector<std::pair<
      cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses);

  virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> GetMatlRequests();

  bool decom;

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).

  #pragma cyclus decl

 ////////// Data accessors /////////

 private:
  bool retired() {
     return exit_time() != -1 && context()->time() >= exit_time();
  }

  void setup_origen_interp_params(OrigenInterface::cyclus2origen&, const cyclus::Material::Ptr);
  void setup_origen_power_history(OrigenInterface::cyclus2origen&, const int, const double);
  void setup_origen_materials(OrigenInterface::cyclus2origen&, const cyclus::Material::Ptr, const int);
  cyclus::CompMap get_origen_discharge_recipe(OrigenInterface::cyclus2origen&);

  /* Module Members */

  /// Level 1 Parameters

  #pragma cyclus var {'default':[],\
                      "tooltip":"Input fuel commodity",\
                      "doc":"Fuel name accepted by this Reactor",\
                      "uitype":["oneormore","incommodity"],\
                      "uilabel":"Fuel Commodity",\
                      "userlevel":1}
  std::vector<std::string> fuel_incommods;

  #pragma cyclus var {'default':[],\
                      "tooltip":"Input fuel recipe",\
                      "doc":"Fuel recipe accepted by Reactor",\
                      "uitype":["oneormore","recipe"],\
                      "uilabel":"Fuel Recipe",\
                      "userlevel":1}
  std::vector<std::string> fuel_recipes;                      

  #pragma cyclus var { \
    "default": [], \
    "uitype":"oneormore",\
    "uilabel": "Fresh Fuel Preference List", \
    "doc": "The preference for each type of fresh fuel requested corresponding"\
           " to each input commodity (same order).  If no preferences are " \
           "specified, 1.0 is used for all fuel requests (default).", \
  }
  std::vector<double> fuel_prefs;

  #pragma cyclus var {'default':"spent_fuel",\
                      "tooltip":"Spent fuel commodity",\
                      "doc":"Name of spent fuel commodity",\
                      "uitype":"outcommodity","uilabel":"Spent Fuel Commodity",\
                      "userlevel":1}
  std::string spent_fuel;

  #pragma cyclus var {'default':"power",\
                      "tooltip":"Reactor power name",\
                      "doc":"Name of commodity Reactor produces",\
                      "uilabel":"Power Name"}
  std::string power_name;

  #pragma cyclus var {"tooltip":"Thermal Power capacity",\
                      "doc":"Reactor thermal power capacity (MW)",\
                      "units":"MW",\
                      "uilabel":"Power Capacity",\
                      "userlevel":1}
  double power_cap;

 //////////// inventory and core params ////////////
  #pragma cyclus var { \
    "doc": "Mass (kg) of a single assembly", \
    "uilabel": "Assembly Mass", \
    "units": "kg", \
    "userlevel":1 }
  double assem_size;

  #pragma cyclus var { \
    "uilabel": "Number of Assemblies in Core", \
    "doc": "Number of assemblies that constitute a full core.", \
    "userlevel":1 }
  int n_assem_core;

  #pragma cyclus var { \
    'default': "n_assem_core",\
    "uilabel": "Number of Assemblies per Batch", \
    "doc": "Number of assemblies that constitute a single batch.\n" \
           "This is the number of assemblies discharged from the core fully " \
           "burned each cycle.\n"\
           "Batch size is equivalent to ``n_assem_batch / n_assem_core``.\n"\
           "Defaults to n_assem_core (single-batch core).", \
    "userlevel":1 }
  int n_assem_batch;

  #pragma cyclus var { \
    "default": 0, \
    "uilabel": "Minimum Fresh Fuel Inventory", \
    "units": "assemblies", \
    "doc": "Number of fresh fuel assemblies to keep on-hand if possible.", \
    "userlevel":1 }
  int n_assem_fresh;

  #pragma cyclus var { \
    "default": 1000000000, \
    "uilabel": "Maximum Spent Fuel Inventory", \
    "units": "assemblies", \
    "doc": "Number of spent fuel assemblies that can be stored on-site before" \
           " Reactor operation stalls.", \
    "userlevel":1 }
  int n_assem_spent;


  ///////// cycle params ///////////
  #pragma cyclus var { \
     "doc": "The duration of a full operational cycle (excluding refueling " \
            "time) in time steps.", \
     "uilabel": "Cycle Length", \
     "units": "time steps" }
  int cycle_time;

  #pragma cyclus var { \
     "doc": "The duration of a full refueling period - the minimum time between"\
            " the end of a cycle and the start of the next cycle.", \
     "uilabel": "Refueling Outage Duration", \
     "units": "time steps" }
  int refuel_time;

  #pragma cyclus var { \
     "default": 0, \
     "doc": "Number of time steps since the start of the last cycle." \
            " Only set this if you know what you are doing", \
     "uilabel": "Time Since Start of Last Cycle", \
     "units": "time steps"}
  int cycle_step;
 
  /// fuel_type is used to determine which fuel recipe parameters 
  /// to extract for Reactor data library interpolation.
  /// e.g., UOX uses U-235 enrichment
  ///       MOX uses Pu-239 enrichment and total Pu fraction
  ///       "Other" does not use recipe-based values for interpolation
  /// \TODO: Add option for one-or-more selection (i.e., one per input recipe...)
  #pragma cyclus var {'default': "UOX",\
                      "tooltip":"Reactor fuel type",\
                      "uilabel":"Fuel type",\
                      "doc":"Fuel type (UOX/MOX/other) used for the reactor."\
                            " Used for reactor data library interpolation from"\
                            " initial fuel composition.",\
                      'uitype':"combobox",\
                      'categorical':['UOX','MOX','other'],\
                      "userlevel":1}
  std::string fuel_type;

  /// Level 2 Parameters  
  /// Default set by environment / build argument in orglib_default_location.h
  /// Extract the ORIGEN library default location from auto-generated header, 
  /// because cycpp needs a string, not a preprocessor directive
  #pragma cyclus exec import os; import re; \
                      orglib_prefix = "./" if os.path.isfile('cyborg/orglib_default_location.h') else "../";\
                      orglib_header = orglib_prefix + 'cyborg/orglib_default_location.h'; \
                      orglib_loc = re.search('ORIGEN_LIBS_DEFAULT\s+\"(.*)\"', open(orglib_header,'r').read()).group(1)
  #pragma cyclus var {'default': orglib_loc,\
                      "tooltip":"Path to ORIGEN Libraries",\
                      "doc":"Path to ORIGEN Libraries",\
                      "uilabel":"Path to ORIGEN Libraries",\
                      "userlevel":2}
  std::string lib_path; 

  #pragma cyclus var {'default':"w17x17",\
                      'tooltip':"Assembly type",\
                      'doc':"ORIGEN library type to be used",\
                      'uilabel':"Assembly Type",\
                      'uitype':"combobox",\
                      'userlevel':2,\
                      'categorical':['ce14x14', 'ce16x16', 'w14x14', 's14x14', 'w15x15', \
                      'w17x17', 'w17x17_ofa', 'ge7x7-0', 'ge8x8-4', 'abb8x8-1', 'ge9x9-7', 'ge10x10-8', \
                      'atrium9-9', 'atrium10-9', 'svea64-1', 'svea100-0']} 
  std::string assembly_type;
      
  /// Level 3 Parameters

  /// Interpolation tags to be used directly with ORIGEN TagManager
  /// Reactor data library interpolation parameters expressed as key-value pairs
  /// Example tags include moderator density, fuel temperature, etc.  
  
  #pragma cyclus var {'default':{},\
                      'tooltip':"ORIGEN library interpolation tag/value pairs",\
                      'uilabel':"Interp. tags",\
                      'alias': ["tags", "tag", "value"],\
                      'uitype': ["oneormore", "string","double"],\
                      'userlevel':3,\
                      'doc':"Tag/value pairs used for ORIGEN reactor data library"\
                            " interpolation. Valid tags depend on the individual"\
                            " library, but may include aspects such as moderator"\
                            " density, fuel temperature, etc.\n\n"\
                            "Note that fuel composition tags should NOT be input"\
                            " here, as they are covered by the fuel_type input and"\
                            " input composition recipe."}
  std::map<std::string,double> interp_tags;

  #pragma cyclus var {'default':[],\
                      'uilabel':'Batch core power frac.',\
                      'tooltip':"Core power fraction for each fuel batch (0,1])",\
                      'uitype': ["oneormore","double"],\
                      'doc': "Core power fraction for each fuel batch; expressed"\
                             " as a value from (0,1], totaling to 1. Unnormalized"\
                             " values are automatically renormalized. Defaults to"\
                             " equal power fraction per cycle if unspecified. MUST"\
                             " be the same number of entries as batches in the core,"\
                             " as determined by n_assem_core / n_assem_batch.",\
                      'userlevel':3\
                     }
  std::vector<double> core_power_frac;

  // should be hidden in ui (internal only). 
  // True if fuel has already been discharged this cycle.
  #pragma cyclus var {"default": 1, "doc": "This should NEVER be set manually",\
                      "internal": True \
  }
  bool discharged;

  // This variable should be hidden/unavailable in ui.  Maps resource object
  // id's to the index for the incommod through which they were received.
  #pragma cyclus var {"default": {}, "doc": "This should NEVER be set manually", \
                      "internal": True }
  std::map<int, int> res_indexes;

  // should be hidden in ui (internal only). 
  // State variable to store the current discharge burnup of the fuel
  #pragma cyclus var {"default": 0.0, "doc": "This should NEVER be set manually",\
                      "internal": True, "units":"MWd/MTHM" }
  double burnup;
  
  /// Material Flow Parameters                                                      
                                                                     
  // referenced (e.g. n_batch_fresh, assem_size, etc.).
  #pragma cyclus var {"capacity": "n_assem_fresh * assem_size",\
                      "tooltip": "Incoming fresh fuel buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> fresh;
  #pragma cyclus var {"capacity":"n_assem_core * assem_size",\
                      "tooltip":"Fuel held in core"}
  cyclus::toolkit::ResBuf<cyclus::Material> core;
  #pragma cyclus var {"capacity": "n_assem_spent * assem_size",\
                      "tooltip":"Discharge fuel holding capacity / output buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> spent;
  
  // Cached composition of spent fuel (calculated from ORIGEN)
  cyclus::Composition::Ptr spentFuelComp;
                                                                   
  //// A policy for requesting material
  //cyclus::toolkit::MatlBuyPolicy buy_policy;
  //// A policy for sending material
  cyclus::toolkit::MatlSellPolicy sell_policy;

  /// Allow test access to private members
  friend class cyborg::ReactorTests::ReactorTest;
  
};

// Helper functions

/// \brief Get heavy metal mass fraction of an element 
/// \param Z      Atomic number of element (e.g., U = 92)
/// \param comp   Cyclus composition from which to extract HM mass fraction
/// \param hm_ele Starting atomic number for summing heavy metal masses (default = 90)
double get_ele_hm_mass_frac(const int, const cyclus::Composition::Ptr, const int Z_HM = 90);

/// \brief Get isotopic mass fraction for a given isotope (e.g., % U-235 in U)
/// \param Z    Atomic number of parent element
/// \param A    Mass number for isotope mass fraction
/// \param comp Cyclus composition to extract isotopic mass fraction 
double get_iso_mass_frac(const int, const int, const cyclus::Composition::Ptr);

}  // namespace cyborg

#endif  // CYCLUS_CYBORG_REACTOR_H_
