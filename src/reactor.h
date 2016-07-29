#ifndef CYCLUS_CYBORG_REACTOR_H_
#define CYCLUS_CYBORG_REACTOR_H_

#include <string>

#include "cyclus.h"
#include "orglib_default_location.h"

namespace cyborg {

/// @class reactor
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The reactor class inherits from the Facility class and is
/// dynamically loaded by the Agent class when requested.
///
/// @section intro Introduction
/// Place an introduction to the agent here.
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
//class reactor : public cyclus::Facility  {
class reactor : public cyclus::Facility,
    public cyclus::toolkit::CommodityProducer {
 public:
  /// Constructor for reactor Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit reactor(cyclus::Context* ctx);

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)

  #pragma cyclus

  #pragma cyclus note {"doc": "A stub facility is provided as a skeleton for the design of new facility agents."}

  /// A verbose printer for the reactor
  virtual std::string str();

  /// Sets up the Reactor's trade requests
  virtual void EnterNotify();

  /// The handleTick function specific to the reactor.
  /// @param time the time of the tick
  virtual void Tick();

  /// The handleTick function specific to the reactor.
  /// @param time the time of the tock
  virtual void Tock();

  virtual void Load_();


  /// Transmute the batch that is about to be discharged from the core to its
  /// fully burnt state 
  void Transmute_();
  
  /// Transmute the specified number of assemblies in the core to their
  /// fully burnt state 
  void Transmute_(int n_assem);
 

  /// Records a reactor event to the output db with the given name and note val.
  void Record(std::string name, std::string val);

  /// Discharge a batch from the core 
  bool Discharge_();

  /// Discharge the specified number of assemblies from the core 
  bool Discharge_(const int);

  /// Deplete the material for this cycle
  /// @param mat   The material to be depleted 
  /// @param power Depletion power of the material in MW
  cyclus::Composition::Ptr Deplete_(cyclus::Material::Ptr, double);

  bool decom;


 ////////// Data accessors /////////
 int get_cycle_time() { return this->cycle_time; }
 int get_cycle_step() { return this->cycle_step; }
 int get_refuel_time() { return this->refuel_time; }
 int get_n_assem_fresh() { return this->n_assem_fresh; }
 int get_n_assem_spent() { return this->n_assem_spent; }
 int get_n_assem_batch() { return this->n_assem_batch; }
 int get_n_assem_core() { return this->n_assem_core; }

 private:
   bool retired() {
      return exit_time() != -1 && context()->time() >= exit_time();
  }

  /* Module Members */

  /// Level 1 Parameters

  #pragma cyclus var {"tooltip":"Input fuel commodity",\
                      "doc":"Fuel name accepted by this reactor",\
                      "uitype":"incommodity","uilabel":"Fuel Commodity"}
  std::string fresh_fuel;
  // Add multiple fuel options later

  #pragma cyclus var {"tooltip":"Input fuel recipe",\
                      "doc":"Fuel recipe accepted by reactor",\
                      "uitype":"recipe","uilabel":"Fuel Recipe"}
  std::string fuel_recipe;                      

  #pragma cyclus var {"tooltip":"Spent fuel commodity",\
                      "doc":"Name of spent fuel commodity",\
                      "uitype":"outcommodity","uilabel":"Spent Fuel Commodity"}
  std::string spent_fuel;

  #pragma cyclus var {"tooltip":"Reactor power name",\
                      "doc":"Name of commodity reactor produces",\
                      "uilabel":"Power Name"}
  std::string power_name;

  #pragma cyclus var {"tooltip":"Thermal Power capacity",\
                      "doc":"Reactor thermal power capacity (MW)",\
                      "units":"MW",\
                      "uilabel":"Power Capacity"}
  double power_cap;

 //////////// inventory and core params ////////////
  #pragma cyclus var { \
    "doc": "Mass (kg) of a single assembly.",   \
    "uilabel": "Assembly Mass", \
    "units": "kg", \
  }
  double assem_size;

  #pragma cyclus var { \
    "uilabel": "Number of Assemblies per Batch", \
    "doc": "Number of assemblies that constitute a single batch.  " \
           "This is the number of assemblies discharged from the core fully " \
           "burned each cycle."                                         \
           "Batch size is equivalent to ``n_assem_batch / n_assem_core``.", \
  }
  int n_assem_batch;

  #pragma cyclus var { \
    "uilabel": "Number of Assemblies in Core", \
    "doc": "Number of assemblies that constitute a full core.", \
  }
  int n_assem_core;

  #pragma cyclus var { \
    "default": 0, \
    "uilabel": "Minimum Fresh Fuel Inventory", \
    "units": "assemblies", \
    "doc": "Number of fresh fuel assemblies to keep on-hand if possible.", \
  }
  int n_assem_fresh;

  #pragma cyclus var { \
    "default": 1000000000, \
    "uilabel": "Maximum Spent Fuel Inventory", \
    "units": "assemblies", \
    "doc": "Number of spent fuel assemblies that can be stored on-site before" \
           " reactor operation stalls.", \
  }
  int n_assem_spent;


  ///////// cycle params ///////////
  #pragma cyclus var { \
     "doc": "The duration of a full operational cycle (excluding refueling " \
            "time) in time steps.", \
     "uilabel": "Cycle Length", \
     "units": "time steps", \
  }
  int cycle_time;

  #pragma cyclus var { \
     "doc": "The duration of a full refueling period - the minimum time between"\
            " the end of a cycle and the start of the next cycle.", \
     "uilabel": "Refueling Outage Duration", \
     "units": "time steps", \
  }
  int refuel_time;

  #pragma cyclus var { \
     "default": 0, \
     "doc": "Number of time steps since the start of the last cycle." \
            " Only set this if you know what you are doing", \
     "uilabel": "Time Since Start of Last Cycle", \
     "units": "time steps", \
  }
  int cycle_step;

/*
  #pragma cyclus var {'default': 0.90,\
                      "tooltip":"Capacity factor (%)",\
                      "doc":"Reactor capacity factor; governs downtime between cycles",\
                      "uilabel":"Capacity Factor"}
  double cap_factor; // default to 90% capacity factor
*/
  #pragma cyclus var {'default': 480,\
                      "tooltip":"Reactor lifetime",\
                      "doc":"Reactor lifetime",\
                      "uilabel":"Reactor Lifetime",\
                      "units":"time steps"}
  int reactor_lifetime; 

  #pragma cyclus var {"tooltip":"Path to ORIGEN Libraries",\
                      "doc":"Path to ORIGEN Libraries",\
                      "uilabel":"Path to ORIGEN Libraries"}
  // Default set by environment / build argument in orglib_default_location.h
  std::string lib_path = ORIGEN_LIBS_DEFAULT; 
  
  #pragma cyclus var {"tooltip":"Fresh fuel enrichment",\
                      "doc":"Fresh fuel enrichment",\
                      "uilabel":"Fresh Fuel Enrichment"}
  double enrichment;

  /// Level 2 Parameters

  #pragma cyclus var {'default':"w17x17",\
                      'tooltip':"Assembly type",\
                      'doc':"ORIGEN library type to be used",\
                      'uilabel':"Assembly Type",\
                      'userlevel':1,\
                      'categorical':['ce_facility14x14', 'ce16x16', 'w14x14', 's14x14', 'w15x15', \
                      'w17x17', 'w17x17_ofa', 'ge7x7-0', 'ge8x8-4', 'abb8x8-1', 'ge9x9-7', 'ge10x10-8', \
                      'atrium9-9', 'atrium10-9', 'svea64-1', 'svea100-0']} 
  std::string assembly_type;

// 'default':0.72,
  #pragma cyclus var {"units":"g/cc",\
                      "tooltip":"Moderator Density",\
                      "doc":"Reactor moderator density",\
                      "uilabel":"Moderator Density",\
                      "userlevel":1}
  double mod_density;

  /// Level 3 Parameters

  // should be hidden in ui (internal only). 
  // True if fuel has already been discharged this cycle.
  #pragma cyclus var {"default": 0, "doc": "This should NEVER be set manually",\
                      "internal": True \
  }
  bool discharged;

  // should be hidden in ui (internal only). 
  // True if conditions to perturb output recipe (burnup, initial composition, etc.) have changed
  #pragma cyclus var {"default": 0, "doc": "This should NEVER be set manually",\
                      "internal": True \
  }
  bool refreshSpentRecipe;


  /// Material Flow Parameters                                                      
                                                                     
  // referenced (e.g. n_batch_fresh, assem_size, etc.).
  #pragma cyclus var {"capacity": "n_assem_fresh * assem_size",\
                      "tooltip": "Incoming fresh fuel buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> fresh;
  #pragma cyclus var {"capacity":"n_assem_core * assem_size",\
                      "tooltip":"Fuel held in core"}
  cyclus::toolkit::ResBuf<cyclus::Material> core;
  #pragma cyclus var {"capacity": "n_assem_spent * assem_size",\
                      "tooltip":"Discharged fuel holding capacity / output buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> spent;
  
  // Cached composition of spent fuel (calculated from ORIGEN)
  cyclus::Composition::Ptr spentFuelComp;
                                                                   
  //// A policy for requesting material
  cyclus::toolkit::MatlBuyPolicy buy_policy;
  //// A policy for sending material
  cyclus::toolkit::MatlSellPolicy sell_policy;

  /// Allow test access to private members
  friend class ReactorTest;
  
};


}  // namespace cyborg

#endif  // CYCLUS_CYBORG_REACTOR_H_
