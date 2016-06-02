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

  void Discharge_(double);

  /// Deplete the material for this cycle
  /// @param mat   The material to be depleted 
  /// @param power Depletion power of the material in MW
  cyclus::Material::Ptr Deplete_(cyclus::Material::Ptr, double);

  int reactor_time;
  bool decom;
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

  #pragma cyclus var {"tooltip":"Fuel capacity",\
                      "doc":"Total reactor fuel capacity (MT)",\
                      "units":"MT",\
                      "uilabel":"Fuel Capacity"}
  double fuel_capacity;

  #pragma cyclus var {"tooltip":"Cycle length",\
                      "doc":"Time to complete an entire cycle",\
                      "uilabel":"Cycle Length",\
                      "units":"time steps"}
  int cycle_length;

  #pragma cyclus var {"tooltip":"Capacity factor",\
                      "doc":"Reactor capacity factor",\
                      "uilabel":"Capacity Factor"}
  double cap_factor; 

  #pragma cyclus var {"tooltip":"Reactor lifetime",\
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

  #pragma cyclus var {'default':0.72,\
                      'units':'g/cc',\
                      "tooltip":"Moderator Density",\
                      "doc":"Reactor moderator density",\
                      "uilabel":"Moderator Density",\
                      "userlevel":1}
  double mod_density;

  #pragma cyclus var {'default':0,\
                      "tooltip":"Burnup",\
                      "doc":"Reactor burnup (MWd/tHM)",\
                      "uilabel":"Burnup",\
                      "userlevel":1}
  double burnup;

  /// Level 3 Parameters


  /// Material Flow Parameters                                                      
                                                                     
  #pragma cyclus var {"tooltip":"Incoming material buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> fresh_inventory;

  #pragma cyclus var {"tooltip":"Fuel buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> fuel;

  #pragma cyclus var {"tooltip":"Outgoing material buffer"}
  cyclus::toolkit::ResBuf<cyclus::Material> spent_inventory;

  //// A policy for requesting material
  cyclus::toolkit::MatlBuyPolicy buy_policy;
  //// A policy for sending material
  cyclus::toolkit::MatlSellPolicy sell_policy;
  
};

}  // namespace cyborg

#endif  // CYCLUS_CYBORG_REACTOR_H_
