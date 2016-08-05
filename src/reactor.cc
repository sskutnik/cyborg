#include "reactor.h"
#include "orglib_default_location.h"
#include "error.h"
#include <math.h>

//using cyclus::Composition;
namespace cyborg {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// pragmas

#pragma cyclus def schema cyborg::reactor

#pragma cyclus def annotations cyborg::reactor

#pragma cyclus def initinv cyborg::reactor

#pragma cyclus def snapshotinv cyborg::reactor

#pragma cyclus def infiletodb cyborg::reactor

#pragma cyclus def snapshot cyborg::reactor

#pragma cyclus def clone cyborg::reactor

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::InitFrom(reactor* m) {
#pragma cyclus impl initfromcopy cyborg::reactor
 cyclus::toolkit::CommodityProducer::Copy(m);
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::InitFrom(cyclus::QueryableBackend* b) {
  #pragma cyclus impl initfromdb cyborg::reactor

  namespace tk = cyclus::toolkit;
  tk::CommodityProducer::Add(tk::Commodity(power_name),
                             tk::CommodInfo(power_cap, power_cap));


}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reactor::reactor(cyclus::Context* ctx) : cyclus::Facility(ctx), decom(false), 
                                         power_cap(0.0), assem_size(0.0),
                                         n_assem_batch(0), n_assem_core(0),
                                         lib_path(ORIGEN_LIBS_DEFAULT),
                                         spent_fuel("spent_fuel")  {
  cyclus::Warn<cyclus::EXPERIMENTAL_WARNING>("The CyBORG reactor is highly experimental.");
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string reactor::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::EnterNotify(){

  cyclus::Facility::EnterNotify();

  // Check for consistency in fuel inputs first
  size_t n = fuel_incommods.size();
  std::stringstream ss;
  if(fuel_prefs.size() > 0 && fuel_prefs.size() != n) {
     ss << "cyborg::reactor has " << fuel_prefs.size() 
        << " fuel preferences, expected " << n << "\n";
  }
  if(fuel_recipes.size() != n) {
     ss << "cyborg::reactor has " << fuel_recipes.size() 
        << " input recipes, expected " << n << "\n";
  }

  // Check for consistency of batch power fractions with number of batches
  size_t n_batch = n_assem_core / n_assem_batch;
  if ( core_power_frac.size() > 0) {
     if( core_power_frac.size() == n_batch ) {
        // Normalize sum of batch power fractions to 1.0
        double powNorm = std::accumulate(core_power_frac.begin(),core_power_frac.end(), 0.0);
        std::transform(core_power_frac.begin(), core_power_frac.end(), core_power_frac.begin(), 
                   std::bind1st(std::multiplies<double>(),1.0/powNorm));    
     }
     else {
        ss << "cyborg::reactor has " << core_power_frac.size() 
           << " values for batch power fraction, however " 
           << " # of batches (n_assem_core / n_assem_batch) = " 
           << n_batch << "\n";
     }
 
  }
  else {
     core_power_frac.resize(n_batch, 1.0/static_cast<double>(n_batch));
  }

  // Handle any input errors
  if(ss.str().size() > 0) {
     throw cyclus::ValueError(ss.str());
  }

  // Initialize fuel_prefs to default values if not specified by user
  if (fuel_prefs.size() == 0) {
     for (size_t i=0; i < fuel_incommods.size(); ++i) {
        this->fuel_prefs.push_back(cyclus::kDefaultPref);
     }
  }

  cyclus::CompMap v;
  cyclus::Composition::Ptr comp, nullComp;
  // dummy comp, use in_recipe if provided
  nullComp = cyclus::Composition::CreateFromAtom(v);

  //TODO: There has to be a way to load trades directly into the core, right?
  buy_policy.Init(this, &fresh, "fresh fuel", this->fuel_capacity(), 1.0, 1.0, this->assem_size);
  for(size_t i=0; i < fuel_incommods.size(); ++i) {
     comp = nullComp; 
     if (fuel_recipes[i] != "") {
        comp = context()->GetRecipe(fuel_recipes[i]); 
     }   
     buy_policy.Set(fuel_incommods[i], comp, fuel_prefs[i]);
  }
  buy_policy.Start();   

  sell_policy.Init(this, &spent, spent_fuel);
  sell_policy.Set(spent_fuel);
  sell_policy.Start();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tick() {
   if(retired()) {
     Record("RETIRED", "");
    
     // Record the last time series entry if the reactor was operating at 
     // the time of retirement
     if( context()->time() == exit_time() ) {
       if( cycle_step > 0 && cycle_step <= cycle_time &&
           core.count() == n_assem_core) {
         cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, power_cap); 
       }
       else {
         cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, 0.0); 
       }

       // Handle transmutation of each batch in the last core
       // i.e., oldest fuel transmuted to "full" burnup,
       // next-oldest cut short by one cycle, etc.

       // Account for potentially abbreviated final cycle length
       double last_cycle = static_cast<double>(cycle_step) / static_cast<double>(cycle_time);

       for(size_t i = core_power_frac.size(); i > 0; --i) {
          Transmute_(n_assem_batch, i, last_cycle);
       }       
     }
     // Dump remaining fresh inventory into spent fuel to be traded away
     while(fresh.count() > 0 && spent.space() >= assem_size) {
        spent.Push(fresh.Pop());
     }     
     if(fresh.count() == 0) fresh.capacity(0); 

     // Attempt to discharge all transmuted assemblies from the core
     discharged = Discharge_(n_assem_core);     
     return;
   } // end retired() check

   // Transmute & discharge if necessary     
   if (cycle_step == cycle_time) {
     Transmute_();
     Record("CYCLE_END", "");
   }
   if(cycle_step >= cycle_time && !discharged) {
     discharged = Discharge_();
     //std::cerr << "Called discharge: result = " << discharged << "  core.count() = " << core.count() << std::endl;
   }
   if(discharged) { 
     // Only load core once we've fully cleared out the fully-burnt assemblies
     Load_();
   }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tock() {
  if (retired()) {
    return;
  }

  if (cycle_step >= cycle_time + refuel_time || cycle_step == 0) {   
    if( core.count() == n_assem_core) {
      // Restart once the core is fully reloaded
      discharged = false;
      cycle_step = 0;
    }    
  }

  if (cycle_step == 0 && core.count() == n_assem_core) {
    Record("CYCLE_START", "");
  }

  if (cycle_step >= 0 && cycle_step < cycle_time &&
      core.count() == n_assem_core) {
    cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, power_cap);
  } else {
    cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, 0);
  }

  // "if" prevents starting cycle after initial deployment until core is full
  // even though cycle_step is its initial zero.
  if (cycle_step > 0 || core.count() == n_assem_core) {
    ++cycle_step;
  }
}


void reactor::Load_() {
  int n = std::min(n_assem_core - core.count(), fresh.count());
/*
  std::cerr << "Attempting to load: " << n << " assemblies." << std::endl
            << "fresh.count() = " << fresh.count() << "  core.count = " << core.count() << std::endl
            << "n_assem_core = " << n_assem_core << "  n_assem_batch = " << n_assem_batch << std::endl;
*/
  if (n == 0) {
    return;
  }

  std::stringstream ss;
  ss << n << " assemblies";
  Record("LOAD", ss.str());
  core.Push(fresh.PopN(n));

}

bool reactor::Discharge_() { return Discharge_(this->n_assem_batch); }

bool reactor::Discharge_(int n_assem_discharged) {

  int npop = std::min(n_assem_discharged, core.count());
  if (n_assem_spent - spent.count() < npop) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }

  std::stringstream ss;
  ss << npop << " assemblies";
  Record("DISCHARGE", ss.str());
      
  // Discharge fuel to spent fuel buffer
  spent.Push(core.PopN(n_assem_discharged));
  return true;
}

void reactor::Transmute_() { Transmute_(n_assem_batch); }

void reactor::Transmute_(int n_assem, int n_cycles, double last_cycle) {
  using cyclus::toolkit::MatVec;
 
  if(n_cycles == -1) {
     n_cycles = round(this->n_assem_core / this->n_assem_batch);
  }

  // Instead of doing a PopN for all assemblies, peek at comps and pop until we've hit the right # of assemblies?
  MatVec old = core.PopN(std::min(n_assem, core.count()));

 /*
 * TODO: Handle multiple depletion recipes, assuming non-homogeneous batches?
 */

  // Check that all assemblies in the batch have the same composition
  // TODO make this a loop until we've traversed to the end of the batch?

  std::vector<int> matIDs(old.size());
/*
     for(auto & mat : old) matIDs.push_back(mat->comp()->id());
     auto idx = std::adjacent_find( matIDs.begin(), matIDs.end(), std::not_equal_to<int>() );
     if(idx != old.end()) {
       // Split batch at idx, re-evaluate
     }   
*/     
  bool isHomogenous = (std::adjacent_find( matIDs.begin(), matIDs.end(), std::not_equal_to<int>() ) == matIDs.end());
  // FOR NOW: Assume everything in the batch is the same composition
  spentFuelComp = this->Deplete_(old[0],n_assem_batch, n_cycles, last_cycle);

  if(!spentFuelComp) {
     throw cyclus::StateError("Spent fuel composition is not set!");
  } 

  for (int i = 0; i < old.size(); i++) {
     // TODO: Add in multiple spent_comps, map assemblies to compositions each time we do a depletion...
     old[i]->Transmute(spentFuelComp);
  }

  // Push transmuted materials back to the bottom of the core buffer
  core.Push(old);
  if (core.count() > old.size()) {
    // rotate untransmuted mats back to back of buffer
    core.Push(core.PopN(core.count() - old.size()));
  }
   

  std::stringstream ss;
  ss << old.size() << " assemblies" << " to discharge burnup " << this->burnup << " MWd/MTHM";
  Record("TRANSMUTE", ss.str());
}


cyclus::Composition::Ptr reactor::Deplete_(cyclus::Material::Ptr mat, const int n_assem, const int n_cycles, const double last_cycle) {
   
    OrigenInterface::cyclus2origen react;
     
    // Set ORIGEN library path
    react.set_lib_path(lib_path);
    
    // Set ID tags
    if(this->assembly_type == "") {
       std::stringstream ss;
       ss << "Cyborg::reactor::Deplete_() - assembly_type unspecified!" << std::endl; 
       throw cyclus::StateError(ss.str());
    }    
    react.set_id_tag("Assembly Type",assembly_type);
        
    this->setup_origen_interp_params(react, mat);
    this->setup_origen_power_history(react, n_cycles, last_cycle);
/*
    std::cerr << "ID TAGS\n=======" << std::endl;
    react.list_id_tags();
    std::cerr << "\nINTERP PARAMS\n=======" << std::endl;
    react.list_parameters();
*/

/*
 *  Note: Do we need to consider non-interpolated parts of the recipe too? (e.g., Pu-240 content, U-234 content...?) 
 *        Store these as ID tags as well to trigger a recipe update? 
 */            

    // If this input state (input composition, library, powers, cycle lengths) 
    // has been previously depleted in the Cyclus context, retrieved the cached
    // recipe. Otherwise, generate the new recipe and store it. 
    //
    // The input state is "hashed" via the TagManager on the Cyclus/ORIGEN interface layer

    std::string depletion_state = react.get_tag_manager_string();
    cyclus::Composition::Ptr comp_out = NULL;
    try {
       comp_out = context()->GetRecipe(depletion_state); 
       if(!comp_out) {
         throw cyclus::KeyError("Empty recipe returned.");
       }
       else {
          //TODO: Need to manually update burnup via calculation if we're using a cached recipe...
          return comp_out;
       }
    }
    catch(cyclus::KeyError ke) {      
       cyclus::Warn<cyclus::KEY_WARNING>("Recipe lookup failed: " + std::string(ke.what()));
    }

    // Create cross-section library
    react.interpolate();

    // Set input compositions to deplete
    this->setup_origen_materials(react, mat, n_assem);
    
    // Run ORIGEN depletion calculation
    react.solve();

    this->burnup = react.burnup_last();  // Burnup in units of MWd/MTHM

    // Update recipe based on discharge compositions from ORIGEN
    cyclus::CompMap dischargeRecipe = this->get_origen_discharge_recipe(react);
    comp_out = cyclus::Composition::CreateFromAtom(dischargeRecipe);    
    if(!comp_out) {
       std::stringstream ss;
       ss << "Error in Deplete_(); could not create discharge fuel recipe!\n";
       throw cyclus::StateError(ss.str());
    }
 
    // Store hashed discharge recipe in Cyclus global context for future reuse
    context()->AddRecipe(depletion_state, comp_out);

    return comp_out;
}

void reactor::setup_origen_interp_params(OrigenInterface::cyclus2origen& react, const cyclus::Material::Ptr mat) {
    // Set Interpolable parameters   
    if(boost::to_upper_copy(this->fuel_type) == "UOX") {
       double enrich = get_iso_mass_frac(92, 235, mat->comp()) * 100.0;
       if(enrich <= 0.0 || enrich > 100.0) {         
          std::stringstream ss;
          ss << "Cyborg::reactor::Deplete_(); invalid U-235 enrichment!"
             << " Calculated enrichment = " << enrich << "\n";
          throw cyclus::ValueError(ss.str());
       }      
       react.add_parameter("Enrichment",enrich);
    }
    else if(boost::to_upper_copy(this->fuel_type) == "MOX") {
       double fr_pu239 = get_iso_mass_frac(94, 239, mat->comp()) * 100.0;
       double fr_pu = get_ele_hm_mass_frac(94, mat->comp()) * 100.0;

       if(fr_pu239 <= 0.0 || fr_pu239 > 100.0) {         
          std::stringstream ss;
          ss << "Cyborg::reactor::Deplete_(); invalid Pu-239 enrichment!"
             << " Calculated enrichment = " << fr_pu239 << "\n";
          throw cyclus::ValueError(ss.str());
       }             
       if(fr_pu <= 0.0 || fr_pu > 100.0) {         
          std::stringstream ss;
          ss << "Cyborg::reactor::Deplete_(); invalid Pu heavy metal fraction!"
             << " Calculated Pu fraction = " << fr_pu << "\n";
          throw cyclus::ValueError(ss.str());
       }      
       react.add_parameter("Plutonium Content",fr_pu);
       react.add_parameter("Plutonium-239 Content",fr_pu239);
    }
   
    //TODO: Do we do any tag vetting here, or just let ORIGEN do it?
    for(auto &tag : interp_tags) {
       react.add_parameter(tag.first,tag.second);
    } 
    //if(this->mod_density > 0.0) react.add_parameter("Moderator Density",this->mod_density);   
}

void reactor::setup_origen_power_history(OrigenInterface::cyclus2origen& react, const int n_cycles, const double last_cycle) {

    std::vector<double> dp_time, dp_pow;
    dp_time.push_back(0.0);
    // Number of cycles is by default the number of batches - i.e., 3 batches => 3 cycles
    // However, it is user-configurable to allow for partially-burnt assemblies 
    // (i.e., for reactor decommissioning behavior)
    for(size_t i=0; i < n_cycles; ++i) {
       // Cycle timestep is in months; use years for ORIGEN for simplicity
       double time_tmp = static_cast<double>(cycle_time)/12.0;
       if(i == (n_cycles - 1)) time_tmp *= last_cycle;
       dp_time.push_back(time_tmp + dp_time.back());

       // Convert power to MWt 
       double cyclePower = power_cap * core_power_frac[i] * 1E6;
       dp_pow.push_back(cyclePower);
       
       // Decay fuel during reload
       if(refuel_time > 0) {
          dp_time.push_back(dp_time.back() + static_cast<double>(refuel_time)/12.0);
          dp_pow.push_back(0.0);
       }
    }
    react.set_time_units("y");
    react.set_time_steps(dp_time); 
    
    react.set_power_units("W");
    react.set_powers(dp_pow);  
}

void reactor::setup_origen_materials(OrigenInterface::cyclus2origen& react, const cyclus::Material::Ptr mat, const int n_assem) {
    // Pass nuclide IDs and masses to ORIGEN
    std::vector<int> in_ids;
    std::vector<double> mass_fraction;

    // Get fuel recipe and convert to ORIGEN format
    cyclus::CompMap comp_in = mat->comp()->mass();
    for(std::map<int,double>::iterator it = comp_in.begin(); it!=comp_in.end(); it++){
       int id = it->first;
       mass_fraction.push_back(it->second);
       // convert id to ORIGEN ZAID format 
       in_ids.push_back(pyne::nucname::zzaaam(id));
    }

    // normalize mass fractions
    std::vector<double> norm_mass(mass_fraction.size());

    // Convert each isotopic mass to absolute isotopic mass in kg (from unnormalized mass fraction)
    double massNorm = std::accumulate(mass_fraction.begin(),mass_fraction.end(), 0.0);    
    std::transform(mass_fraction.begin(), mass_fraction.end(), norm_mass.begin(), 
                   std::bind1st(std::multiplies<double>(),mat->quantity()/massNorm * n_assem));
    //for(auto mass : norm_mass) { std::cerr << "Normed mass: " << mass/(mat->quantity()) <<  std::endl; }
    react.set_materials(in_ids,norm_mass);
}

// Get materials and convert nuclide ids back to Cyclus format
cyclus::CompMap reactor::get_origen_discharge_recipe(OrigenInterface::cyclus2origen& react) {

    std::vector<int> org_id;
    react.get_ids_zzzaaai(org_id);
     
    std::for_each(org_id.begin(), org_id.end(), [](int &nucID){ pyne::nucname::id(nucID); });
   
    // Get mass data from ORIGEN
    std::vector<double> org_atom;
    react.get_masses_final(org_atom,"GATOMS");

    // Normalize to atom fractions
    double atomNorm = std::accumulate(org_atom.begin(),org_atom.end(), 0.0);    
    std::transform(org_atom.begin(), org_atom.end(), org_atom.begin(), 
                   std::bind1st(std::multiplies<double>(),1.0/atomNorm));
    
    cyclus::CompMap v;
    for(int j=0; j!=org_id.size(); ++j){       
       if(org_atom[j] > 0.) { 
          v[org_id[j]] = org_atom[j];
          //if(org_atom[j] > 1.E-4) std::cerr << "Setting v[" << org_id[j] << "] to: " << org_atom[j] << std::endl;
       }
    }
    return v;
}

void reactor::Record(std::string name, std::string val) {
  context()
      ->NewDatum("ReactorEvents")
      ->AddVal("AgentId", id())
      ->AddVal("Time", context()->time())
      ->AddVal("Event", name)
      ->AddVal("Value", val)
      ->Record();
}

double get_ele_hm_mass_frac(const int Z, const cyclus::Composition::Ptr comp, const int Z_HM) {
   cyclus::CompMap mass_map = comp->mass();
   double hm_mass = 0.0;
   double ele_mass = 0.0;

   for(auto &nuc : mass_map) {
     if( static_cast<int>(floor(nuc.first / 1E7)) >= Z_HM ) hm_mass += nuc.second;     
     if( static_cast<int>((nuc.first / 1E7)) == Z) ele_mass += nuc.second;
   }
   if (hm_mass <= 0.0) return -1;
   return (ele_mass / hm_mass);
}

double get_iso_mass_frac(const int Z, const int A, const cyclus::Composition::Ptr comp) {
   cyclus::CompMap mass_map = comp->mass();
   double ele_mass = 0.0;
   double iso_mass = 0.0;

   for(auto &nuc : mass_map) {
     if( static_cast<int>(floor(pyne::nucname::id(nuc.first) / 1E7)) == Z ) {
       ele_mass += nuc.second;
       if( static_cast<int>((pyne::nucname::id(nuc.first) - Z*1E7)/1E4) == A) iso_mass += nuc.second;
     }  
   }
   if(ele_mass <= 0.0) return -1;
   return (iso_mass / ele_mass);   
}


extern "C" cyclus::Agent* Constructreactor(cyclus::Context* ctx) {
  return new reactor(ctx);
}

}  // namespace cyborg
