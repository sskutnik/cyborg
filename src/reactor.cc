#include "reactor.h"
#include "orglib_default_location.h"
#include "cyclus_origen_interface.h"
#include "error.h"
#include "composition.h"
#include <math.h>

//using cyclus::Composition;
namespace cyborg {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// pragmas

#pragma cyclus def schema cyborg::Reactor

#pragma cyclus def annotations cyborg::Reactor

#pragma cyclus def initinv cyborg::Reactor

#pragma cyclus def snapshotinv cyborg::Reactor

#pragma cyclus def infiletodb cyborg::Reactor

#pragma cyclus def snapshot cyborg::Reactor

#pragma cyclus def clone cyborg::Reactor

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Reactor::InitFrom(Reactor* m) {
#pragma cyclus impl initfromcopy cyborg::Reactor
 cyclus::toolkit::CommodityProducer::Copy(m);
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Reactor::InitFrom(cyclus::QueryableBackend* b) {
  #pragma cyclus impl initfromdb cyborg::Reactor

  namespace tk = cyclus::toolkit;
  tk::CommodityProducer::Add(tk::Commodity(power_name),
                             tk::CommodInfo(power_cap, power_cap));


}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Reactor::Reactor(cyclus::Context* ctx) : cyclus::Facility(ctx), decom(false), 
                                         power_cap(0.0), assem_size(0.0),
                                         n_assem_batch(0), n_assem_core(0),
                                         lib_path(ORIGEN_LIBS_DEFAULT),
                                         spent_fuel("spent_fuel")  {
  cyclus::Warn<cyclus::EXPERIMENTAL_WARNING>("The CyBORG Reactor is highly experimental.");
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string Reactor::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Reactor::EnterNotify(){

  cyclus::Facility::EnterNotify();

  // Check for consistency in fuel inputs first
  size_t n = fuel_incommods.size();
  std::stringstream ss;
  if(fuel_prefs.size() > 0 && fuel_prefs.size() != n) {
     ss << "cyborg::Reactor has " << fuel_prefs.size() 
        << " fuel preferences, expected " << n << "\n";
  }
  if(fuel_recipes.size() != n) {
     ss << "cyborg::Reactor has " << fuel_recipes.size() 
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
        ss << "cyborg::Reactor has " << core_power_frac.size() 
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

  // TODO: Add a spent fuel throughput parameter
  sell_policy.Init(this, &spent, spent_fuel, std::numeric_limits<double>::max(), false, this->assem_size);
  sell_policy.Set(spent_fuel);
  sell_policy.Start();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Reactor::Tick() {
   if(retired()) {
     Record("RETIRED", "");
    
     // Record the last time series entry if the Reactor was operating at 
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
          // Need to manually rotate assemblies to ensure we deplete each batch 
          // in the core. This is due to the fact that we don't "pop" each one 
          // into spent immmediatley after trasmuting, unlike the normal case.
          cyclus::toolkit::MatVec depleted = core.PopN(std::min(n_assem_batch, core.count()));
          core.Push(depleted);
       }       
     }
     
     // Attempt to discharge all transmuted assemblies from the core  
     /*
     while ( core.count() > 0) {
        if(!Discharge_()) break;
     }     
     if(core.count() == 0) discharged = true;
     */
     if(core.count() > 0) discharged = Discharge_(n_assem_core);     

     // Dump remaining fresh inventory into spent fuel to be traded away
     while(fresh.count() > 0 && spent.space() >= assem_size) {
        spent.Push(fresh.Pop());
     }     

     return;
   } // end retired() check

   // Transmute & discharge if necessary     
   if (cycle_step == cycle_time) {
     Transmute_();
     Record("CYCLE_END", "");
   }
   if((cycle_step >= cycle_time) && !discharged) {
     discharged = Discharge_();
   }
   if(discharged) { 
     // Only load core once we've fully cleared out the fully-burnt assemblies
     Load_();
   }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Reactor::Tock() {

  if (retired()) {
    // Can retire sell policy once the last remaining assemblies are traded out
    if(spent.count() == 0) sell_policy.Stop();
    return;
  }

  if (cycle_step >= cycle_time + refuel_time || cycle_step == 0) {   
    // Make sure we actually discharged all of the burnt fuel; i.e., the core is 
    // "full" even when we fail to discharge (due to lack of space)
    if( discharged && core.count() == n_assem_core) {
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


void Reactor::Load_() {
  int n = std::min(n_assem_core - core.count(), fresh.count());
  if (n == 0) return;
  
  std::stringstream ss;
  ss << n << " assemblies";
  Record("LOAD", ss.str());
  core.Push(fresh.PopN(n));

}

bool Reactor::Discharge_() { return Discharge_(this->n_assem_batch); }

bool Reactor::Discharge_(int n_assem_discharged) {
  
  int npop = std::min(n_assem_discharged, core.count());
  if (n_assem_spent - spent.count() < npop) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }

  std::stringstream ss;
  ss << npop << " assemblies";
  Record("DISCHARGE", ss.str());
 
  // Discharge fuel to spent fuel buffer
  spent.Push(core.PopN(npop));

  return true;
}

void Reactor::Transmute_() { Transmute_(n_assem_batch); }

void Reactor::Transmute_(int n_assem, int n_cycles, double last_cycle) {
  using cyclus::toolkit::MatVec;
 
  if(n_cycles == -1) {
     n_cycles = round(this->n_assem_core / this->n_assem_batch);
  }


  MatVec old = core.PopN(std::min(n_assem, core.count()));

  // Check whether all assemblies in the batch have the same composition
  // Note: Need to actually check by mass composition, as composition IDs 
  //       can be different for the same composition...

  //typedef std::reference_wrapper<const cyclus::CompMap> CompMapRef;
  std::vector<cyclus::CompMap> depCompMaps;

  //for(auto & mat : old) matIDs.push_back(mat->comp()->id());  
  for(auto & mat : old)  depCompMaps.push_back(mat->comp()->mass()); 
  
  std::vector<cyclus::CompMap>::iterator idx = depCompMaps.begin();
  std::vector<cyclus::CompMap>::iterator idx_old = idx; 

  int n_subbatch = 0;
  int index = 0;

  while(idx != depCompMaps.end()) {
     idx = std::adjacent_find(idx_old, depCompMaps.end(), std::not_equal_to<cyclus::CompMap>() ); 

     index = std::distance(depCompMaps.begin(), idx);
     if(idx == depCompMaps.end()) {
       // Reached the end of the IDs array; use the last entry
       index = depCompMaps.size() - 1;
     }
     else if(idx == idx_old) {
       // Unique element in difference seqeuence; i.e., 1-2-3. Need to kick iterator along to continue
       idx++;
     }
     n_subbatch = idx - idx_old;
  
     spentFuelComp = this->Deplete_(old[index], n_subbatch, n_cycles, last_cycle); 
     if(!spentFuelComp) {
        throw cyclus::StateError("Spent fuel composition is not set!");
     } 
     
     // Transmute all assemblies in the homogeneous sub-batch
     for(auto it = idx_old; it != idx; it++) {
        old[std::distance(depCompMaps.begin(), it)]->Transmute(spentFuelComp);
     }
     idx_old = idx;
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


cyclus::Composition::Ptr Reactor::Deplete_(cyclus::Material::Ptr mat, const int n_assem, 
                                           const int n_cycles, const double last_cycle) { 

    OrigenInterface::cyclus2origen react;
     
    // Set ORIGEN library path
    react.set_lib_path(lib_path);
    
    // Set ID tags
    if(this->assembly_type == "") {
       std::stringstream ss;
       ss << "Cyborg::Reactor::Deplete_() - assembly_type unspecified!" << std::endl; 
       throw cyclus::StateError(ss.str());
    }    
    react.set_id_tag("Assembly Type",assembly_type);
    
    this->setup_origen_interp_params(react, mat);
    this->setup_origen_power_history(react, n_cycles, n_assem, last_cycle);


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

void Reactor::setup_origen_interp_params(OrigenInterface::cyclus2origen& react, const cyclus::Material::Ptr mat) {
    // Set Interpolable parameters   

    if(boost::to_upper_copy(this->fuel_type) == "UOX") {
       double enrich = get_iso_mass_frac(92, 235, mat->comp()) * 100.0;
       if(enrich <= 0.0 || enrich > 100.0) {         
          std::stringstream ss;
          ss << "Cyborg::Reactor::Deplete_(); invalid U-235 enrichment!"
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
          ss << "Cyborg::Reactor::Deplete_(); invalid Pu-239 enrichment!"
             << " Calculated enrichment = " << fr_pu239 << "\n";
          throw cyclus::ValueError(ss.str());
       }             
       if(fr_pu <= 0.0 || fr_pu > 100.0) {         
          std::stringstream ss;
          ss << "Cyborg::Reactor::Deplete_(); invalid Pu heavy metal fraction!"
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

void Reactor::setup_origen_power_history(OrigenInterface::cyclus2origen& react, const int n_cycles, 
                                         const int n_assem,  const double last_cycle) {

    std::vector<double> dp_time, dp_pow;
    dp_time.push_back(0.0);
    // Number of cycles is by default the number of batches - i.e., 3 batches => 3 cycles
    // However, it is user-configurable to allow for partially-burnt assemblies 
    // (i.e., for Reactor decommissioning behavior)
    for(size_t i=0; i < n_cycles; ++i) {
       // Cycle timestep is in months; use years for ORIGEN for simplicity
       double time_tmp = static_cast<double>(cycle_time)/12.0;
       if(i == (n_cycles - 1)) time_tmp *= last_cycle;
       dp_time.push_back(time_tmp + dp_time.back());

       // Convert power to MWt 
       double cyclePower = power_cap * core_power_frac[i] * n_assem / n_assem_batch * 1E6;
       dp_pow.push_back(cyclePower);

       // Decay fuel during reload; don't include decay from last reload cycle
       if(refuel_time > 0 and i < (n_cycles-1)) {
          dp_time.push_back(dp_time.back() + static_cast<double>(refuel_time)/12.0);
          dp_pow.push_back(0.0);
       }   
    }
    react.set_time_units("y");
    react.set_time_steps(dp_time); 
    
    react.set_power_units("W");
    react.set_powers(dp_pow);  
}

void Reactor::setup_origen_materials(OrigenInterface::cyclus2origen& react, const cyclus::Material::Ptr mat, 
                                     const int n_assem) {
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

    react.set_mat_units("KILOGRAMS");
    react.set_materials(in_ids,norm_mass);
}

// Get materials and convert nuclide ids back to Cyclus format
cyclus::CompMap Reactor::get_origen_discharge_recipe(OrigenInterface::cyclus2origen& react) {
   
    // Get mass data from ORIGEN
    std::map<int,double> org_atom;
    react.get_masses_final_map(org_atom,"GATOMS");   

    cyclus::CompMap v;
    int pyneID;
    for(auto const &nucl : org_atom) {
       if(nucl.second > 0) {
          pyneID = pyne::nucname::zzaaam_to_id(nucl.first);
          v[pyneID] = nucl.second;
          //if(nucl.second > 1E-4) std::cerr << "Setting v[" << pyneID << "] -> " << nucl.second << std::endl;
       }
    }
    // Force normalization
    //cyclus::compmath::Normalize(&v,1.0);

    return v;
}

void Reactor::Record(std::string name, std::string val) {
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

   if(comp->mass().size() == 0) {
     std::stringstream ss;
     ss << "get_iso_mass_frac(): null composition found!\n";
     throw cyclus::StateError(ss.str());
   }

   for(auto &nuc : mass_map) {
     if( static_cast<int>(floor(pyne::nucname::id(nuc.first) / 1E7)) == Z ) {
       ele_mass += nuc.second;
       if( static_cast<int>((pyne::nucname::id(nuc.first) - Z*1E7)/1E4) == A) iso_mass += nuc.second;
     }  
   }
   if(ele_mass <= 0.0) return -1;
   return (iso_mass / ele_mass);   
}

/// Buy policy; manually configured (from cycamore::Reactor) to allow for JIT fuel trading
std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> Reactor::GetMatlRequests() {
  using cyclus::RequestPortfolio;
  using cyclus::Request;

  std::set<RequestPortfolio<cyclus::Material>::Ptr> ports;
  cyclus::Material::Ptr m;

  // second min expression reduces assembles to amount needed until
  // retirement if it is near.
  int n_assem_order = n_assem_core - core.count() + n_assem_fresh - fresh.count();

  if (exit_time() != -1) {
    // the +1 accounts for the fact that the Reactor is alive and gets to
    // operate during its exit_time time step.
    int t_left = exit_time() - context()->time() + 1;
    int t_left_cycle = cycle_time + refuel_time - cycle_step;
    double n_cycles_left = static_cast<double>(t_left - t_left_cycle) /
                         static_cast<double>(cycle_time + refuel_time);
    n_cycles_left = ceil(n_cycles_left);
    int n_need = std::max(0.0, n_cycles_left * n_assem_batch - n_assem_fresh + n_assem_core - core.count());
    n_assem_order = std::min(n_assem_order, n_need);
  }

  if (n_assem_order == 0 || this->retired()) {
    return ports;
  }

  for (int i = 0; i < n_assem_order; i++) {
    RequestPortfolio<cyclus::Material>::Ptr port(new RequestPortfolio<cyclus::Material>());
    std::vector<Request<cyclus::Material>*> mreqs;
    for (int j = 0; j < fuel_incommods.size(); j++) {
      std::string commod = fuel_incommods[j];
      double pref = fuel_prefs[j];
      cyclus::Composition::Ptr recipe = context()->GetRecipe(fuel_recipes[j]);

      m = cyclus::Material::CreateUntracked(assem_size, recipe);
      Request<cyclus::Material>* r = port->AddRequest(m, this, commod, pref, true);
      mreqs.push_back(r);
    }
    port->AddMutualReqs(mreqs);
    ports.insert(port);
  }

  return ports;
}

void Reactor::AcceptMatlTrades(const std::vector<
    std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses) {
  std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                        cyclus::Material::Ptr> >::const_iterator trade;

  std::stringstream ss;
  int nload = std::min((int)responses.size(), n_assem_core - core.count());
  if (nload > 0) {
    ss << nload << " assemblies";
    Record("LOAD", ss.str());
  }

  for (trade = responses.begin(); trade != responses.end(); ++trade) {
    std::string commod = trade->first.request->commodity();
    cyclus::Material::Ptr m = trade->second;
    index_res(m, commod);

    if (core.count() < n_assem_core) {
      core.Push(m);
    } else {
      fresh.Push(m);
    }
  }
}

void Reactor::index_res(cyclus::Resource::Ptr m, std::string incommod) {
  for (int i = 0; i < fuel_incommods.size(); i++) {
    if (fuel_incommods[i] == incommod) {
      res_indexes[m->obj_id()] = i;
      return;
    }
  }
  throw cyclus::ValueError("cyborg::Reactor - received unsupported incommod material");
}


extern "C" cyclus::Agent* ConstructReactor(cyclus::Context* ctx) {
  return new Reactor(ctx);
}

}  // namespace cyborg
