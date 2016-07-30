#include "reactor.h"
#include "cyclus_origen_interface.h"
#include "error.h"

//using cyclus::Composition;
namespace cyborg {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reactor::reactor(cyclus::Context* ctx) : cyclus::Facility(ctx), decom(false), 
                                         discharged(false), refreshSpentRecipe(true),
                                         power_cap(0.0), assem_size(0.0),
                                         n_assem_batch(0), n_assem_core(0),
                                         enrichment(0.0),
                                         fresh_fuel("fresh_fuel"), spent_fuel("spent_fuel")  {
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string reactor::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::EnterNotify(){
    Facility::EnterNotify();
    buy_policy.Init(this, &fresh, fresh_fuel);
    buy_policy.Set(fresh_fuel).Start();

    sell_policy.Init(this, &spent, spent_fuel);
    sell_policy.Set(spent_fuel).Start();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tick() {
  std::cerr << "Calling Tick(): decom = " << decom << std::endl;
  if (!decom) {
    // Transmute & Discharge if necessary     
    //std::cerr << "cycle_time = " << cycle_time << "  cycle_length = " << cycle_length << std::endl;
    if (cycle_step == cycle_time) {
      std::cerr << "Calling transmute" << std::endl;
      Transmute_();
    }
    if(cycle_step >= cycle_time && !discharged) {
      discharged = Discharge_();
      //std::cerr << "Called discharge: result = " << discharged << std::endl;
    }
  }
  else {
    // Should ideally push out all fresh fuel first...
    fresh.capacity(0);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tock() {
  if (retired()) {
    return;
  }

  if (cycle_step >= cycle_time + refuel_time) {   
    if( core.count() == n_assem_core) {
      // Restart once the core is fully reloaded
      discharged = false;
      cycle_step = 0;
    }
    else if(discharged) { 
      // Only load core once we've fully cleared out the fully-burnt assemblies
      Load_();
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
  
  //std::cerr << "Finished Tock()" << std::endl;
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
  //std::cerr << "LOADING " << ss.str() << std::endl;
  core.Push(fresh.PopN(n));
  //std::cerr << "PUSHED " << ss.str() << std::endl;
}

bool reactor::Discharge_() { return Discharge_(this->n_assem_batch); }

bool reactor::Discharge_(int n_assem_discharged) {

  int npop = std::min(n_assem_batch, core.count());
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

void reactor::Transmute_(int n_assem) {
  using cyclus::toolkit::MatVec;

  MatVec old = core.PopN(std::min(n_assem, core.count()));
  core.Push(old);
  if (core.count() > old.size()) {
    // rotate untransmuted mats back to back of buffer
    core.Push(core.PopN(core.count() - old.size()));
  }
 
 /*
 * TODO: Call recipe update if needed; otherwise, use old recipe
 * TODO: Alternative: Generate new recipe every time new fuel / power / etc. conditions
 *       come about, push this onto the stack?
 * TODO: Handle multiple depletion recipes, assuming non-homogeneous batches?
 */
  
  if(refreshSpentRecipe) {
     //TODO: Handle cycle powers by batch, send in powers vector

     // Calculate total thermal power of depleted assemblies; 
     // Convert from MWt => W
     double cyclePower = power_cap * static_cast<double>(n_assem_batch)/static_cast<double>(n_assem_core) * 1E6;

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
     spentFuelComp = this->Deplete_(old[0],n_assem_batch,cyclePower);
     refreshSpentRecipe = false;  //TODO: Add code to check when this needs to be turned back on
  }
  for (int i = 0; i < old.size(); i++) {
     // TODO: Add in multiple spent_comps, map assemblies to compositions each time we do a depletion...
     old[i]->Transmute(spentFuelComp);
  }
  std::stringstream ss;
  ss << old.size() << " assemblies" << " to discharge burnup " << this->burnup << " MWd/MTHM";
  Record("TRANSMUTE", ss.str());
  std::cerr << ss.str() << std::endl;
}


//TODO - Separate transmute behaviors from ORIGEN behaviors
cyclus::Composition::Ptr reactor::Deplete_(cyclus::Material::Ptr mat, int n_assem, double power) {
    
    // TODO: Store interface as a persistent member but free unneeded components after deplete?
    // TODO: Store burnup somewhere and record discharge burnup each time Transmute is called?
    OrigenInterface::cyclus2origen react;
     
    // Set ORIGEN library path
    react.set_lib_path(lib_path);
    
    // Set ID tags
    if(this->assembly_type == "") {
       std::stringstream ss;
       ss << "Cyborg::reactor::Deplete_() - assembly_type unspecified!" << std::endl; 
       throw cyclus::ValueError(ss.str());
    }    
    react.set_id_tag("Assembly Type",assembly_type);
    
    // Set Interpolable parameters   
    // TODO: Eventually determine enrichment / interpolable dimensions dynamically...
    if(this->enrichment <= 0.0 || this->enrichment > 100.0) {
       std::stringstream ss;
       ss << "Cyborg::reactor::Deplete_() - invalid enrichment specified! this->enrichment = "
          << this->enrichment << std::endl; 
       throw cyclus::ValueError(ss.str());
    }
    react.add_parameter("Enrichment",enrichment);
   
    if(this->mod_density > 0.0) react.add_parameter("Moderator Density",this->mod_density);   
    
    // Set depletion time & power
    std::vector<double> dp_time, dp_pow;
    dp_time.push_back(0.0);
    // Number of cycles is assumed to be proportional to core fraction per batch
    // i.e., 1/3 fraction => 3 cycles
    //for(size_t i=1; i <= round(this->fuel_capacity*1.E3/mat->quantity()); ++i) {
    for(size_t i=0; i < round(this->n_assem_core / this->n_assem_batch); ++i) {
       // Cycle timestep is in months; use years for ORIGEN for simplicity
       dp_time.push_back(static_cast<double>(cycle_time)/12.0 + dp_time.back());        
       // SES TODO: Eventually handle non-uniform cycle powers
       dp_pow.push_back(power);
       //std::cerr << "Pushing back time: " << cycle_length*i*1.0/12.0 << "  power: " << power << std::endl;
       
       // Decay fuel during reload
       dp_time.push_back(dp_time.back() + static_cast<double>(refuel_time)/12.0);
       dp_pow.push_back(0.0);
    }
    react.set_time_units("y");
    react.set_time_steps(dp_time); 
    
    react.set_power_units("watt");
    react.set_powers(dp_pow);  

    // Create cross-section library
    react.interpolate();

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
    

    // Run Calculation
    react.solve();

    // Store burnup (units of MWd/MTU)
    this->burnup = react.burnup_last();  

    // Get materials and convert nuclide ids back to Cyclus format
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
    cyclus::Composition::Ptr comp_out = cyclus::Composition::CreateFromAtom(v);
    return comp_out;
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

extern "C" cyclus::Agent* Constructreactor(cyclus::Context* ctx) {
  return new reactor(ctx);
}

}  // namespace cyborg
