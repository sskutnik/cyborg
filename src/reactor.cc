#include "reactor.h"
#include "cyclus_origen_interface.h"
#include "error.h"

namespace cyborg {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reactor::reactor(cyclus::Context* ctx) : cyclus::Facility(ctx), reactor_time(1), decom(false) {
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string reactor::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::EnterNotify(){
    Facility::EnterNotify();
    buy_policy.Init(this, &fresh_inventory, std::string("fresh_inventory"));
    buy_policy.Set(fresh_fuel).Start();

    sell_policy.Init(this, &spent_inventory, std::string("spent_inventory")).Set(spent_fuel).Start();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tick() {
    if (!decom) {
        fuel.capacity(fuel_capacity*1000);
        fresh_inventory.capacity(fuel.space());
    }
    else if (decom) {
        fresh_inventory.capacity(0);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void reactor::Tock() {
    // Only continue to operate if not exceeding reactor lifetime
    if (!decom) {
        // Load Fuel
        Load_();
        // Transmute & Discharge if necessary
        if (reactor_time % cycle_length == 0 && reactor_time != reactor_lifetime) {
            Discharge_(3.0);
        }
        else if (reactor_time == reactor_lifetime){
            Discharge_(1.0);
            decom = true;
        }
        ++reactor_time;
    }
    std::cerr << "Finished Tock()" << std::endl;
}

void reactor::Load_() {
    if (fuel.space() > 0){
        // Push material to fuel buffer from fresh inventory
        //fuel.Push(fresh_inventory.Pop(fresh_inventory.quantity()));
        double toLoad = std::min(fuel.space(),fresh_inventory.quantity());
        fuel.Push(fresh_inventory.Pop( toLoad ));
    }
}

void reactor::Discharge_(double core_fraction) {
    // Pop 1/3 of core from fuel buffer (except on retiring time step)
    cyclus::Material::Ptr to_burn = fuel.Pop(fuel.capacity()/core_fraction);

    std::cerr << "Calling deplete..." << std::endl;
    // Transmute material ready for discharge
    // Assume even power between cycles (for now)
    double cyclePower = power_cap*core_fraction;
    to_burn = Deplete_(to_burn, cyclePower);
    std::cerr << "Finished deplete." << std::endl;

    // Discharge fuel to spent fuel buffer
    spent_inventory.Push(to_burn);
}

cyclus::Material::Ptr reactor::Deplete_(cyclus::Material::Ptr mat, double power) {
    
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
    if(this->enrichment <= 0.0 || this->enrichment > 100.0) {
       std::stringstream ss;
       ss << "Cyborg::reactor::Deplete_() - invalid enrichment specified! this->enrichment = "
          << this->enrichment << std::endl; 
       throw cyclus::ValueError(ss.str());
    }
    std::cerr << "Enrichent = " << this->enrichment << std::endl; 
    react.add_parameter("Enrichment",enrichment);
   
    std::cerr << "Moderator Density = " << this->mod_density << std::endl;
    if(this->mod_density > 0.0) react.add_parameter("Moderator Density",this->mod_density);   

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

    // Convert each isotopic mass to absolute isotopic mass in grams (from unnormalized mass fraction)
    double massNorm = std::accumulate(mass_fraction.begin(),mass_fraction.end(), 0.0);    
    std::transform(mass_fraction.begin(), mass_fraction.end(), norm_mass.begin(), 
                   std::bind1st(std::multiplies<double>(),mat->quantity()/massNorm*1000.0));

    for(auto mass : norm_mass) { std::cerr << "Normed mass: " << mass << std::endl; }
    react.set_materials(in_ids,norm_mass);
    
    // Set depletion time
    // SES: ASSUMING TIME IN MONTHS; NEED TO VALIDATE THIS
    // SES: Also, need to make this a loop based on # of timesteps?
    std::vector<double> dp_time;
    dp_time.push_back(0.0);
    dp_time.push_back(cycle_length);
    
    //for(auto dt : dp_time) { std::cerr << "Pushing back time: " << dt << std::endl; }

    react.set_time_steps(dp_time); 
    react.set_time_units("years");
    
    // Set reactor power 
    std::vector<double> dp_pow;
    dp_pow.push_back(power);
    react.set_powers(dp_pow);  
 
    // Run Calculation
    react.solve();

    // Get materials and convert nuclide ids back to Cyclus format
    std::vector<int> org_id;
    react.get_ids_zzzaaai(org_id);
     
    std::for_each(org_id.begin(), org_id.end(), [](int &nucID){ pyne::nucname::id(nucID); });
   
    // Get mass data from ORIGEN
    std::vector<double> org_atom;
    //react.get_concentrations_final(org_atom);
    react.get_masses_final(org_atom,"atoms ppm");

    cyclus::CompMap v;
    for(int j=0; j!=org_id.size(); ++j){
       // ***TEMPORARY*** suppression of small nuclide inventories to allow for testing
       if(org_atom[j] > 1E-12) { v[org_id[j]] = org_atom[j]*1.E6;
          std::cerr << "Setting v[" << org_id[j] << "] to: " << org_atom[j]*1.E6 << std::endl;
       }
    }
    cyclus::Composition::Ptr comp_out = cyclus::Composition::CreateFromAtom(v);
    mat->Transmute(comp_out);

    return mat;
}

extern "C" cyclus::Agent* Constructreactor(cyclus::Context* ctx) {
  return new reactor(ctx);
}

}  // namespace cyborg
