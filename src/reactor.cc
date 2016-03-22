#include "reactor.h"

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

	// Transmute material ready for discharge
	to_burn = Deplete_(to_burn);

	// Discharge fuel to spent fuel buffer
	spent_inventory.Push(to_burn);
}

cyclus::Material::Ptr reactor::Deplete_(cyclus::Material::Ptr mat) {
	cyclus2origen react;
    
    // Set ORIGEN library path
    react.set_lib_path(lib_path);
    
    // Set ID tags 
    react.set_id_tag("Assembly Type",assembly_type);
    
    // Set Interpolable parameters
    react.add_parameter("Enrichment",enrichment);
    
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
            // convert id 
        in_ids.push_back(pyne::nucname::zzaaam(id));
	}
        // normalize mass fractions
    std::vector<double>::iterator big = std::max_element(mass_fraction.begin(),mass_fraction.end());
    std::vector<double>::iterator little = std::min_element(mass_fraction.begin(),mass_fraction.end());
    std::vector<double> norm_mass;
    double mat_mass = mat->quantity();
    double norm;
    for(std::vector<double>::iterator it = mass_fraction.begin(); it!=mass_fraction.end(); it++){
        norm = (*it - *little)/(*big - *little) * (mat_mass);
        norm_mass.push_back(norm);
    }
    react.set_materials_with_masses(in_ids,norm_mass);
    
    // Set depletion time
    std::vector<double> dp_time;
    dp_time.push_back(cycle_length);
    react.set_time_steps(dp_time);
    react.set_time_units("months");
    
    // Set reactor power 
    std::vector<double> dp_pow;
    dp_pow.push_back(power_cap*1000000);
    react.set_powers(dp_pow);
    
    // Run Calculation
    react.solve();

	// Get materials and convert nuclide ids back to Cyclus format
	std::vector<int> org_id, out_id;
    react.get_ids_zzzaaai(org_id);
	for(std::vector<int>::iterator i=org_id.begin(); i!=org_id.end(); i++){
	out_id.push_back(pyne::nucname::id(*i));
	}
	// Get mass data from ORIGEN
    std::vector<double> org_atom;
    react.get_concentrations_final(org_atom);
	cyclus::CompMap v;
	for(int j=0; j!=out_id.size(); ++j){
		v[out_id[j]] = org_atom[j];
	}
    
	cyclus::Composition::Ptr comp_out = cyclus::Composition::CreateFromAtom(v);
	mat->Transmute(comp_out);

	return mat;
}

extern "C" cyclus::Agent* Constructreactor(cyclus::Context* ctx) {
  return new reactor(ctx);
}

}  // namespace cyborg
