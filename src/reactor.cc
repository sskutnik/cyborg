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
		fuel.Push(fresh_inventory.Pop(fresh_inventory.quantity()));
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
	std::vector<int> in_ids;

	// Get fuel recipe and convert to ORIGEN format
	cyclus::Composition::Ptr comp_in = context()->GetRecipe(fuel_recipe);
	cyclus::CompMap mappy = comp_in->mass();
	for(std::map<int,double>::iterator it = mappy.begin(); it!=mappy.end(); it++){
	int id = it->first;
	double mass_fraction = it->second;
	// convert id 
	in_ids.push_back(pyne::nucname::zzaaam(id));
	}

	//dummy org_id and org_atom until ORIGEN functions implemented
	std::vector<int> org_id;
	std::vector<double> org_atom;

	// After ORIGEN depletion - convert nuclide ids back to Cyclus format
	std::vector<int> out_id;
	for(std::vector<int>::iterator i=org_id.begin(); i!=org_id.end(); i++){
	out_id.push_back(pyne::nucname::id(*i));
	}
	// need to get mass fraction data from ORIGEN
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
