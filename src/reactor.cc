#include "reactor.h"

namespace cyborg {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
reactor::reactor(cyclus::Context* ctx) : cyclus::Facility(ctx), reactor_time(1), decom(false) {}

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
			Discharge_(3);
		}
		else if (reactor_time == reactor_lifetime){
			Discharge_(1);
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

void reactor::Discharge_(int core_fraction) {
	// Pop 1/3 of core from fuel buffer (except on retiring time step)
	cyclus::Material::Ptr to_burn = fuel.Pop(fuel.capacity()/core_fraction);

	// Transmute material ready for discharge
	to_burn = Deplete_(to_burn);

	// Discharge fuel to spent fuel buffer
	spent_inventory.Push(to_burn);
}

cyclus::Material::Ptr reactor::Deplete_(cyclus::Material::Ptr mat) {
	cyclus::CompMap v;
	v[922350000] = 1;
	v[922380000] = 10;
	cyclus::Composition::Ptr comp = cyclus::Composition::CreateFromAtom(v);

	mat->Transmute(comp);
	return mat;
}

extern "C" cyclus::Agent* Constructreactor(cyclus::Context* ctx) {
  return new reactor(ctx);
}

}  // namespace cyborg
