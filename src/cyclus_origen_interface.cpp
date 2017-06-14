#include "cyclus_origen_interface.h"
#include <math.h>
#include <boost/algorithm/string.hpp>
#include "error.h"
#include "Origen/Core/dc/Library.h"
#include "Origen/Core/dc/Material.h"
#include "Origen/Core/dc/TagManager.h"
#include "Origen/Core/dc/ConcentrationConverter.h"
#include "Origen/Core/fn/io.h"
#include "Origen/Core/fn/interp.h"
#include "Origen/Core/io/LibraryIO.h"
#include "Origen/Solver/SolverSelector.h"
#include "ScaleUtils/IO/DB.h"
#include "ScaleUtils/IO/Utils.h"

namespace OrigenInterface {

void cyclus2origen::set_lib_names(const std::vector<std::string> &lib_names){
  b_lib_names.clear();
  b_lib_names.resize(lib_names.size());
  b_lib_names=lib_names;
}

void cyclus2origen::set_lib_path(const std::string lib_path){
  using cyclus::IOError;
  std::string tmp_path = boost::trim_copy(lib_path);

  auto dr = opendir(tmp_path.c_str());
  if(dr==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_lib_path(" << __LINE__ << ") : File '" << tmp_path 
       << "' is not a directory!" << std::endl;
    throw IOError(ss.str());
  }
  closedir(dr);
  b_lib_path=tmp_path;
}

void cyclus2origen::add_lib_names(const std::vector<std::string> &lib_names){
  for(auto lib : lib_names) b_lib_names.push_back(lib);
}

void cyclus2origen::remove_lib_names(const std::vector<std::string> &lib_names){
  for(size_t i = 0; i < lib_names.size(); i++){
    for(size_t j = 0; j < b_lib_names.size(); j++){
      if(b_lib_names[j]==lib_names[i]){
        b_lib_names.erase(b_lib_names.begin()+j);
        break;
      }
    }
  }
}

void cyclus2origen::list_lib_names() const{
  for(auto lib : b_lib_names) std::cout << "Library name: " << lib << '\n';
}

void cyclus2origen::get_lib_names(std::vector<std::string> &lib_names) const{
  using cyclus::StateError;
  if(!lib_names.empty()){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_lib_names(" << __LINE__ 
       << ") : Return vector for lib_names not empty upon function call!" << std::endl;
    throw StateError(ss.str());
  }
  lib_names=b_lib_names;
}

void cyclus2origen::set_id_tag(const std::string idname, const std::string idvalue){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  b_tm->setIdTag(idname,idvalue);
}

void cyclus2origen::set_id_tags(const std::map<std::string,std::string> &tags){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  for(auto tag : tags){
    this->set_id_tag(tag.first,tag.second);
  }
}

void cyclus2origen::remove_id_tag(const std::string idname){
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::remove_id_tag(" << __LINE__ 
       << ") : No tag manager found on this interface object!" << std::endl;
    throw StateError(ss.str());
  }
  if(!b_tm->hasTag(idname)){
    std::stringstream ss;
    ss << "Cyborg::reactor::remove_id_tag(" << __LINE__ 
       << ") : Tag manager does not have a tag with name = " << idname 
       << "!" << std::endl;
    throw StateError(ss.str());
  }
  b_tm->deleteTag(idname);
}

void cyclus2origen::list_id_tags() const{
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::list_id_tags(" << __LINE__ 
       << ") : No tag manager found on this interface object!"  << std::endl;
    throw StateError(ss.str());
  }
  if(b_tm->listIdTags().size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::list_id_tags(" << __LINE__ 
       << ") : No ID tags found on this interface object!  Use set_id_tags()."  
       << std::endl;
    cyclus::Warn<cyclus::WARNING>(ss.str());
    //throw StateError(ss.str());
  }
  for(auto tags : b_tm->listIdTags()){
    std::cout << "Tag name: " << tags << ", value: " << b_tm->getIdTag(tags) 
              << "." << std::endl;
  }
}

void cyclus2origen::get_id_tags(std::vector<std::string> &names, std::vector<std::string> &values) const{
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_id_tags(" << __LINE__ 
       << ") : No tag manager found on this interface object!" << std::endl;
    throw StateError(ss.str());
  }
  if(!names.empty()){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_id_tags(" << __LINE__ 
       << ") : Return vector for ID tag names not emtpy upon function call!" << std::endl;
    throw StateError(ss.str());
  }
  if(!values.empty()){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_id_tags(" << __LINE__ 
       << ") : Return vector for ID tag values not empty upon function call!" << std::endl;
    throw StateError(ss.str());
  }
  for(auto tag : b_tm->listIdTags()){
    names.push_back(tag);
    values.push_back(b_tm->getIdTag(tag));
  }
}

void cyclus2origen::set_materials(const std::vector<int> &ids, const std::vector<double> &concs){
  using cyclus::StateError;
  using cyclus::ValueError;

  if(b_lib_interp==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_materials(" << __LINE__ 
       << ") : No library found on this interface object!  Use interpolate() first." 
       << std::endl;
    throw StateError(ss.str());
  }
  if(b_vol<=0){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_materials(" << __LINE__ 
       << ") : Volume should be implicitly set...." << std::endl;
    throw ValueError(ss.str());
  }
  if(ids.size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_materials(" << __LINE__ 
       << ") : No IDs provided in ID vector!" << std::endl;
    throw StateError(ss.str());
  }
  if(concs.size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_materials(" << __LINE__ 
       << ") : No concentrations were provided in the concentrations vector!" << std::endl;
    throw StateError(ss.str());
  }
  if(concs.size()!=ids.size()){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_materials(" << __LINE__ 
       << ") : Size mismatch between ID and concentrations vectors!" << std::endl;
    throw ValueError(ss.str());
  }
  std::string name = "cyclus_";
  if(std::string::npos==b_interp_name.find(name)){
    name.append(b_interp_name);
  }

  std::vector<int> tmp_ids;
  tmp_ids.assign(ids.begin(),ids.end());

  //  Looping through IDs provided to ensure it is valid and translate to pizzzaaa as necessary.
  for( auto &zaid : tmp_ids ) {     
     if(ScaleData::Utils::is_valid_pizzzaaa(zaid)) continue;    

     if(ScaleData::Utils::is_valid_zzzaaai(zaid)) zaid = ScaleData::Utils::zzzaaai_to_pizzzaaa(zaid);
     if(!ScaleData::Utils::is_valid_pizzzaaa(zaid)){
        std::stringstream ss;
        ss << "Cyborg::reactor::set_materials(" << __LINE__ 
           << ") : Unrecognizeable nuclide ID name! Use zzzaaai or pizzzaaa format." 
           << std::endl;
        throw ValueError(ss.str());
     }
  }

  // Group together nuclide IDs and masses under a concentrations object to handle unit conversions
  // Replace any prior concentrations (if initialized)
  // Units should be set before concentration values are set.
  b_nucset = std::make_shared<Origen::NuclideSet>(tmp_ids);
  b_concs = std::make_shared<Origen::Concentrations>();

  b_concs->set_nuclide_set(*b_nucset);
  b_concs->set_units(b_concUnits);  // Units first. Setting in constructor potentially wonky.
  b_concs->set_vals(concs);         // Then values.


  int id = 1001; // Arbitrary.

  b_mat = Origen::SP_Material(new Origen::Material(b_lib_interp,name,id,b_vol));
  b_mat->set_concs_at(*b_concs,0);
/*
  std::cerr << "Material::initial_mass() = " << b_mat->initial_mass() << " grams." 
            << std::endl;
  std::cerr << "Material::initial_hm_mass() = " << b_mat->initial_hm_mass() 
            << " grams." << std::endl;
*/
  // Takes b_lib and interpolates to new burnups based on b_times and b_powers.
  prob_spec_lib(b_lib_interp,b_times,b_fluxes,b_powers); 

}

void cyclus2origen::reset_material(){
  if(b_mat!=nullptr) b_mat->clear();
}

void cyclus2origen::set_mat_units(const std::string mat_units){
  b_concUnits = Origen::convertStringToConcUnit(mat_units);
}

void cyclus2origen::set_time_units(const char* time_units){
  b_timeUnits = Origen::Time::units(time_units);
}

void cyclus2origen::set_power_units(const char* power_units){
  b_powerUnits = Origen::Power::units(power_units);
}

void cyclus2origen::add_time_step(const double time){
   using cyclus::StateError;
   if(time<0){
      std::stringstream ss;
      ss << "Cyborg::reactor::add_time_step(" << __LINE__ 
         << ") : Time value provided is non-physical (i.e. - less than zero)!" 
         << std::endl;
      throw StateError(ss.str());
   }
   b_times.push_back(time);
}

void cyclus2origen::set_fluxes(const std::vector<double> &fluxes){
  using cyclus::StateError;
  if(fluxes.size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_flux(" << __LINE__ << ") : Vector of fluxes provided is empty!" 
       << std::endl;
    throw StateError(ss.str());
  }
  b_fluxes = fluxes;
}

void cyclus2origen::add_flux(const double flux){
  using cyclus::StateError;
  if(flux<0){
    std::stringstream ss;
    ss << "Cyborg::reactor::add_flux(" << __LINE__ 
       << ") : Flux value provided is non-physical (<0)!" << std::endl;
    throw StateError(ss.str());
  }
  b_fluxes.push_back(flux);
}

void cyclus2origen::delete_fluxes(){
   b_fluxes.clear();
}

void cyclus2origen::set_powers(const std::vector<double> &powers){
  using cyclus::StateError;
  if(powers.size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::set_powers(" << __LINE__ 
       << ") : Vector of powers provided is empty!" << std::endl;
    throw StateError(ss.str());
  }
  b_powers = powers;
}

void cyclus2origen::add_power(const double power){
  using cyclus::StateError;
  if(power<0){
    std::stringstream ss;
    ss << "Cyborg::reactor::add_power(" << __LINE__ 
       << ") : Power provided is non-physical (<0)!" << std::endl;
    throw StateError(ss.str());
  }
  b_powers.push_back(power);
}

void cyclus2origen::delete_powers(){
   b_powers.clear();
}

void cyclus2origen::add_parameter(const std::string name, const double value){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  b_tm->setInterpTag(name,value);
}

void cyclus2origen::set_parameters(const std::map<std::string,double> &params){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  for(auto param : params){
    this->add_parameter(param.first,param.second);
  }
}

void cyclus2origen::remove_parameter(const std::string name){
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::remove_parameter(" << __LINE__ 
       << ") : No tag manager found on this interface object!" << std::endl;
    throw StateError(ss.str());
  }
  if(!b_tm->hasTag(name)){
    std::stringstream ss;
    ss << "Cyborg::reactor::remove_parameter(" << __LINE__ << ") : No tag with name " 
       << name << " found on this interface object!" << std::endl;
    throw StateError(ss.str());
  }

  b_tm->deleteTag(name);
}

void cyclus2origen::list_parameters() const{
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::list_parameters(" << __LINE__ 
       << ") : No tag manager found on this interface object!" << std::endl;
    throw StateError(ss.str());
  }
  if(b_tm->listInterpTags().size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::list_parameters(" << __LINE__ 
       << ") : No parameters found on this tag manager!" << std::endl;
    cyclus::Warn<cyclus::WARNING>(ss.str());
    //throw StateError(ss.str());
  }
  for(auto tag : b_tm->listInterpTags()){
    std::cout << "Interp tag name: " << tag << ", value: " << b_tm->getInterpTag(tag) 
              << "." << std::endl;
  }
}

void cyclus2origen::get_parameters(std::vector<std::string> &names, std::vector<double> &values) const{
  using cyclus::StateError;
  if(!names.empty()){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_parameters(" << __LINE__ 
       << ") : Return vector for names not empty upon function call!" << std::endl;
    throw StateError(ss.str());
  }
  if(!values.empty()){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_parameters(" << __LINE__ 
       << ") : Return vector for values not empty upon function call!" << std::endl;
    throw StateError(ss.str());
  }
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::get_parameters(" << __LINE__ 
       << ") : No tag manager present with parameters!" << std::endl;
    throw StateError(ss.str());
  }
  for(auto tag : b_tm->listInterpTags()){
    names.push_back(tag);
    values.push_back(b_tm->getInterpTag(tag));
  }
}

void cyclus2origen::interpolate() {
  using cyclus::ValueError;
  using cyclus::StateError;
  if(b_tm==NULL){
    std::stringstream ss;
    ss << "Cyborg::reactor::interpolate(" << __LINE__ << ") : No tag manager found!" 
       << std::endl;
    throw StateError(ss.str());
  }
  if(b_tm->listIdTags().size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::interpolate(" << __LINE__ << ") : No ID tags found!"
       << std::endl;
    throw StateError(ss.str());
  }

  if(b_tagman_list.size() > 0)
  {
    for(auto& tm : b_tagman_list)
    {
      if(!tm->interpolationCompare(*b_tm))
      {
        b_tagman_list.clear();
        b_lib_names.clear();
        break;
      }
    }
  }
/*
  if(b_lib_names.size()==0 && b_lib_path.size()==0 && b_tagman_list.size()==0){
    std::stringstream ss;
    ss << "Cyborg::reactor::interpolate(" << __LINE__ << ") : No library names or path specified!" << std::endl;
    throw cyclus::ValueError(ss.str());
  }
*/
  std::vector<Origen::TagManager> libTMs;
  if(b_lib_names.size()==0 && b_tagman_list.size()==0){

    if(b_lib_path.size()==0) {
      std::stringstream ss;
      ss << "Cyborg::reactor::interpolate(" << __LINE__  
         << ") : No library names or path specified!" << std::endl;
      throw cyclus::ValueError(ss.str());
    }

    struct dirent *drnt;
    auto dr = opendir(b_lib_path.c_str());
    std::string midstring = "";
    if(&(b_lib_path.back()) != "/") midstring = "/";
    while(true){
      drnt=readdir(dr);
      if(!drnt) break;
      std::string lib_name (drnt->d_name);
      if(lib_name=="." || lib_name=="..") continue;
      lib_name = b_lib_path + midstring + lib_name;
      struct stat buffer;
      if(stat(lib_name.c_str(), &buffer) != 0){
        std::cout << lib_name << " doesn't exist!" << std::endl;
      }else{
        b_lib_names.push_back(lib_name);
      }
    }
    closedir(dr);

    // Bail if no libraries specified
    if(b_lib_names.size() == 0) {
      std::stringstream ss;
      ss << "Cyborg::reactor::interpolate(" << __LINE__ 
         << ") : No libraries specified or found!" << std::endl;
      throw ValueError(ss.str());
    }

    //std::vector<Origen::SP_TagManager> tms = Origen::collectLibrariesParallel(b_lib_names);
    // Serial for now until I can get the repo update working
    std::vector<Origen::SP_TagManager> tms = Origen::collectLibraries(b_lib_names); 
    for(auto& tm : tms) libTMs.push_back(*tm);

    // Bail if no libraries found
    if(libTMs.size() == 0) {
      std::stringstream ss;
      ss << "Cyborg::reactor::interpolate(" << __LINE__ 
         << ") : No libraries found that have tag managers!" << std::endl;
      throw ValueError(ss.str());
    }

    // Down-select to libraries matching specified ID tags
    libTMs = Origen::selectLibraries(libTMs,*b_tm);
    b_tagman_list.clear();
    for(auto tm : libTMs) b_tagman_list.push_back(std::make_shared<Origen::TagManager>(tm));
    if(b_tagman_list.size() == 0){
      std::stringstream ss;
      ss << "Cyborg::reactor::interpolate(" << __LINE__ 
         << ") : No libraries found that match specified ID tags!" << std::endl;
      throw ValueError(ss.str());
    }
  }
  else
  {
    // We already have a list of identified libraries to interpolate
    for(auto& tm : b_tagman_list) libTMs.emplace_back(*tm);
  }
  b_lib_interp = Origen::interpLibraryND(libTMs,*b_tm);
  b_interp_name = (b_lib_interp->scp_tag_manager())->getIdTag("Filename");
}

void cyclus2origen::solve() {
   using cyclus::StateError;
   using cyclus::ValueError;
   
   if(b_powers.size()==0 && b_fluxes.size()==0) {
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__ 
         << ") : No powers or fluxes found on this interface object!  Use set_fluxes or set_powers."
         << std::endl;
      throw StateError(ss.str());
   }

   this->solve(b_times, b_fluxes, b_powers);
}

void cyclus2origen::solve(std::vector<double>& times, std::vector<double>& fluxes, std::vector<double>& powers){
   using cyclus::StateError;
   using cyclus::ValueError;
   if(b_mat==NULL){
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__
         << ") : No material object found on this interface object!"
         << "  Use set_materials() or set_materials_with_masses()." << std::endl;
      throw StateError(ss.str());
   }
   if(b_mat->library()==NULL) {
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__ 
         << ") : No library object found on the material object on this interface object!"
         << std::endl;
      throw StateError(ss.str());
   }
   if(powers.size()==0 && fluxes.size()==0){
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__ 
         << ") : No powers or fluxes specified for depletion!" << std::endl;
      throw StateError(ss.str());
   }
   if(powers.size()!=(times.size()-1) && powers.size()!=0){
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__ 
         << ") : Powers vector must be exactly 1 element shorter than times vector!" 
         << std::endl << "Power vector size is " << powers.size() 
         << " and times vector size is " << times.size() << "." << std::endl;
      throw ValueError(ss.str());
   }

   if(b_mat->amount_at(0)->size() == 0){
      std::stringstream ss;
      ss << "Cyborg::reactor::solve(" << __LINE__ 
         << ") : Materials object has no material masses set."
         << "  Run set_materials first." << std::endl;
      throw StateError(ss.str());
   }

   // Initialize the solver
   Origen::SP_Solver solver;
   ScaleUtils::IO::DB db;
   db.set<std::string>("solver","cram");
   solver = Origen::SolverSelector::get_solver(db);
   b_mat->set_solver(solver);
    
   size_t num_steps = powers.size()>fluxes.size() ? powers.size() : fluxes.size();
   size_t libPos = 0;
   std::vector<double> dt_rel(4, 0.25);

   auto powFluxIter = (powers.size() > 0) ? powers.begin() : fluxes.begin();
   for(size_t i = 0; i < num_steps; i++){             
     b_mat->add_step((times[i+1]-times[i])*Origen::Time::factor<Origen::Time::SECONDS>(this->b_timeUnits));
     b_mat->set_transition_matrix(b_lib->newsp_transition_matrix_at(libPos));
      
     if(b_fluxes.size()==0) {
        b_mat->set_power(*powFluxIter);
        std::cerr << "(*powFluxIter) = " << (*powFluxIter) << " - libPos = " << libPos <<  std::endl;
     } else {
        b_mat->set_flux(*powFluxIter);
     }
     std::vector<double> tmpFlux, tmpPower;
     b_mat->solve(dt_rel, &tmpFlux, &tmpPower);

     // Use the next library position if the next time has a non-zero power
     if(*std::next(powFluxIter,1) > 0) ++libPos;
     std::advance(powFluxIter,1);
   }
}

void cyclus2origen::get_masses(std::vector<std::vector<double> > &masses_out, const std::string units) const{

  masses_out.clear();
  masses_out.resize(b_times.size());

  for(size_t i=0; i < b_times.size(); ++i) {
     this->get_masses_at(i,masses_out[i], units);
  }
}

void cyclus2origen::get_masses_at(int p, std::vector<double> &masses_out, const std::string units) const{
   using cyclus::ValueError;

   if( p < 0 || p >= b_mat->ntimes() ){
     std::stringstream ss;
     ss << "Cyborg::reactor::get_masses_at(" << __LINE__ << ") : Step requested " 
        << p << " falls outside the bounds [0," << b_mat->ntimes() << ")!" << std::endl;
     throw ValueError(ss.str());
   }

   Origen::ConcentrationUnit concUnits = Origen::convertStringToConcUnit(units);

   Origen::SP_NuclideSet tmpNucSet = std::make_shared<Origen::NuclideSet>(*(b_mat->sizzzaaa_list()));
   Origen::SP_Concentrations tmpConcs = std::make_shared<Origen::Concentrations>(); // Get values off of materials object, then set_units to desired units.
   tmpConcs->set_nuclide_set(*tmpNucSet);

   masses_out.clear();
   b_mat->get_concs_at(tmpConcs.get(), p);                    // Get values.
//   tmpConcs->set_units(Origen::ConcentrationUnit::KILOGRAMS); // Set units.
   tmpConcs->set_units(concUnits); // Set units.
 
   tmpConcs->get_vals(masses_out); 
}

void cyclus2origen::get_masses_at_map(int p, std::map<int,double> &masses_out, const std::string id_type, const std::string units) const{
   using cyclus::ValueError;

   if( p < 0 || p >= b_mat->ntimes() ){
       std::stringstream ss;
       ss << "Cyborg::reactor::get_masses_at_map(" << __LINE__ 
          << "): Step requested " << p << " falls outside the bounds [0," 
          << b_mat->ntimes() << ")!" << std::endl;
       throw ValueError(ss.str());
   }

   std::vector<int> ids;   
   if(id_type == "sizzzaaa"){
      this->get_ids(ids);
   }else if(id_type == "zzzaaai"){
       this->get_ids_zzzaaai(ids);
   }else{
       std::stringstream ss;
       ss << "Cyborg::reactor::get_masses_at_map(" << __LINE__ 
          << ") : Type of nuclide ids requested is not recognized."
          << "  Must be 'sizzzaaa' or 'zzzaaai'." << std::endl;
       throw ValueError(ss.str());
   }

   std::vector<double> concs;
   this->get_masses_at(p, concs, units);

   if(ids.size() != concs.size())
   {
       std::stringstream ss;
       ss << "Cyborg::reactor::get_masses_at_map(" << __LINE__ 
          << ") : Number of nuclide ids (" << ids.size() 
          << ") does not match the number of concentrations (" << concs.size()  
          << ")." << std::endl;
       throw ValueError(ss.str());
   }

   for(size_t i = 0; i < ids.size(); i++) masses_out[ids[i]] = concs[i];
}

void cyclus2origen::get_masses_final(std::vector<double> &masses_out, const std::string units) const{
  this->get_masses_at(b_mat->nsteps(), masses_out, units);
}

void cyclus2origen::get_masses_final_map(std::map<int,double> &masses_out, const std::string id_type, const std::string units) const
{
   using cyclus::ValueError;

   std::vector<int> ids;   
   if(id_type == "sizzzaaa"){
      this->get_ids(ids);
   }else if(id_type == "zzzaaai"){
       this->get_ids_zzzaaai(ids);
   }else{
       std::stringstream ss;
       ss << "Cyborg::reactor::get_masses_final_map(" << __LINE__ 
          << ") : Type of nuclide ids requested is not recognized."
          << "  Must be 'sizzzaaa' or 'zzzaaai'." << std::endl;
       throw ValueError(ss.str());
   }

   std::vector<double> concs;
   this->get_masses_final(concs, units);

   if(ids.size() != concs.size())
   {
       std::stringstream ss;
       ss << "Cyborg::reactor::get_masses_final_map(" << __LINE__ 
          << ") : Number of nuclide ids (" << ids.size() 
          << ") does not match the number of concentrations (" << concs.size() 
          << ")." << std::endl;
       throw ValueError(ss.str());
   }

   for(size_t i = 0; i < ids.size(); i++) masses_out[ids[i]] = concs[i];
}

void cyclus2origen::get_ids(std::vector<int> &ids_out) const{
  ids_out = *(b_mat->sizzzaaa_list());
}

void cyclus2origen::get_ids_zzzaaai(std::vector<int> &ids_out) const{
  for(auto id : *(b_mat->sizzzaaa_list())){
    ids_out.push_back(ScaleData::Utils::pizzzaaa_to_zzzaaai(id));
  }
}

double cyclus2origen::burnup_last() const {    
   return this->burnup_at(b_mat->nsteps());
}

double cyclus2origen::burnup_at(const int stepNum) const {
   using cyclus::ValueError;
   using cyclus::StateError;

   if(b_mat->nsteps() == 0) {
     std::stringstream ss;
     ss << "Cyborg::reactor::burnup_at(" << __LINE__  
        << ") : No burnup steps found!" << std::endl;
     throw StateError(ss.str());
     return -1.0;
   }
   if(stepNum < 0 || stepNum > b_mat->nsteps()) {
      std::stringstream ss;
      ss << "Cyborg::reactor::burnup_at(" << __LINE__ << ") : Step requested " 
         << stepNum << " falls outside the bounds [0," << b_mat->nsteps() 
         << ")!" << std::endl;
      throw ValueError(ss.str());
      return -1.0;
   }   
   return b_mat->burnup_at(stepNum ); 
}

std::vector<double> cyclus2origen::get_burnups() const { 
   
   using cyclus::StateError;

   std::vector<double> burnups;
   if(b_mat->nsteps() == 0) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_burnups(" << __LINE__ 
        << ") : No burnup steps found!" << std::endl;
     throw StateError(ss.str());
     return burnups;
   }
   for(size_t i=0; i < b_mat->nsteps(); ++i) {
      burnups.push_back(this->burnup_at(i));
   }
}

std::vector<double> cyclus2origen::get_times(std::string units) const { 
   
   using cyclus::StateError;

   std::vector<double> times;
   if(this->b_times.size() == 0) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_times(" << __LINE__ 
        << ") : No times found!" << std::endl;
     throw StateError(ss.str());
     return times;
   }
  
   Origen::Time::UNITS tmpUnits = Origen::Time::units(units.c_str());
   if(tmpUnits == Origen::Time::UNITS::UNKNOWN) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_times(" << __LINE__ << ") : Unknown time units: "
        << units << "; unable to convert times!" << std::endl;
     throw StateError(ss.str());
     return times;
   }
   // Convert power units if requested
   if( tmpUnits != b_timeUnits && tmpUnits != Origen::Time::UNITS::UNKNOWN) {
     for(size_t i=0; i < b_times.size(); ++i) {
       times.push_back( b_times.at(i) / Origen::Time::factor(tmpUnits, b_timeUnits));
      }
   }
   else {
     times = b_times;
   }
   return times;
}

std::vector<double> cyclus2origen::get_powers(std::string units) const { 
   
   using cyclus::StateError;

   std::vector<double> powers;
   if(this->b_powers.size() == 0) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_powers(" << __LINE__ << ") : No powers found!" << std::endl;
     throw StateError(ss.str());
     return powers;
   }
  
   Origen::Power::UNITS tmpUnits = Origen::Power::units(units.c_str());
   if(tmpUnits == Origen::Power::UNITS::UNKNOWN) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_powers(" << __LINE__ << ") : Unknown power units: " 
        << units << "; unable to convert powers!" << std::endl;
     throw StateError(ss.str());
     return powers;
   }
   // Convert power units if requested
   if( tmpUnits != b_powerUnits) {
     for(size_t i=0; i < b_powers.size(); ++i) {
       powers.push_back( b_powers.at(i) / Origen::Power::factor(tmpUnits, b_powerUnits));
      }
   }
   else {
     powers = b_powers;
   }
   return powers;
}

void cyclus2origen::prob_spec_lib(Origen::SP_Library lib, const std::vector<double> &times, 
                                  const std::vector<double> &fluxes, const std::vector<double> &powers) {
   std::vector<double> timeTmp = times;
   std::vector<double> powTmp; 

   if(b_timeUnits!=Origen::Time::DAYS){
      //for(auto &time : timeTmp) time *= Origen::Time::factor(Origen::Time::DAYS, b_timeUnits);
      double timeFactor = Origen::Time::factor(Origen::Time::DAYS, b_timeUnits);
      std::transform(timeTmp.begin(), timeTmp.end(), timeTmp.begin(),
                     std::bind1st(std::multiplies<double>(),timeFactor));
   }
   if(!fluxes.empty() && powers.empty() ) {
      for(auto& flux : fluxes) powTmp.push_back(flux*b_mat->power_factor_bos());
   } else if(!(fluxes.empty() || powers.empty()) ) {
      std::stringstream ss;
      ss << "Cyborg::reactor::prob_spec_lib(" << __LINE__ 
         << ") : Both the fluxes and powers vectors have values! Choose one!" << std::endl;
      throw cyclus::ValueError(ss.str());
   }
   else {
      powTmp = powers;
   }

   if(times.size()!=powers.size()+1){
      std::stringstream ss;
      ss << "Cyborg::reactor::prob_spec_lib(" << __LINE__ 
         << ") : Powers or fluxes vectors not exactly 1 element shorter than times vector!" 
         << std::endl << "Powers or fluxes vector has " << powers.size() 
         << " elements and times vector has " << times.size() << " elements." << std::endl;
      throw cyclus::ValueError(ss.str());
   }

   std::vector<double> cycBU, interpBU;
   double deltaBU;
   cycBU.push_back(0.0);

   size_t nonZeroPowers = std::count_if(powers.begin(), powers.end(), [](double p){return p > 0;});

   for(size_t i = 0; i < powers.size(); i++){
      if(! (powers[i] > 0)) continue;

      // Powers in watts, times in days, and hm mass in grams => buTmp in MWd/MTU (equiv. to W*d/g)
      deltaBU = powTmp[i]*(timeTmp[i+1]-timeTmp[i])/b_mat->initial_hm_mass();
      // NOTE: We're interpolating to cycle MIDPOINT burnup, not end burnup
      interpBU.push_back(deltaBU/2 + cycBU.back());
      cycBU.push_back(deltaBU + cycBU.back());
   }

   if(interpBU.size() != nonZeroPowers){
      std::stringstream ss;
      ss << "Cyborg::reactor::prob_spec_lib(" << __LINE__ 
         << ") : Calculated burnup vector does not have same size as provided powers vector!" 
         << std::endl;
      throw cyclus::StateError(ss.str());
   }
   b_lib = lib->interpolate_Interp1D(interpBU);
}

const std::string cyclus2origen::get_tag_manager_string() const
{
  if(!b_tm) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_tag_manager_string(" << __LINE__ 
        << ") : TagManager not initialized!" << std::endl;
     throw cyclus::StateError(ss.str());
  }
/*
 * SES: Shouldn't need an interpolated library to hash the state; everything needed available pre-interpolation
  if(!b_lib) {
     std::stringstream ss;
     ss << "Cyborg::reactor::get_tag_manager_string(" << __LINE__ << ") : ORIGEN library not yet initialized!" << std::endl;
     throw cyclus::StateError(ss.str());
  }
  if(b_mat->nsteps()==0)
  {
    std::stringstream ss;
    ss << "Cyborg::reactor::get_tag_manager_string(" << __LINE__ << ") : Depletion calculation has not yet occurred." << std::endl;
    throw cyclus::StateError(ss.str());
  }
  //b_lib->get_tag_manager(tm);
*/
  if(b_times.size() == 0 || b_powers.size() == 0) {
     std::stringstream ss;
     ss << "WARNING: Cyborg::reactor::get_tag_manager_string(" << __LINE__ 
        << ") : Power history not initialized!" << std::endl;
     cyclus::Warn<cyclus::WARNING>(ss.str());
  }

  Origen::TagManager tm(*b_tm);

  std::stringstream times;
  for(auto time : b_times)
  {
    times << time << ",";
  }

  std::string times_val = times.str();
  times_val.pop_back();
  std::stringstream times_name;
  times_name << "Times (" << Origen::Time::name(b_timeUnits) << ")";
  tm.setIdTag(times_name.str(),times_val);

  std::stringstream powers;
  for(auto power : b_powers)
  {
    powers << power << ",";
  }

  std::string powers_val = powers.str();
  powers_val.pop_back();
  std::stringstream powers_name;
  powers_name << "Powers (" << Origen::Power::name(b_powerUnits) << ")";
  tm.setIdTag(powers_name.str(),powers_val);

  return tm.to_string();
}


}//namespace
