#include "cyclus_origen_interface.h"
#include "Origen/Core/dc/ConcentrationConverter.h"
#include "Origen/Core/fn/io.h"
#include "Origen/Core/fn/interp.h"
#include "Origen/Core/io/LibraryIO.h"
#include "ScaleUtils/IO/DB.h"
#include "ScaleUtils/IO/Utils.h"
#include "Nemesis/harness/DBC.hh"

namespace OrigenInterface {

void cyclus2origen::set_lib_names(const std::vector<std::string> &lib_names){
  b_lib_names.clear();
  b_lib_names.resize(lib_names.size());
  b_lib_names=lib_names;
}

void cyclus2origen::set_lib_path(const std::string lib_path){
  Check(ScaleUtils::IO::directoryExists(lib_path));
  b_lib_path=lib_path;
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
  for(auto lib : b_lib_names) std::cout << "Library name: " << lib << std::endl;
}

void cyclus2origen::get_lib_names(std::vector<std::string> &lib_names) const{
  Check(lib_names.empty());
  lib_names=b_lib_names;
}

void cyclus2origen::set_id_tag(const std::string idname, const std::string idvalue){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  b_tm->setIdTag(idname,idvalue);
}

void cyclus2origen::set_id_tags(const std::map<std::string,std::string> &tags){
  for(auto tag : tags){
    this->set_id_tag(tag.first,tag.second);
  }
}

void cyclus2origen::remove_id_tag(const std::string idname){
  Check(b_tm!=NULL);
  Check(b_tm->hasTag(idname));
  b_tm->deleteTag(idname);
}

void cyclus2origen::list_id_tags() const{
  Check(b_tm!=NULL);
  Check((b_tm->listIdTags()).size()!=0);
  for(auto tags : b_tm->listIdTags()){
    std::cout << "Tag name: " << tags << ", value: " << b_tm->getIdTag(tags) << "." << std::endl;
  }
}

void cyclus2origen::get_id_tags(std::vector<std::string> &names, std::vector<std::string> &values) const{
  Check(names.empty());
  Check(values.empty());
  Check(b_tm!=NULL);
  for(auto tag : b_tm->listIdTags()){
    names.push_back(tag);
    values.push_back(b_tm->getIdTag(tag));
  }
}

void cyclus2origen::set_materials_with_masses(std::vector<int> &ids, const std::vector<double> &masses){
  Check(b_lib!=NULL);
  Check(ids.size()>0);
  Check(ids.size()==masses.size());
  std::vector<double> concs;
  Origen::ConcentrationConverter cv;
  for(size_t i = 0; i < ids.size(); i++){
    concs.push_back(cv.convert_to(Origen::ConcentrationUnit::CM_2_BARN,ids[i],Origen::ConcentrationUnit::KILOGRAMS,masses[i],b_vol));
  }
  this->set_materials(ids,concs);
}

void cyclus2origen::set_materials(std::vector<int> &ids, const std::vector<double> &concs){
  Check(b_lib!=NULL);
  Check(b_vol>0);
  Check(ids.size()>0);
  Check(concs.size()==ids.size());
  std::string name = "cyclus_";
  name.append(b_interp_name);
  int id = 1001;
  for(size_t i = 0; i < ids.size(); i++){
    if(ScaleData::Utils::is_valid_zzzaaai(ids[i])) ids[i] = ScaleData::Utils::zzzaaai_to_pizzzaaa(ids[i]);
    Check(ScaleData::Utils::is_valid_pizzzaaa(ids[i]));
  }
  Origen::SP_Material mat(new Origen::Material(b_lib,name,id,b_vol));
  mat->set_numden_bos(concs,ids,b_vol);
  b_mat = mat;
}

void cyclus2origen::set_mat_units(const std::string mat_units){
  b_concs_units = Origen::convertStringToConcUnit(mat_units);
}

void cyclus2origen::set_time_units(const char* time_units){
  b_time_units = Origen::Time::units(time_units);
}

void cyclus2origen::set_flux(const std::vector<double> &fluxes){
  Check(fluxes.size()!=0);
  b_fluxes = fluxes;
}

void cyclus2origen::set_powers(const std::vector<double> &powers){
  Check(powers.size()!=0);
  b_powers = powers;
}

void cyclus2origen::add_parameter(const std::string name, const double value){
  if(b_tm==NULL) b_tm = Origen::SP_TagManager(new Origen::TagManager());
  b_tm->setInterpTag(name,value);
}

void cyclus2origen::set_parameters(const std::map<std::string,double> &params){
  for(auto param : params){
    this->add_parameter(param.first,param.second);
  }
}

void cyclus2origen::remove_parameter(const std::string name){
  Check(b_tm!=NULL);
  Check(b_tm->hasTag(name));
  b_tm->deleteTag(name);
}

void cyclus2origen::list_parameters() const{
  Check(b_tm!=NULL);
  Check((b_tm->listInterpTags()).size()!=0);
  for(auto tag : b_tm->listInterpTags()){
    std::cout << "Interp tag name: " << tag << ", value: " << b_tm->getInterpTag(tag) << "." << std::endl;
  }
}

void cyclus2origen::get_parameters(std::vector<std::string> &names, std::vector<double> &values) const{
  Check(names.empty());
  Check(values.empty());
  Check(b_tm!=NULL);
  for(auto tag : b_tm->listInterpTags()){
    names.push_back(tag);
    values.push_back(b_tm->getInterpTag(tag));
  }
}

//void set_solver(const int svlr=0){
//  Defaulting to CRAM.  Should be this easy to swap in
//  matrex if desired, but not including it for now.
//
//  if(slvr==1){
//    b_slv=Origen::Solver_matrex();
//  } else {
//    b_slv=Origen::Solver_cram();
//  }
//}

void cyclus2origen::interpolate(){
  Check(b_tm!=NULL);
  Check(b_tm->listIdTags().size()!=0);
  Check(b_lib_names.size()!=0 || b_lib_path.size()!=0);
  if(b_lib_names.size()==0){
    // figure out how to grab filenames from a directory.
    struct dirent *drnt;
    auto dr = opendir(b_lib_path.c_str());
    std::string midstring = "";
    if(&(b_lib_path.back()) != "/") midstring = "/";
    while(true){
      drnt=readdir(dr);
      if(!drnt) break;
      std::string lib_name (drnt->d_name);
      if(lib_name=="." || lib_name==".." || lib_name=="arp2obl.py" || lib_name=="arpdata.txt") continue;
      lib_name = b_lib_path + midstring + lib_name;
      struct stat buffer;
      if(stat(lib_name.c_str(), &buffer) != 0){
        std::cout << lib_name << " doesn't exist!" << std::endl;
      }else{
        b_lib_names.push_back(lib_name);
      }
    }
    closedir(dr);
  }

  std::vector<Origen::SP_TagManager> tms = Origen::collectLibraries(b_lib_names);
  std::vector<Origen::TagManager> tagman;
  for(auto& tm : tms) tagman.push_back(*tm);
  tagman = Origen::selectLibraries(tagman,*b_tm);
  b_lib = Origen::interpLibraryND(tagman,*b_tm);
  b_interp_name = (b_lib->scp_tag_manager())->getIdTag("Filename");
}

void cyclus2origen::solve(){
  Check(b_mat!=NULL);
  Check(b_powers.size()!=0||b_fluxes.size()!=0);
  Check(b_powers.size()==(b_times.size()-1)||b_fluxes.size()==(b_times.size()-1));
  Origen::SP_Solver solver;
  ScaleUtils::IO::DB db;
  db.set<std::string>("solver","cram");
  solver = Origen::SolverSelector::get_solver(db);
  b_mat->set_solver(solver);
  prob_spec_lib(b_lib,b_times,b_fluxes,b_powers);
  size_t num_steps = b_powers.size()>b_fluxes.size() ? b_powers.size() : b_fluxes.size();
  for(size_t i = 0; i < num_steps; i++){
    b_mat->add_step(b_times[i+1]-b_times[i]);
    b_mat->set_transition_matrix(b_mat->library()->newsp_transition_matrix_at(0));
    if(b_fluxes.size()==0){
      b_mat->set_flux(b_powers[i]/b_mat->power_factor_bos());
    }else if(b_powers.size()==0){
      b_mat->set_flux(b_fluxes[i]);
    }
    solver->set_transition_matrix( &*b_mat->transition_matrix() );
    Origen::SP_Vec_Dbl n0 = b_mat->amount_bos();
    Origen::SP_Vec_Dbl n1 = b_mat->amount_eos();
    solver->solve(*n0,b_mat->flux(),b_mat->dt(),&*n1);
    solver->clear();
  }
}

void cyclus2origen::solve(std::vector<double>& times, std::vector<double>& fluxes, std::vector<double>& powers){
  Check(b_mat!=NULL);
  Check(fluxes.size()>0||powers.size()>0);
  Check(fluxes.size()==(times.size()-1)||powers.size()==(times.size()-1));
  Origen::SP_Solver solver;
  ScaleUtils::IO::DB db;
  db.set<std::string>("solver","cram");
  solver = Origen::SolverSelector::get_solver(db);
  b_mat->set_solver(solver);
  prob_spec_lib(b_lib,times,fluxes,powers);
  size_t num_steps = powers.size()>fluxes.size() ? powers.size() : fluxes.size();
  for(size_t i = 0; i < num_steps; i++){
    b_mat->allocate_step();
    b_mat->set_transition_matrix(b_mat->library()->newsp_transition_matrix_at(0));
    if(fluxes.size()==0){
      b_mat->set_flux(powers[i]/b_mat->power_factor_bos());
    }else if(powers.size()==0){
      b_mat->set_flux(fluxes[i]);
    }
    b_mat->set_dt(times[i+1]-times[i]);

    solver->set_transition_matrix( &*b_mat->transition_matrix() );
    Origen::SP_Vec_Dbl n0 = b_mat->amount_bos();
    Origen::SP_Vec_Dbl n1 = b_mat->amount_eos();
    solver->solve(*n0,b_mat->flux(),b_mat->dt(),&*n1);
    solver->clear();
  }
}

void cyclus2origen::get_concentrations(std::vector<std::vector<double>> &concs_out) const{
  for(size_t i = 0; i < b_mat->ntimes(); i++){
    std::cout << "Getting concentrations at time " << b_mat->dt() << "." << std::endl;
    Origen::SP_DoubleList vals = b_mat->amount_at(i);
    std::vector<double> vals_vec;
    for(size_t j = 0; j < vals->size(); j++){
      vals_vec.push_back(vals->at(j));
    }
    concs_out.push_back(vals_vec);
  }
}

void cyclus2origen::get_concentrations_at(int p, std::vector<double> &concs_out) const{
  Check(p>=0);
  Check(p<=b_mat->ntimes());
  Origen::SP_DoubleList vals = b_mat->amount_at(p);
  for(size_t i = 0; i < vals->size(); i++){
    concs_out.push_back(vals->at(i));
  }
}

void cyclus2origen::get_concentrations_final(std::vector<double> &concs_out) const{
  Origen::SP_DoubleList vals = b_mat->amount_at(b_mat->ntimes()-1);
  for(size_t i = 0; i < vals->size(); i++){
    concs_out.push_back(vals->at(i));
  }
}

void cyclus2origen::get_masses(std::vector<std::vector<double> > &masses_out) const{
  Check(masses_out.empty());
  std::vector<std::vector<double> > concs;
  this->get_concentrations(concs);
  std::vector<int> ids;
  this->get_ids(ids);
  Check(ids.size()==concs[0].size());
  Origen::ConcentrationConverter cv;
  for(size_t i = 0; i < concs.size(); i++){
    std::vector<double> tmp;
    for(size_t j = 0; j < concs[0].size(); j++){
      tmp.push_back(cv.convert_to(Origen::ConcentrationUnit::KILOGRAMS,ids[i],Origen::ConcentrationUnit::CM_2_BARN,concs[i][j],b_vol));
    }
    masses_out.push_back(tmp);
  }
}

void cyclus2origen::get_masses_at(int p, std::vector<double> &masses_out) const{
  Check(masses_out.empty());
  std::vector<double> concs;
  this->get_concentrations_at(p,concs);
  std::vector<int> ids;
  this->get_ids(ids);
  Origen::ConcentrationConverter cv;
  for(size_t i = 0; i < concs.size(); i++){
    masses_out.push_back(cv.convert_to(Origen::ConcentrationUnit::KILOGRAMS,ids[i],Origen::ConcentrationUnit::CM_2_BARN,concs[i],b_vol));
  }
}

void cyclus2origen::get_masses_final(std::vector<double> &masses_out) const{
  Check(masses_out.empty());
  std::vector<double> concs;
  this->get_concentrations_final(concs);
  std::vector<int> ids;
  this->get_ids(ids);
  Origen::ConcentrationConverter cv;
  for(size_t i = 0; i < concs.size(); i++){
    masses_out.push_back(cv.convert_to(Origen::ConcentrationUnit::KILOGRAMS,ids[i],Origen::ConcentrationUnit::CM_2_BARN,concs[i],b_vol));
  }
}

void cyclus2origen::get_ids(std::vector<int> &ids_out) const{
  ids_out = *(b_mat->sizzzaaa_list());
}

void cyclus2origen::get_ids_zzzaaai(std::vector<int> &ids_out) const{
  for(auto id : *(b_mat->sizzzaaa_list())){
    ids_out.push_back(ScaleData::Utils::pizzzaaa_to_zzzaaai(id));
  }
}

void cyclus2origen::prob_spec_lib(Origen::SP_Library lib,std::vector<double> &times,std::vector<double> &fluxes,std::vector<double> &powers){
  if(b_time_units!=Origen::Time::DAYS){
    for(auto& time : times) time /= Origen::Time::factor(b_time_units,Origen::Time::DAYS);
  }
  if(fluxes.size()>0&&powers.size()==0){
    for(auto& flux : fluxes) powers.push_back(flux*b_mat->power_factor_bos());
  }else if(fluxes.size()>0&&powers.size()>0){
    Check(FALSE);
  }
  std::vector<double> burnups;
  Check(times.size()==powers.size()+1);
  for(size_t i = 0; i < powers.size(); i++){
// 1e3 factor arises from converting powers from watts to megawatts and mass from g to MT.
    burnups.push_back(1e3*powers[i]*(times[i+1]-times[i])/b_mat->initial_hm_mass());
  }
  Check(burnups.size()==powers.size());
  lib->interpolate_Interp1D(burnups);
}

}//namespace
