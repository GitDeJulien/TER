#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <toml/toml.hpp>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& other = toml::find(config, "other");
   this->_mesh_name = toml::find<std::string>(other, "mesh");
   // this->_mu = toml::find<double>(other, "mu");
   // this->_numerical_flux_choice = toml::find<std::string>(other, "numerical_flux");
   this->_results = toml::find<std::string>(other, "results");
   this->_xb = toml::find<double>(other, "xb");
   this->_h0g = toml::find<double>(other, "h0g");
   this->_h0d = toml::find<double>(other, "h0d");

   // Physics
   const auto& physics = toml::find(config, "physics");
   this->_g = toml::find<double>(other, "g");
   this->_rho = toml::find<double>(other, "rho");

   // Time
   const auto& time = toml::find(config, "time");
   this->_t0 = toml::find<double>(time, "t0");
   this->_tfinal = toml::find<double>(time, "tfinal");
   this->_dt = toml::find<double>(time, "dt");
   this->_cfl = toml::find<double>(time, "cfl");
   this->_scheme = toml::find<std::string>(time, "scheme");

   // Boundary conditions
   const auto& BC = toml::find(config, "BC");
   this->_BC_ref = toml::find<std::vector<int> >(BC, "ref");
   this->_BC_type = toml::find<std::vector<std::string> >(BC, "BC");

   // Créer le dossier de résultats
   system(("mkdir -p ./" +this->_results).c_str());
   // Supprimer les anciens résultats
   system(("rm -f ./" +this->_results + "/*.vtk").c_str());
   // Copier le fichier de données dans le dossier résultats
   system(("cp -r ./" + this->_file_name + " ./"
   + this->_results + "/params.txt").c_str());
}

#define _DATA_FILE_CPP
#endif
