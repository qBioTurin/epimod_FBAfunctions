// ----------------------------------------------------------------------------
// general_functions_template.cpp
//
// Template for hypernode-specific C++ code.
// Use generate_cpp_from_arcs() in R to inject:
//   • V                (culture volume in mL)
//   • delta            (cell density in cell/mL)
//   • bacteria_names   (vector of "n_<id>")
//   • bacteriaBiomass_names (vector of "biomass_e_<id>")
// Everything else below is generic and should not be edited.
// ----------------------------------------------------------------------------

//
// Explanation and Implications:
// Role of general_functions_template.cpp: This file contains the generic C++ implementations for core biological processes like Starvation, Death, and Duplication that operate on the population and biomass levels of your hybrid models. These functions are designed to be general and parameterized (e.g., starvation_rate, half_life, duplication_rate, bacterium index), making them reusable across different microbial species within a multi-species model.
//
// Dynamic Customization: The comments within the template (// Simulation parameters (to be overwritten by R generator)) clearly indicate that certain sections (like V, delta, bacteria_names, bacteriaBiomass_names) are meant to be programmatically updated by an R script before compilation. This is the "injection" process mentioned.
//
// Interaction with R Scripts: R scripts in src/R/utils/ (like generate_cpp_from_arcs.R or setup_models.R) will likely read this template, inject model-specific details (like the actual names of n_ecs, n_cbd1, etc., and global V and delta from config_minimal_doublet.yaml), and then write out a model-specific C++ file (e.g., functions_minimal_doublet.cpp in your hypernodes/minimal_doublet/src/ was an example of such generated code). This generated file is then linked when model_generation compiles the solver.
//
// Clean Separation: Placing it in src/cpp/templates/ keeps it distinct from the directly editable C++ source (fba_automation/, shared_functions/) and from the generated model-specific C++ files (which will reside in models/EcCb/scripts/ or similar per case study).
//

#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>

// Bring the FBGLPK namespace into the current scope
using namespace FBGLPK;

// Simulation parameters (to be overwritten by R generator)
double V = 0.0;                // (mL)
long long int delta = 0;       // (cell/mL)
long long int max_total_bacteria = V * delta;
double total_bacteria = 0;     // total bacteria currently in simulation

// ----------------------------------------------------------------------------
// SharedMetabolite transition (stub implementation)
// ----------------------------------------------------------------------------
double SharedMetabolite(double *Value,
                        std::vector<class FBGLPK::LPprob>& vec_fluxb,
                        std::map<std::string,int>& NumTrans,
                        std::map<std::string,int>& NumPlaces,
                        const std::vector<std::string>& NameTrans,
                        const struct InfTr* Trans,
                        const int T,
                        const double& time
) {
  double rate = 0;
  return rate;
}

// ----------------------------------------------------------------------------
// Species identifiers (to be overwritten by R generator)
// ----------------------------------------------------------------------------
std::vector<std::string> bacteria_names = {
  // "n_ecs", "n_cbd1", ...
};

std::vector<std::string> bacteriaBiomass_names = {
  // "biomass_e_ecs", "biomass_e_cbd1", ...
};

// ----------------------------------------------------------------------------
// Update the total_bacteria global by summing floor(Value[place]) across species
// ----------------------------------------------------------------------------
void updateTotalBacteria(double *Value,
                         const std::map<std::string,int>& NumPlaces,
                         const std::vector<class FBGLPK::LPprob>& vec_fluxb) {
  total_bacteria = 0;
  for (size_t i = 0; i < bacteria_names.size(); ++i) {
    auto placeIt   = NumPlaces.find(bacteria_names[i]);
    auto biomassIt = NumPlaces.find(bacteriaBiomass_names[i]);
    if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
      total_bacteria += std::floor(Value[placeIt->second]);
    }
  }
}

// ----------------------------------------------------------------------------
// Debug file writers
// ----------------------------------------------------------------------------
static std::map<std::string, std::ofstream> outfileMap;
void WriteRateToFile(const std::string& bacteriumName, double time, double rate) {
  if (outfileMap.find(bacteriumName) == outfileMap.end()) {
    std::string fileName = "../logStarvationRate" + bacteriumName + ".csv";
    outfileMap[bacteriumName].open(fileName, std::ios::out);
    outfileMap[bacteriumName] << "Time,Rate" << std::endl;
  }
  outfileMap[bacteriumName] << std::fixed << std::setprecision(16)
                            << time << "," << rate << std::endl;
}

static std::map<std::string, std::ofstream> deathOutfileMap;
static std::map<std::string, std::ofstream> duplicationOutfileMap;

void WriteDeathRateToFile(const std::string& bacteriumName, double time, double rate) {
  if (deathOutfileMap.find(bacteriumName) == deathOutfileMap.end()) {
    std::string fileName = "../logDeathRate" + bacteriumName + ".csv";
    deathOutfileMap[bacteriumName].open(fileName, std::ios::out);
    deathOutfileMap[bacteriumName] << "Time,Rate" << std::endl;
  }
  deathOutfileMap[bacteriumName] << std::fixed << std::setprecision(16)
                                 << time << "," << rate << std::endl;
}

void WriteDuplicationRateToFile(const std::string& bacteriumName, double time, double rate) {
  if (duplicationOutfileMap.find(bacteriumName) == duplicationOutfileMap.end()) {
    std::string fileName = "../logDuplicationRate" + bacteriumName + ".csv";
    duplicationOutfileMap[bacteriumName].open(fileName, std::ios::out);
    duplicationOutfileMap[bacteriumName] << "Time,Rate" << std::endl;
  }
  duplicationOutfileMap[bacteriumName] << std::fixed << std::setprecision(16)
                                       << time << "," << rate << std::endl;
}

// ----------------------------------------------------------------------------
// Starvation rate function
// ----------------------------------------------------------------------------
double Starvation(double *Value,
                  std::vector<class FBGLPK::LPprob>& vec_fluxb,
                  std::map<std::string,int>& NumTrans,
                  std::map<std::string,int>& NumPlaces,
                  const std::vector<std::string>& NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time,
                  const double starvation_rate,
                  unsigned long bacterium) {
  
  if (bacterium >= bacteria_names.size()) {
    std::cerr << "Invalid bacterium index: " << bacterium << std::endl;
    return 0; 
  }
  
  const std::string& bacteriumName = bacteria_names[bacterium];
  auto placeIt   = NumPlaces.find(bacteriumName);
  auto biomassIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
  if (placeIt == NumPlaces.end() || biomassIt == NumPlaces.end()) {
    std::cerr << "Place index not found for bacterium index: " << bacterium << std::endl;
    return 0; 
  }
  
  double currentBiomass   = Value[biomassIt->second];
  double BioMin           = vec_fluxb[bacterium].getBioMin();
  double num_places_specie = std::floor(Value[placeIt->second]);
  double effectiveBiomass = (currentBiomass > BioMin)
    ? (currentBiomass - BioMin)
    : 0.0;
  double rate = (num_places_specie >= 1)
    ? (starvation_rate * effectiveBiomass)
    : 0.0;
  
  std::cout << "Starvation rate for " << bacteriumName
            << " at time " << time
            << " is " << rate << std::endl;
  // To debug to file:
  // WriteRateToFile(bacteriumName, time, rate);
  
  return rate;
}

// ----------------------------------------------------------------------------
// Death rate function
// ----------------------------------------------------------------------------
double Death(double *Value,
             std::vector<class FBGLPK::LPprob>& vec_fluxb,
             std::map<std::string,int>& NumTrans,
             std::map<std::string,int>& NumPlaces,
             const std::vector<std::string>& NameTrans,
             const struct InfTr* Trans,
             const int T,
             const double& time,
             const double half_life,
             unsigned long bacterium) {
  
  if (bacterium >= bacteria_names.size()) {
    std::cerr << "Invalid bacterium index in Death: " << bacterium << std::endl;
    return 0;
  }
  
  auto placeIt   = NumPlaces.find(bacteria_names[bacterium]);
  auto biomassIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
  double rate = 0;
  if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
    double currentBiomass = Value[biomassIt->second];
    double num_places_specie = std::floor(Value[placeIt->second]);
    double meanBiomass = vec_fluxb[bacterium].getBioMean();
    if (currentBiomass > 0 && num_places_specie >= 1) {
      rate = half_life * num_places_specie * (meanBiomass / currentBiomass);
    }
  } else {
    std::cerr << "Place not found for bacterium index: " << bacterium << std::endl;
  }
  
  const std::string& bacteriumName = bacteria_names[bacterium];
  std::cout << "Death rate for " << bacteriumName
            << " at time " << time
            << " is " << rate << std::endl;
  // WriteDeathRateToFile(bacteriumName, time, rate);
  
  return rate;
}

// ----------------------------------------------------------------------------
// Duplication rate function
// ----------------------------------------------------------------------------
double Duplication(double *Value,
                   std::vector<class FBGLPK::LPprob>& vec_fluxb,
                   std::map<std::string,int>& NumTrans,
                   std::map<std::string,int>& NumPlaces,
                   const std::vector<std::string>& NameTrans,
                   const struct InfTr* Trans,
                   const int T,
                   const double& time,
                   const double duplication_rate,
                   unsigned long bacterium) {
  
  if (bacterium >= bacteria_names.size()) {
    std::cerr << "Invalid bacterium index in Duplication: " << bacterium << std::endl;
    return 0;
  }
  
  auto placeIt   = NumPlaces.find(bacteria_names[bacterium]);
  auto biomassIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
  double rate = 0;
  if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
    double currentBiomass   = Value[biomassIt->second];
    double num_places_specie = std::floor(Value[placeIt->second]);
    double maxBiomass       = vec_fluxb[bacterium].getBioMax();
    
    updateTotalBacteria(Value, NumPlaces, vec_fluxb);
    
    if (num_places_specie >= 1 &&
        total_bacteria < max_total_bacteria) {
      rate = num_places_specie
      * duplication_rate
      * (currentBiomass / maxBiomass)
      * (1.0 - (total_bacteria / static_cast<double>(max_total_bacteria)));
    }
  } else {
    std::cerr << "Place not found for bacterium index: " << bacterium << std::endl;
  }
  
  const std::string& bacteriumName = bacteria_names[bacterium];
  std::cout << "Duplication rate for " << bacteriumName
            << " at time " << time
            << " is " << rate << std::endl;
  // WriteDuplicationRateToFile(bacteriumName, time, rate);
  
  return rate;
}

// ----------------------------------------------------------------------------
// End of general_functions_template.cpp
// ----------------------------------------------------------------------------
