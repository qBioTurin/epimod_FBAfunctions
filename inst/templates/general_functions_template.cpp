// ----------------------------------------------------------------------------
/* general_functions_template.cpp

   Template for hypernode-specific C++ code.
   Use generate_cpp_from_arcs() in R to inject:
     • V                (culture volume in mL)
     • delta            (cell density in cell/mL)
     • bacteria_names   (vector of "n_<id>")
     • bacteriaBiomass_names (vector of "biomass_e_<id>")
   Everything else below is generic and should not be edited.

   Explanation and Implications:
   This file provides generic C++ implementations for core biological processes
   such as Starvation, Death, and Duplication acting on populations and biomass
   in hybrid models. The functions are parameterized (e.g., starvation_rate,
   half_life, duplication_rate, bacterium index) to be reused across multiple
   microbial species within multi-species models.

   Dynamic Customization:
   Sections explicitly marked "// Simulation parameters (to be overwritten by R
   generator)" and "// Species identifiers (to be overwritten by R generator)"
   are meant to be programmatically updated by an R script prior to compilation.
   The R generator injects actual values for V, delta, and the species vectors.

   Numerical Robustness and Positivity Strategy (step-size free):
   ODE solvers (e.g., LSODE) do not enforce state non-negativity. This file
   ensures practical non-negativity for populations by:
     • Soft capping *outflow* rates by the available amount n (Death).
     • Soft capping *inflow* rates by the remaining capacity (K − N_total) (Duplication).
     • Smooth gating near extinction (n ≈ 1) to avoid RHS discontinuities.
     • Optional post-step clamping utilities as a last-resort safety net.

   Notation in comments below:
     • n             := continuous population level stored in Value[place]
     • places        := floor(n + ε) for n > 1 (0 otherwise) — effective count
     • bm            := current biomass (Value at biomass place)
     • BioMin        := minimal biomass threshold (from FBA)
     • meanBm, maxBm := mean / maximal biomass (from FBA)
     • K             := carrying capacity = V * delta (cells)
     • N_total       := sum of per-species effective counts
     • fade(n)       := smooth gating in [0,1], 0 for n≤1 and 1 for n≥2

   Key soft-cap (soft minimum) used for static, step-size-free bounding:
     soft_cap(r, a) = (r * a) / (r + a + ε)  ≈  min(r, a) smoothly
   which ensures the applied rate never exceeds the available "a" but matches
   the raw rate r when r is small compared to a, and does not need the time step.
*/
// ----------------------------------------------------------------------------

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
#include <limits>

// Bring the FBGLPK namespace into the current scope
using namespace FBGLPK;

// ---------------------------------------------------------------------------
// Positivity helpers, clamping, and diagnostics (static, step-size free)
// ---------------------------------------------------------------------------
namespace {

constexpr double EPS = 1e-9;

// Extinction threshold and smoothing width.
// For n ≤ N_MIN the species is considered extinct for rate purposes.
// Between n ∈ (N_MIN, N_MIN+SMOOTH_WIDTH) rates fade smoothly.
constexpr double N_MIN = 1.0;
constexpr double SMOOTH_WIDTH = 1.0; // defines the 1→2 fading range

// Positive part utility: max(x, 0)
inline double pospart(double x) { return (x > 0.0 ? x : 0.0); }

// Smooth cubic fade on n:
//   fade(n) = 0                 for n ≤ 1
//   fade(n) = 1                 for n ≥ 2
//   fade(n) = 3t^2 − 2t^3       for n ∈ (1,2), t = (n−1)/1
inline double fade_above_one(double n) {
    if (n <= N_MIN) return 0.0;
    if (n >= N_MIN + SMOOTH_WIDTH) return 1.0;
    const double t = (n - N_MIN) / SMOOTH_WIDTH; // t ∈ (0,1)
    return t * t * (3.0 - 2.0 * t);
}

inline double clamp01(double x) {
    return (x < 0.0) ? 0.0 : (x > 1.0 ? 1.0 : x);
}

// Non-negative biomass; if negative due to noise, log and clamp to 0.
inline double safeNonNeg(double raw, const std::string& name, const char* kind) {
    if (raw < -EPS) {
        std::cerr << "[Warn][NegValue] " << kind << "=" << name
                  << " raw=" << std::setprecision(16) << raw << std::endl;
        return 0.0;
    }
    return (raw <= EPS) ? 0.0 : raw;
}

// Robust integer-like "places" from continuous n.
// For n ≤ 1, return 0 (extinct); otherwise floor(n + ε).
inline long long safeCountFromPlace(double raw_n, const std::string& species) {
    if (raw_n <= N_MIN) return 0LL;
    long long n = static_cast<long long>(std::floor(raw_n + EPS));
    return (n < 0) ? 0LL : n;
}

// Carrying capacity K = V * delta (overflow-protected).
inline long long compute_capacity(double V, long long delta, long long& K_out) {
    long double prod = static_cast<long double>(V) * static_cast<long double>(delta);
    if (prod < 0) prod = 0;
    const long double LLMAX = static_cast<long double>(std::numeric_limits<long long>::max());
    if (prod > LLMAX) {
        std::cerr << "[Warn][Capacity] V*delta overflow; clamped to LLONG_MAX\n";
        prod = LLMAX;
    }
    K_out = static_cast<long long>(prod);
    return K_out;
}

// Soft minimum (soft cap) that smoothly limits a rate r by availability a:
//   soft_cap(r, a) = (r * a) / (r + a + ε)  ≈ min(r, a)
// Step-size free: no time step appears.
inline double soft_cap(double r, double a) {
    if (r <= 0.0 || a <= 0.0) return 0.0;
    return (r * a) / (r + a + EPS);
}

} // namespace

// ----------------------------------------------------------------------------
// Simulation parameters (to be overwritten by R generator)
// ----------------------------------------------------------------------------
double V = 0.0;                 // culture volume (mL)        [to be injected]
long long int delta = 0;        // cell density (cells/mL)    [to be injected]
long long int max_total_bacteria = V * delta;  // carrying capacity K = V * delta
double total_bacteria = 0.0;    // current total effective count N_total

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
// Debug file writers (optional)
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
// SharedMetabolite transition (stub)
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
  double rate = 0.0;
  return rate;
}

// ----------------------------------------------------------------------------
// updateTotalBacteria
// ----------------------------------------------------------------------------
// Accumulates the total effective count:
//   N_total = Σ_i places_i, where places_i = floor(n_i + ε) for n_i > 1; else 0.
// This is robust to small negative noise and values near the extinction
// threshold.
void updateTotalBacteria(double *Value,
                         const std::map<std::string,int>& NumPlaces,
                         const std::vector<class FBGLPK::LPprob>& vec_fluxb) {
  (void)vec_fluxb; // kept for signature compatibility
  total_bacteria = 0.0;
  for (size_t i = 0; i < bacteria_names.size(); ++i) {
    auto pIt = NumPlaces.find(bacteria_names[i]);
    auto bIt = NumPlaces.find(bacteriaBiomass_names[i]);
    if (pIt != NumPlaces.end() && bIt != NumPlaces.end()) {
      total_bacteria += static_cast<double>(
          safeCountFromPlace(Value[pIt->second], bacteria_names[i])
      );
    }
  }
}

// ----------------------------------------------------------------------------
// Starvation
// ----------------------------------------------------------------------------
// Biological intent: consume surplus biomass relative to a minimal threshold.
// Continuous model (pre-cap):
//   rate_starvation_raw = s * max(0, bm − BioMin) * fade(n)
// Applied rate:
//   rate_starvation = rate_starvation_raw
// Notes:
//   • Starvation affects biomass only; it is gated near extinction by fade(n)
//     to reduce RHS discontinuities around n ≈ 1.
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
        std::cerr << "[Starvation][Error] Invalid index: " << bacterium << std::endl;
        return 0.0;
    }

    const std::string& sp = bacteria_names[bacterium];
    auto pIt = NumPlaces.find(sp);
    auto bIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
    if (pIt == NumPlaces.end() || bIt == NumPlaces.end()) {
        std::cerr << "[Starvation][Error] Place not found for species: " << sp << std::endl;
        return 0.0;
    }

    const double n_cont = Value[pIt->second];
    const double fade   = fade_above_one(n_cont);
    const double bm     = safeNonNeg(Value[bIt->second], bacteriaBiomass_names[bacterium].c_str(), "biomass");
    double BioMin       = vec_fluxb[bacterium].getBioMin();
    if (BioMin < 0.0) BioMin = 0.0;

    const double effBm  = (bm > BioMin) ? (bm - BioMin) : 0.0;
    double rate         = starvation_rate * effBm * fade;

    if (!std::isfinite(rate) || rate < 0.0) rate = 0.0;

    std::cout << "[Starvation] sp="<<sp<<", t="<<time
              <<", n="<<n_cont<<", bm="<<bm<<", BioMin="<<BioMin
              <<", fade="<<fade<<", rate="<<rate<<std::endl;

    return rate;
}

// ----------------------------------------------------------------------------
// Death (with static, step-size-free soft cap)
// ----------------------------------------------------------------------------
// Biological intent (raw):
//   rate_death_raw = λ * places * (meanBm / bm) * fade(n)
// Soft cap by availability (no time step):
//   rate_death = soft_cap(rate_death_raw, n+)
// where n+ = max(n, 0). This ensures the applied outflow cannot exceed the
// available amount, without requiring knowledge of the solver step.
// Definitions:
//   λ        := half_life parameter (model-specific scaling)
//   places   := floor(n + ε) for n > 1, else 0
//   meanBm   := mean biomass from FBA
//   bm       := current biomass
//   fade(n)  := smooth gating as above
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
        std::cerr << "[Death][Error] Invalid index: " << bacterium << std::endl;
        return 0.0;
    }

    const std::string& sp = bacteria_names[bacterium];
    auto pIt = NumPlaces.find(sp);
    auto bIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
    if (pIt == NumPlaces.end() || bIt == NumPlaces.end()) {
        std::cerr << "[Death][Error] Place not found for species: " << sp << std::endl;
        return 0.0;
    }

    const double n_cont = Value[pIt->second];
    if (n_cont <= N_MIN) {
        std::cout << "[Death][Extinct] sp="<<sp<<", t="<<time<<", n="<<n_cont<<" => rate=0\n";
        return 0.0;
    }

    const double fade      = fade_above_one(n_cont);
    const long long places = safeCountFromPlace(n_cont, sp);
    const double bm        = safeNonNeg(Value[bIt->second], bacteriaBiomass_names[bacterium].c_str(), "biomass");
    double meanBm          = vec_fluxb[bacterium].getBioMean();
    if (meanBm < 0.0) meanBm = 0.0;

    // Raw model rate (unbounded)
    double rate_raw = (bm > 0.0 && places >= 1)
        ? (half_life * static_cast<double>(places) * (meanBm / bm))
        : 0.0;

    rate_raw *= fade;

    // Static soft cap by availability n+ (step-size free)
    const double n_avail = pospart(n_cont);
    double rate = soft_cap(rate_raw, n_avail);

    if (!std::isfinite(rate) || rate < 0.0) rate = 0.0;

    std::cout << "[Death] sp="<<sp<<", t="<<time
              <<", n="<<n_cont<<", bm="<<bm<<", meanBm="<<meanBm
              <<", places="<<places<<", fade="<<fade
              <<", rate_raw="<<rate_raw<<", n_avail="<<n_avail
              <<", rate="<<rate<<std::endl;

    return rate;
}

// ----------------------------------------------------------------------------
// Duplication (with static, step-size-free soft cap)
// ----------------------------------------------------------------------------
// Biological intent (raw):
//   rate_dup_raw = places * μ * (bm / maxBm) * (1 − N_total / K) * fade(n)
// Soft cap by remaining capacity (no time step):
//   rate_dup = soft_cap(rate_dup_raw, (K − N_total)+ )
// ensuring the applied inflow cannot exceed remaining capacity in rate form,
// without knowing the solver time step.
// Definitions:
//   μ         := duplication_rate
//   maxBm     := maximal biomass from FBA
//   K         := V * delta (carrying capacity)
//   N_total   := total effective count (sum across species)
//   fade(n)   := smooth gating as above
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
        std::cerr << "[Duplication][Error] Invalid index: " << bacterium << std::endl;
        return 0.0;
    }

    const std::string& sp = bacteria_names[bacterium];
    auto pIt = NumPlaces.find(sp);
    auto bIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
    if (pIt == NumPlaces.end() || bIt == NumPlaces.end()) {
        std::cerr << "[Duplication][Error] Place not found for species: " << sp << std::endl;
        return 0.0;
    }

    const double n_cont = Value[pIt->second];
    if (n_cont <= N_MIN) {
        std::cout << "[Duplication][Extinct] sp="<<sp<<", t="<<time<<", n="<<n_cont<<" => rate=0\n";
        return 0.0;
    }

    const double fade      = fade_above_one(n_cont);
    const long long places = safeCountFromPlace(n_cont, sp);
    const double bm        = safeNonNeg(Value[bIt->second], bacteriaBiomass_names[bacterium].c_str(), "biomass");
    const double maxBm     = vec_fluxb[bacterium].getBioMax();
    if (maxBm <= EPS) return 0.0;

    // Compute capacity K and total effective count
    compute_capacity(V, delta, max_total_bacteria);
    updateTotalBacteria(Value, NumPlaces, vec_fluxb);
    if (max_total_bacteria <= 0) return 0.0;

    const double cap_ratio = clamp01(static_cast<double>(total_bacteria) /
                                     static_cast<double>(max_total_bacteria));
    const double logistic  = 1.0 - cap_ratio; // ∈ [0,1]

    // Raw model rate (unbounded)
    double rate_raw = (places >= 1)
        ? (static_cast<double>(places) * duplication_rate * (bm / maxBm) * logistic * fade)
        : 0.0;

    // Static soft cap by remaining capacity (step-size free)
    const double capacity_left = pospart(static_cast<double>(max_total_bacteria) - total_bacteria);
    double rate = soft_cap(rate_raw, capacity_left);

    if (!std::isfinite(rate) || rate < 0.0) rate = 0.0;

    std::cout << "[Duplication] sp="<<sp<<", t="<<time
              <<", n="<<n_cont<<", bm="<<bm<<", maxBm="<<maxBm
              <<", places="<<places<<", logistic="<<logistic<<", fade="<<fade
              <<", rate_raw="<<rate_raw<<", capacity_left="<<capacity_left
              <<", K="<<max_total_bacteria<<", N_total="<<total_bacteria
              <<", rate="<<rate<<std::endl;

    return rate;
}

// ----------------------------------------------------------------------------
// Clamp utilities (optional; solver-loop safety net)
// ----------------------------------------------------------------------------
// ClampNonNegativeState:
//   For each species i, if n_i < 0 then set n_i := 0, and similarly for biomass.
//   This is a last-resort guard against numerical undershoot after a solver step.
void ClampNonNegativeState(double *Value,
                           const std::map<std::string,int>& NumPlaces) {
    for (size_t i = 0; i < bacteria_names.size(); ++i) {
        auto pIt = NumPlaces.find(bacteria_names[i]);
        if (pIt != NumPlaces.end()) {
            if (Value[pIt->second] < 0.0) {
                std::cerr << "[Fixup] " << bacteria_names[i]
                          << " < 0 -> 0 (raw=" << Value[pIt->second] << ")\n";
                Value[pIt->second] = 0.0;
            }
        }
        auto bIt = NumPlaces.find(bacteriaBiomass_names[i]);
        if (bIt != NumPlaces.end()) {
            if (Value[bIt->second] < 0.0) {
                std::cerr << "[Fixup] " << bacteriaBiomass_names[i]
                          << " < 0 -> 0 (raw=" << Value[bIt->second] << ")\n";
                Value[bIt->second] = 0.0;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// End of general_functions_template.cpp
// ----------------------------------------------------------------------------

