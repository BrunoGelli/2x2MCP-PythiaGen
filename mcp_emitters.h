#ifndef MCP_EMITTERS_H
#define MCP_EMITTERS_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int kMcpPdg = 1000222;

enum class EmitterType : int { Pseudoscalar = 0, Vector = 1, DyReserved = 2 };
enum class ProductionMode : int { LightMesons = 0, Charmonium = 1, DyReserved = 2 };

struct EmitterConfig {
    std::string name;
    int pdg;
    EmitterType type;
    double mass_GeV;
    ProductionMode production_mode;
    std::string forced_decay_description;
    double reference_br_placeholder;
    bool enabled;
};

inline const std::vector<EmitterConfig>& emitterTable() {
    static const std::vector<EmitterConfig> table = {
        {"pi0",   111, EmitterType::Pseudoscalar, 0.1349768, ProductionMode::LightMesons, "pi0 -> gamma chi chibar", 0.98823, true},
        {"eta",   221, EmitterType::Pseudoscalar, 0.5478620, ProductionMode::LightMesons, "eta -> gamma chi chibar", 0.39410, true},
        {"etap",  331, EmitterType::Pseudoscalar, 0.9577800, ProductionMode::LightMesons, "eta' -> gamma chi chibar", 0.0220, true},
        {"rho0",  113, EmitterType::Vector,       0.7752600, ProductionMode::LightMesons, "rho0 -> chi chibar", 4.72e-5, true},
        {"omega", 223, EmitterType::Vector,       0.7826500, ProductionMode::LightMesons, "omega -> chi chibar", 7.36e-5, true},
        {"phi",   333, EmitterType::Vector,       1.0194610, ProductionMode::LightMesons, "phi -> chi chibar", 2.973e-4, true},
        {"jpsi",  443, EmitterType::Vector,       3.0969000, ProductionMode::Charmonium, "J/psi -> chi chibar", 5.971e-2, true}
    };
    return table;
}

inline std::string productionModeName(ProductionMode mode) {
    if (mode == ProductionMode::LightMesons) return "light_mesons";
    if (mode == ProductionMode::Charmonium) return "charmonium";
    return "dy_reserved";
}

inline ProductionMode parseProductionMode(const std::string& s) {
    if (s == "light_mesons" || s == "light") return ProductionMode::LightMesons;
    if (s == "charmonium" || s == "onia") return ProductionMode::Charmonium;
    if (s == "dy_reserved" || s == "dy") return ProductionMode::DyReserved;
    throw std::runtime_error("Unknown production mode: " + s);
}

inline const EmitterConfig* findEmitter(const std::string& name) {
    for (const auto& e : emitterTable()) if (e.name == name || std::to_string(e.pdg) == name) return &e;
    return nullptr;
}

inline bool isKinematicallyOpen(const EmitterConfig& e, double mchi) { return e.mass_GeV > 2.0 * mchi; }
inline double phaseSpacePseudoscalar(double mchi, double m) { if (!std::isfinite(m) || m <= 0 || 2*mchi >= m) return 0.0; return std::pow(std::max(0.0, 1.0 - 4*mchi*mchi/(m*m)), 3); }
inline double phaseSpaceVector(double mchi, double m) { if (!std::isfinite(m) || m <= 0 || 2*mchi >= m) return 0.0; return std::sqrt(std::max(0.0, 1.0 - 4*mchi*mchi/(m*m))); }

#endif
