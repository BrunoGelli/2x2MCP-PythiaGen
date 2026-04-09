#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "Pythia8/Pythia.h"

namespace {

// -----------------------------------------------------------------------------
// Physics and bookkeeping constants
// -----------------------------------------------------------------------------

// We keep using PDG 52 for the MCP state so the same particle ID can later be
// reused consistently if a validated Pythia dark-matter / DY card is plugged in.
constexpr int kMcpPdg = 52;
constexpr int kDyParentPdg = 0;

// Mutex used only for terminal output.
std::mutex gPrintMutex;

// -----------------------------------------------------------------------------
// Shared run state
// -----------------------------------------------------------------------------

enum class ParentType : int {
    Pseudoscalar = 0,
    Vector = 1,
    DrellYan = 2
};

enum class GeometryId : int {
    ArgoNeuT = 0,
    TwoByTwo = 1
};

enum class StopReason : int {
    Running = 0,
    AcceptedTargetReached = 1,
    EventCapReached = 2,
    WallTimeReached = 3
};

struct ProgressState {
    // Global counters for one mass point.
    std::atomic<long long> totalAccepted{0};
    std::atomic<long long> totalEvents{0};
    std::atomic<long long> totalNextFailures{0};
    std::atomic<int> finishedThreads{0};

    // Global stop flag shared across workers.
    std::atomic<bool> stopRequested{false};
    std::atomic<int> stopReason{static_cast<int>(StopReason::Running)};
};

struct ParentSpec {
    int pdg = 0;
    std::string name;
    ParentType type = ParentType::Pseudoscalar;
    std::vector<int> daughters;
    bool enabled = true;
    bool is_dy = false;
};

struct Counter {
    long long n_parent_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;
};

struct ThreadSummary {
    int thread_id = -1;
    int seed = 0;
    double mcp_mass = 0.0;
    std::string geometry_name;
    long long n_events_generated = 0;
    long long n_next_calls = 0;
    long long n_next_failures = 0;
    StopReason stop_reason = StopReason::Running;
    std::map<int, Counter> per_parent;
};

struct OutputRow {
    int thread_id = -1;
    int seed = 0;
    int geometry_id = -1;
    int parent_pdg = 0;
    int parent_type = -1;
    int stop_reason = -1;
    char geometry_name[16] = "";
    char parent_name[24] = "";
    char parent_type_name[16] = "";
    char stop_reason_name[24] = "";
    double mcp_mass = 0.0;
    long long n_events_generated = 0;
    long long n_next_calls = 0;
    long long n_next_failures = 0;
    long long n_parent_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;
    double acceptance_fraction = 0.0;
    double parent_yield_per_event = 0.0;
};

// -----------------------------------------------------------------------------
// Small helpers
// -----------------------------------------------------------------------------

std::string trim(const std::string& s) {
    const auto begin = s.find_first_not_of(" \t\n\r");
    if (begin == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t\n\r");
    return s.substr(begin, end - begin + 1);
}

std::vector<std::string> split(const std::string& input, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(input);
    std::string item;
    while (std::getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::string toCommand(const std::string& lhs, const std::string& rhs) {
    return lhs + " = " + rhs;
}

std::string parentTypeName(ParentType type) {
    switch (type) {
        case ParentType::Pseudoscalar: return "pseudoscalar";
        case ParentType::Vector: return "vector";
        case ParentType::DrellYan: return "dy";
    }
    return "unknown";
}

std::string stopReasonName(StopReason reason) {
    switch (reason) {
        case StopReason::Running: return "running";
        case StopReason::AcceptedTargetReached: return "accepted_target";
        case StopReason::EventCapReached: return "event_cap";
        case StopReason::WallTimeReached: return "wall_time";
    }
    return "unknown";
}

GeometryId parseGeometry(const std::string& geometry) {
    const std::string g = trim(geometry);
    if (g == "argoneut") return GeometryId::ArgoNeuT;
    if (g == "2x2" || g == "twobytwo") return GeometryId::TwoByTwo;
    throw std::runtime_error("Unknown geometry '" + geometry + "'. Use 'argoneut' or '2x2'.");
}

std::string geometryName(GeometryId geometry) {
    switch (geometry) {
        case GeometryId::ArgoNeuT: return "argoneut";
        case GeometryId::TwoByTwo: return "2x2";
    }
    return "unknown";
}

std::string formatDuration(double seconds) {
    if (seconds < 0.0 || !std::isfinite(seconds)) return "--:--:--";
    const long long s = static_cast<long long>(std::llround(seconds));
    const long long hh = s / 3600;
    const long long mm = (s % 3600) / 60;
    const long long ss = s % 60;
    std::ostringstream out;
    out << std::setfill('0') << std::setw(2) << hh << ':'
        << std::setw(2) << mm << ':'
        << std::setw(2) << ss;
    return out.str();
}

// Detector acceptance test.
//
// All dimensions are in meters in the projected detector plane, while z_det is
// expressed in the same effective coordinate system used by the old generator.
// The function remains intentionally simple because this code is for fast
// acceptance studies, not full detector simulation.
bool passesAcceptance(double px, double py, double pz, GeometryId geometry) {
    if (pz <= 0.0) return false;

    if (geometry == GeometryId::ArgoNeuT) {
        const double z_det = 1000.0;
        const double x_at_detector = px / pz * z_det;
        const double y_at_detector = py / pz * z_det;
        if (y_at_detector < -0.2 || y_at_detector > 0.2) return false;
        return (x_at_detector >= -0.24 && x_at_detector <= 0.24);
    }

    const double z_det = 1040.0;
    const double x_at_detector = px / pz * z_det;
    const double y_at_detector = py / pz * z_det;
    if (y_at_detector < -0.7 || y_at_detector > 0.7) return false;
    const bool in_left_box = (x_at_detector >= -0.65 && x_at_detector <= -0.05);
    const bool in_right_box = (x_at_detector >= 0.05 && x_at_detector <= 0.65);
    return in_left_box || in_right_box;
}

std::vector<ParentSpec> defaultParents() {
    return {
        {111, "pi0",   ParentType::Pseudoscalar, {22, kMcpPdg, -kMcpPdg}, true, false},
        {221, "eta",   ParentType::Pseudoscalar, {22, kMcpPdg, -kMcpPdg}, true, false},
        {331, "etap",  ParentType::Pseudoscalar, {22, kMcpPdg, -kMcpPdg}, true, false},
        {113, "rho0",  ParentType::Vector,       {kMcpPdg, -kMcpPdg},     true, false},
        {223, "omega", ParentType::Vector,       {kMcpPdg, -kMcpPdg},     true, false},
        {333, "phi",   ParentType::Vector,       {kMcpPdg, -kMcpPdg},     true, false},
        {443, "jpsi",  ParentType::Vector,       {kMcpPdg, -kMcpPdg},     true, false},
        {kDyParentPdg, "dy", ParentType::DrellYan, {}, false, true}
    };
}

std::map<int, ParentSpec> makeParentMap(const std::vector<ParentSpec>& parents) {
    std::map<int, ParentSpec> out;
    for (const auto& p : parents) out[p.pdg] = p;
    return out;
}

std::vector<ParentSpec> filterParents(const std::vector<ParentSpec>& parents, const std::string& parentListCsv) {
    if (parentListCsv.empty() || parentListCsv == "all") {
        std::vector<ParentSpec> enabled;
        for (const auto& p : parents) {
            if (!p.is_dy) enabled.push_back(p);
        }
        return enabled;
    }

    std::vector<ParentSpec> selected;
    const std::vector<std::string> tokens = split(parentListCsv, ',');
    for (const auto& token : tokens) {
        bool found = false;
        for (const auto& p : parents) {
            if (token == p.name || token == std::to_string(p.pdg)) {
                selected.push_back(p);
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Unknown parent selection token: '" + token + "'");
        }
    }
    return selected;
}

// Configure the millicharged particle placeholder state.
void configureMCP(Pythia8::Pythia& pythia, double mcpMass) {
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 1e12");

    pythia.readString("52:new = chi chibar 2 0 0 0 0");
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(9) << mcpMass;
        pythia.readString(toCommand("52:m0", ss.str()));
    }
    pythia.readString("52:tau0 = 1e12");
    pythia.readString("52:mayDecay = off");
}

// Force the meson decay table into the toy channels used by this acceptance
// study. This is intentionally explicit so future edits do not silently inherit
// unwanted default channels from the base Pythia tune.
void configureMesonParents(Pythia8::Pythia& pythia, const std::vector<ParentSpec>& parents) {
    for (const auto& parent : parents) {
        if (parent.is_dy) continue;
        pythia.readString(std::to_string(parent.pdg) + ":onMode = off");

        std::ostringstream channel;
        channel << parent.pdg << ":addChannel = 1 1.0 0";
        for (const int daughter : parent.daughters) channel << ' ' << daughter;
        pythia.readString(channel.str());
    }
}

// Walk the event ancestry to associate an MCP with one of the tracked parent
// mesons. The "same-PDG mother" clause handles cases where a particle is copied
// through intermediate shower/history bookkeeping steps before the actual parent
// of interest appears.
int identifyTrackedParent(const Pythia8::Event& event, int particleIndex, const std::map<int, ParentSpec>& parentMap) {
    int current = particleIndex;
    for (int depth = 0; depth < 8; ++depth) {
        if (current <= 0 || current >= event.size()) break;
        const int mother1 = event[current].mother1();
        const int mother2 = event[current].mother2();

        const auto testMother = [&](int idx) -> int {
            if (idx <= 0 || idx >= event.size()) return 0;
            const int pdg = event[idx].id();
            if (parentMap.count(pdg)) return pdg;
            if (parentMap.count(-pdg)) return -pdg;
            return 0;
        };

        int found = testMother(mother1);
        if (found != 0) return found;
        found = testMother(mother2);
        if (found != 0) return found;

        if (mother1 > 0 && mother1 < event.size() && std::abs(event[mother1].id()) == std::abs(event[current].id())) {
            current = mother1;
            continue;
        }
        break;
    }
    return 0;
}

std::string makeProgressBar(double fraction, int width = 28) {
    fraction = std::max(0.0, std::min(1.0, fraction));
    const int filled = static_cast<int>(std::round(fraction * width));
    std::string bar;
    bar.reserve(width + 2);
    bar.push_back('[');
    for (int i = 0; i < width; ++i) bar.push_back(i < filled ? '#' : '-');
    bar.push_back(']');
    return bar;
}

void requestGlobalStop(ProgressState* progress, StopReason reason) {
    bool expected = false;
    if (progress->stopRequested.compare_exchange_strong(expected, true)) {
        progress->stopReason.store(static_cast<int>(reason), std::memory_order_relaxed);
    }
}

void maybeUpdateGlobalStopFromGoals(ProgressState* progress,
                                    long long acceptedGoalTotal,
                                    long long eventCapTotal) {
    if (acceptedGoalTotal > 0 &&
        progress->totalAccepted.load(std::memory_order_relaxed) >= acceptedGoalTotal) {
        requestGlobalStop(progress, StopReason::AcceptedTargetReached);
        return;
    }

    if (eventCapTotal > 0 &&
        progress->totalEvents.load(std::memory_order_relaxed) >= eventCapTotal) {
        requestGlobalStop(progress, StopReason::EventCapReached);
    }
}

void printConfigurationSummary(int seed,
                               int nThreads,
                               long long acceptedTargetTotal,
                               int maxEventsPerThread,
                               int maxWallSeconds,
                               double mcpMass,
                               GeometryId geometry,
                               const std::vector<ParentSpec>& parents,
                               const std::string& outFileName,
                               const std::string& beamConfig,
                               const std::string& momentumConfig) {
    std::lock_guard<std::mutex> lock(gPrintMutex);
    std::cout << "Configuration summary\n"
              << "  seed base:              " << seed << '\n'
              << "  threads:                " << nThreads << '\n'
              << "  accepted target (mass): " << acceptedTargetTotal << '\n'
              << "  max events / thread:    " << maxEventsPerThread << '\n'
              << "  max wall time [s]:      " << maxWallSeconds << '\n'
              << "  MCP mass [GeV]:         " << std::fixed << std::setprecision(6) << mcpMass << '\n'
              << "  geometry:               " << geometryName(geometry) << '\n'
              << "  parents:                ";
    for (size_t i = 0; i < parents.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << parents[i].name;
    }
    std::cout << '\n'
              << "  output file:            " << outFileName << '\n'
              << "  beam config:            " << beamConfig << '\n'
              << "  momentum config:        " << momentumConfig << "\n\n";
}

void progressMonitor(
    ProgressState* progress,
    int nThreads,
    long long acceptedGoalTotal,
    long long eventCapTotal,
    int maxWallSeconds,
    double mcpMass,
    const std::string& geometry,
    const std::string& parentsLabel,
    std::chrono::steady_clock::time_point startTime) {

    using namespace std::chrono_literals;

    while (progress->finishedThreads.load(std::memory_order_relaxed) < nThreads) {
        const auto now = std::chrono::steady_clock::now();
        const double elapsedSec = std::chrono::duration<double>(now - startTime).count();

        if (maxWallSeconds > 0 && elapsedSec >= static_cast<double>(maxWallSeconds)) {
            requestGlobalStop(progress, StopReason::WallTimeReached);
        }

        maybeUpdateGlobalStopFromGoals(progress, acceptedGoalTotal, eventCapTotal);

        const long long accepted = progress->totalAccepted.load(std::memory_order_relaxed);
        const long long events = progress->totalEvents.load(std::memory_order_relaxed);
        const long long failures = progress->totalNextFailures.load(std::memory_order_relaxed);
        const int finished = progress->finishedThreads.load(std::memory_order_relaxed);

        const double acceptedFrac = (acceptedGoalTotal > 0)
            ? std::min(1.0, static_cast<double>(accepted) / static_cast<double>(acceptedGoalTotal))
            : 0.0;
        const double eventsPerSec = (elapsedSec > 0.0) ? static_cast<double>(events) / elapsedSec : 0.0;
        const double acceptedPerSec = (elapsedSec > 0.0) ? static_cast<double>(accepted) / elapsedSec : 0.0;
        const double etaAcceptedSec = (acceptedPerSec > 0.0 && acceptedGoalTotal > accepted)
            ? static_cast<double>(acceptedGoalTotal - accepted) / acceptedPerSec
            : std::numeric_limits<double>::quiet_NaN();
        const double etaWallSec = (maxWallSeconds > 0)
            ? std::max(0.0, static_cast<double>(maxWallSeconds) - elapsedSec)
            : std::numeric_limits<double>::quiet_NaN();

        std::string etaText = "--:--:--";
        if (std::isfinite(etaAcceptedSec) && std::isfinite(etaWallSec)) {
            etaText = formatDuration(std::min(etaAcceptedSec, etaWallSec));
        } else if (std::isfinite(etaAcceptedSec)) {
            etaText = formatDuration(etaAcceptedSec);
        } else if (std::isfinite(etaWallSec)) {
            etaText = formatDuration(etaWallSec);
        }

        {
            std::lock_guard<std::mutex> lock(gPrintMutex);
            std::cout << '\r'
                      << "mass=" << std::fixed << std::setprecision(6) << mcpMass << " GeV"
                      << " | geom=" << geometry
                      << " | parents=" << parentsLabel
                      << " | acc " << makeProgressBar(acceptedFrac)
                      << ' ' << accepted << '/' << acceptedGoalTotal
                      << " | evt " << events
                      << " | fail " << failures
                      << " | rate " << std::fixed << std::setprecision(1) << eventsPerSec << " evt/s"
                      << " | accRate " << acceptedPerSec << " /s"
                      << " | elapsed " << formatDuration(elapsedSec)
                      << " | ETA " << etaText
                      << " | done " << finished << '/' << nThreads
                      << std::flush;
        }

        std::this_thread::sleep_for(500ms);
    }

    const auto now = std::chrono::steady_clock::now();
    const double elapsedSec = std::chrono::duration<double>(now - startTime).count();

    const long long accepted = progress->totalAccepted.load(std::memory_order_relaxed);
    const long long events = progress->totalEvents.load(std::memory_order_relaxed);
    const long long failures = progress->totalNextFailures.load(std::memory_order_relaxed);
    const double acceptedFrac = (acceptedGoalTotal > 0)
        ? std::min(1.0, static_cast<double>(accepted) / static_cast<double>(acceptedGoalTotal))
        : 0.0;

    std::lock_guard<std::mutex> lock(gPrintMutex);
    std::cout << '\r'
              << "mass=" << std::fixed << std::setprecision(6) << mcpMass << " GeV"
              << " | geom=" << geometry
              << " | parents=" << parentsLabel
              << " | acc " << makeProgressBar(acceptedFrac)
              << ' ' << accepted << '/' << acceptedGoalTotal
              << " | evt " << events
              << " | fail " << failures
              << " | elapsed " << formatDuration(elapsedSec)
              << " | done " << nThreads << '/' << nThreads
              << " | stop=" << stopReasonName(static_cast<StopReason>(progress->stopReason.load(std::memory_order_relaxed)))
              << "  done.\n";
}

ThreadSummary runThread(
    int baseSeed,
    long long acceptedGoalTotal,
    long long eventCapTotal,
    int maxWallSeconds,
    int threadId,
    double mcpMass,
    GeometryId geometry,
    const std::vector<ParentSpec>& parents,
    const std::string& beamConfig,
    const std::string& momentumConfig,
    ProgressState* progress,
    std::chrono::steady_clock::time_point startTime) {

    ThreadSummary summary;
    summary.thread_id = threadId;
    summary.seed = baseSeed + threadId;
    summary.mcp_mass = mcpMass;
    summary.geometry_name = geometryName(geometry);

    const std::map<int, ParentSpec> parentMap = makeParentMap(parents);
    for (const auto& p : parents) summary.per_parent[p.pdg] = Counter{};

    Pythia8::Pythia pythia("", false);
    pythia.readFile(beamConfig);
    pythia.readFile(momentumConfig);
    pythia.readString("Random:setSeed = on");
    pythia.readString(toCommand("Random:seed", std::to_string(summary.seed)));

    // Silence the default periodic Pythia progress chatter. We provide our own
    // monitor thread instead.
    pythia.readString("Next:numberCount = 0");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");

    configureMCP(pythia, mcpMass);
    configureMesonParents(pythia, parents);
    pythia.init();

    while (!progress->stopRequested.load(std::memory_order_relaxed)) {
        const auto now = std::chrono::steady_clock::now();
        const double elapsedSec = std::chrono::duration<double>(now - startTime).count();
        if (maxWallSeconds > 0 && elapsedSec >= static_cast<double>(maxWallSeconds)) {
            requestGlobalStop(progress, StopReason::WallTimeReached);
            break;
        }

        maybeUpdateGlobalStopFromGoals(progress, acceptedGoalTotal, eventCapTotal);
        if (progress->stopRequested.load(std::memory_order_relaxed)) break;

        ++summary.n_next_calls;
        if (!pythia.next()) {
            ++summary.n_next_failures;
            progress->totalNextFailures.fetch_add(1, std::memory_order_relaxed);
            continue;
        }

        ++summary.n_events_generated;
        const long long newTotalEvents = progress->totalEvents.fetch_add(1, std::memory_order_relaxed) + 1;

        for (int i = 0; i < pythia.event.size(); ++i) {
            const auto& particle = pythia.event[i];
            const int id = particle.id();

            auto parentIt = summary.per_parent.find(id);
            if (parentIt != summary.per_parent.end()) {
                parentIt->second.n_parent_total++;
            }

            if (std::abs(id) != kMcpPdg) continue;
            const int parentPdg = identifyTrackedParent(pythia.event, i, parentMap);
            if (parentPdg == 0) continue;

            Counter& counter = summary.per_parent[parentPdg];
            counter.n_mcp_total++;
            if (passesAcceptance(particle.px(), particle.py(), particle.pz(), geometry)) {
                counter.n_mcp_accepted++;
                const long long newAccepted = progress->totalAccepted.fetch_add(1, std::memory_order_relaxed) + 1;
                if (acceptedGoalTotal > 0 && newAccepted >= acceptedGoalTotal) {
                    requestGlobalStop(progress, StopReason::AcceptedTargetReached);
                }
            }
        }

        if (eventCapTotal > 0 && newTotalEvents >= eventCapTotal) {
            requestGlobalStop(progress, StopReason::EventCapReached);
        }
    }

    summary.stop_reason = static_cast<StopReason>(progress->stopReason.load(std::memory_order_relaxed));
    progress->finishedThreads.fetch_add(1, std::memory_order_relaxed);
    return summary;
}

void connectOrCreateBranch(TTree* tree, const char* name, void* address, const char* leaflist) {
    if (TBranch* branch = tree->GetBranch(name)) {
        branch->SetAddress(address);
    } else {
        tree->Branch(name, address, leaflist);
    }
}

void writeSummaryTree(
    const std::vector<ThreadSummary>& summaries,
    const std::vector<ParentSpec>& parents,
    const std::string& outFileName,
    GeometryId geometry) {

    const std::map<int, ParentSpec> parentMap = makeParentMap(parents);

    TFile outFile(outFileName.c_str(), "UPDATE");
    TTree* tree = dynamic_cast<TTree*>(outFile.Get("mcp_summary"));
    if (!tree) {
        tree = new TTree("mcp_summary", "MCP acceptance summary by thread and parent");
    }

    OutputRow row;
    connectOrCreateBranch(tree, "thread_id", &row.thread_id, "thread_id/I");
    connectOrCreateBranch(tree, "seed", &row.seed, "seed/I");
    connectOrCreateBranch(tree, "geometry_id", &row.geometry_id, "geometry_id/I");
    connectOrCreateBranch(tree, "geometry_name", row.geometry_name, "geometry_name/C");
    connectOrCreateBranch(tree, "parent_pdg", &row.parent_pdg, "parent_pdg/I");
    connectOrCreateBranch(tree, "parent_type", &row.parent_type, "parent_type/I");
    connectOrCreateBranch(tree, "parent_name", row.parent_name, "parent_name/C");
    connectOrCreateBranch(tree, "parent_type_name", row.parent_type_name, "parent_type_name/C");
    connectOrCreateBranch(tree, "stop_reason", &row.stop_reason, "stop_reason/I");
    connectOrCreateBranch(tree, "stop_reason_name", row.stop_reason_name, "stop_reason_name/C");
    connectOrCreateBranch(tree, "mcp_mass", &row.mcp_mass, "mcp_mass/D");
    connectOrCreateBranch(tree, "n_events_generated", &row.n_events_generated, "n_events_generated/L");
    connectOrCreateBranch(tree, "n_next_calls", &row.n_next_calls, "n_next_calls/L");
    connectOrCreateBranch(tree, "n_next_failures", &row.n_next_failures, "n_next_failures/L");
    connectOrCreateBranch(tree, "n_parent_total", &row.n_parent_total, "n_parent_total/L");
    connectOrCreateBranch(tree, "n_mcp_total", &row.n_mcp_total, "n_mcp_total/L");
    connectOrCreateBranch(tree, "n_mcp_accepted", &row.n_mcp_accepted, "n_mcp_accepted/L");
    connectOrCreateBranch(tree, "acceptance_fraction", &row.acceptance_fraction, "acceptance_fraction/D");
    connectOrCreateBranch(tree, "parent_yield_per_event", &row.parent_yield_per_event, "parent_yield_per_event/D");

    for (const auto& summary : summaries) {
        for (const auto& [parentPdg, counter] : summary.per_parent) {
            auto it = parentMap.find(parentPdg);
            if (it == parentMap.end()) continue;
            const ParentSpec& parent = it->second;

            row = OutputRow{};
            row.thread_id = summary.thread_id;
            row.seed = summary.seed;
            row.geometry_id = static_cast<int>(geometry);
            std::snprintf(row.geometry_name, sizeof(row.geometry_name), "%s", summary.geometry_name.c_str());
            row.parent_pdg = parentPdg;
            row.parent_type = static_cast<int>(parent.type);
            std::snprintf(row.parent_name, sizeof(row.parent_name), "%s", parent.name.c_str());
            std::snprintf(row.parent_type_name, sizeof(row.parent_type_name), "%s", parentTypeName(parent.type).c_str());
            row.stop_reason = static_cast<int>(summary.stop_reason);
            std::snprintf(row.stop_reason_name, sizeof(row.stop_reason_name), "%s", stopReasonName(summary.stop_reason).c_str());
            row.mcp_mass = summary.mcp_mass;
            row.n_events_generated = summary.n_events_generated;
            row.n_next_calls = summary.n_next_calls;
            row.n_next_failures = summary.n_next_failures;
            row.n_parent_total = counter.n_parent_total;
            row.n_mcp_total = counter.n_mcp_total;
            row.n_mcp_accepted = counter.n_mcp_accepted;
            row.acceptance_fraction = (counter.n_mcp_total > 0)
                ? static_cast<double>(counter.n_mcp_accepted) / static_cast<double>(counter.n_mcp_total)
                : 0.0;
            row.parent_yield_per_event = (summary.n_events_generated > 0)
                ? static_cast<double>(counter.n_parent_total) / static_cast<double>(summary.n_events_generated)
                : 0.0;
            tree->Fill();
        }
    }

    outFile.cd();
    tree->Write("", TObject::kOverwrite);
    outFile.Close();
}

void printUsage(const char* argv0) {
    std::cout << "Usage:\n"
              << "  " << argv0
              << " [seed] [nThreads] [acceptedTargetTotal] [mcpMass_GeV] [geometry]"
              << " [parentList=all] [maxEventsPerThread=1000000] [outFile=mcp_acceptance_results.root]"
              << " [beamConfig=beam.config] [momentumConfig=momentum.config] [maxWallSeconds=0]\n\n"
              << "Notes:\n"
              << "  - acceptedTargetTotal is the total accepted-event target for the whole mass point,\n"
              << "    not per thread.\n"
              << "  - maxWallSeconds = 0 disables the wall-time safety stop.\n\n"
              << "Examples:\n"
              << "  " << argv0 << " 1000 8 500 0.02 2x2\n"
              << "  " << argv0 << " 1000 8 500 0.02 argoneut pi0,eta,etap 1000000 results.root beam.config momentum.config 1800\n"
              << "  " << argv0 << " 1000 4 2000 0.50 2x2 rho0,omega,phi,jpsi\n";
}

std::string joinParentNames(const std::vector<ParentSpec>& parents) {
    std::ostringstream ss;
    for (size_t i = 0; i < parents.size(); ++i) {
        if (i) ss << ',';
        ss << parents[i].name;
    }
    return ss.str();
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 6) {
        printUsage(argv[0]);
        return 1;
    }

    try {
        const int seed = std::stoi(argv[1]);
        const int nThreads = std::stoi(argv[2]);
        const long long acceptedTargetTotal = std::stoll(argv[3]);
        const double mcpMass = std::stod(argv[4]);
        const GeometryId geometry = parseGeometry(trim(argv[5]));
        const std::string parentList = (argc > 6) ? trim(argv[6]) : "all";
        const int maxEventsPerThread = (argc > 7) ? std::stoi(argv[7]) : 1000000;
        const std::string outFileName = (argc > 8) ? argv[8] : "mcp_acceptance_results.root";
        const std::string beamConfig = (argc > 9) ? argv[9] : "beam.config";
        const std::string momentumConfig = (argc > 10) ? argv[10] : "momentum.config";
        const int maxWallSeconds = (argc > 11) ? std::stoi(argv[11]) : 0;

        if (nThreads <= 0 || acceptedTargetTotal <= 0 || maxEventsPerThread <= 0) {
            throw std::runtime_error("nThreads, acceptedTargetTotal, and maxEventsPerThread must be positive.");
        }
        if (mcpMass <= 0.0) {
            throw std::runtime_error("mcpMass must be positive.");
        }
        if (maxWallSeconds < 0) {
            throw std::runtime_error("maxWallSeconds must be non-negative.");
        }

        const std::vector<ParentSpec> parents = filterParents(defaultParents(), parentList);
        if (parents.empty()) {
            throw std::runtime_error("No parents enabled after parsing parent list.");
        }

        const long long eventCapTotal = static_cast<long long>(nThreads) * static_cast<long long>(maxEventsPerThread);
        const std::string parentsLabel = joinParentNames(parents);

        printConfigurationSummary(seed,
                                  nThreads,
                                  acceptedTargetTotal,
                                  maxEventsPerThread,
                                  maxWallSeconds,
                                  mcpMass,
                                  geometry,
                                  parents,
                                  outFileName,
                                  beamConfig,
                                  momentumConfig);

        std::vector<std::thread> workers;
        std::vector<ThreadSummary> summaries(nThreads);
        ProgressState progress;
        const auto startTime = std::chrono::steady_clock::now();

        std::thread monitor(progressMonitor,
                            &progress,
                            nThreads,
                            acceptedTargetTotal,
                            eventCapTotal,
                            maxWallSeconds,
                            mcpMass,
                            geometryName(geometry),
                            parentsLabel,
                            startTime);

        for (int i = 0; i < nThreads; ++i) {
            workers.emplace_back([&, i]() {
                summaries[i] = runThread(seed,
                                         acceptedTargetTotal,
                                         eventCapTotal,
                                         maxWallSeconds,
                                         i,
                                         mcpMass,
                                         geometry,
                                         parents,
                                         beamConfig,
                                         momentumConfig,
                                         &progress,
                                         startTime);
            });
        }

        for (auto& worker : workers) worker.join();
        if (monitor.joinable()) monitor.join();

        writeSummaryTree(summaries, parents, outFileName, geometry);

        std::cout << "Wrote summary tree 'mcp_summary' to " << outFileName << "\n";
        std::cout << "Final stop reason: "
                  << stopReasonName(static_cast<StopReason>(progress.stopReason.load(std::memory_order_relaxed)))
                  << "\n";
        std::cout << "Note: DY is already supported in the output schema and plotting layer, but generator-side Pythia DY steering is still intentionally left off until you choose a validated local process card.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return 2;
    }
}
