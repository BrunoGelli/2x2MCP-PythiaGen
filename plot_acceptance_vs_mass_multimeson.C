#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

namespace {

struct ParentInfo {
    int pdg;
    int type;  // 0 = pseudoscalar, 1 = vector
    int mode;
    const char* name;
    double mass;
    double br;
};

std::map<int, ParentInfo> make_parent_table()
{
    return {
        {111, {111, 0, 0, "#pi^{0}", 0.1349768, 0.98823}},
        {221, {221, 0, 0, "#eta", 0.547862, 0.39410}},
        {331, {331, 0, 0, "#eta'", 0.95778, 0.022}},
        {113, {113, 1, 0, "#rho^{0}", 0.77526, 4.72e-5}},
        {223, {223, 1, 0, "#omega", 0.78265, 7.36e-5}},
        {333, {333, 1, 0, "#phi", 1.019461, 2.973e-4}},
        {443, {443, 1, 1, "J/#psi", 3.0969, 5.971e-2}},
    };
}

double phase_space(int type, double m_mcp, double m_parent)
{
    if (m_parent <= 0.0 || 2.0 * m_mcp >= m_parent) {
        return 0.0;
    }

    const double x = 4.0 * m_mcp * m_mcp / (m_parent * m_parent);
    const double one_minus_x = std::max(0.0, 1.0 - x);

    if (type == 0) {
        return std::pow(one_minus_x, 3);
    }

    return std::sqrt(one_minus_x);
}

struct Accumulator {
    long long n_events_generated = 0;
    long long n_emitter_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;
};

void update_x_range(double x, double& x_min, double& x_max)
{
    if (x <= 0.0) {
        return;  // needed because both canvases use log-x
    }

    x_min = std::min(x_min, x);
    x_max = std::max(x_max, x);
}

void expand_log_range(double& x_min, double& x_max, double pad_fraction = 0.05)
{
    if (x_min <= 0.0 || x_max <= 0.0) {
        return;
    }

    if (x_max <= x_min) {
        x_min /= 1.2;
        x_max *= 1.2;
    }

    const double log_min = std::log10(x_min);
    const double log_max = std::log10(x_max);
    const double pad = pad_fraction * (log_max - log_min);

    x_min = std::pow(10.0, log_min - pad);
    x_max = std::pow(10.0, log_max + pad);
}

TGraph* make_graph(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int marker_style,
    int color
)
{
    TGraph* graph = new TGraph(static_cast<int>(x.size()), x.data(), y.data());

    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);

    return graph;
}

TGraph* make_positive_y_graph(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int marker_style,
    int color
)
{
    std::vector<double> x_positive;
    std::vector<double> y_positive;

    x_positive.reserve(x.size());
    y_positive.reserve(y.size());

    for (std::size_t i = 0; i < x.size(); ++i) {
        if (x[i] > 0.0 && y[i] > 0.0) {
            x_positive.push_back(x[i]);
            y_positive.push_back(y[i]);
        }
    }

    if (x_positive.empty()) {
        return nullptr;
    }

    return make_graph(x_positive, y_positive, marker_style, color);
}

}  // namespace

void plot_acceptance_vs_mass_multimeson(
    const char* inputFile = "outputs/aggregate_summary.root",
    double pot = 1.5e19,
    double epsilon = 1e-2,
    bool sum_charmonium_with_light = false
)
{
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << "\n";
        return;
    }

    TTree* tree = static_cast<TTree*>(file->Get("mcp_summary"));
    if (!tree) {
        tree = static_cast<TTree*>(file->Get("aggregate_summary"));
    }

    if (!tree) {
        std::cerr
            << "Need mcp_summary/aggregate_summary tree. "
            << "CSV-only aggregate can be plotted after conversion or use raw ROOT files.\n";
        return;
    }

    double mass = 0.0;
    int pdg = 0;
    int production_mode = 0;
    long long n_events_generated = 0;
    long long n_emitter_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;

    const char* mass_branch =
        tree->GetBranch("mcp_mass_GeV") ? "mcp_mass_GeV" : "mcp_mass";

    const char* parent_branch =
        tree->GetBranch("emitter_pdg") ? "emitter_pdg" : "parent_pdg";

    const char* emitter_count_branch =
        tree->GetBranch("n_emitter_total") ? "n_emitter_total" : "n_parent_total";

    tree->SetBranchAddress(mass_branch, &mass);
    tree->SetBranchAddress(parent_branch, &pdg);

    if (tree->GetBranch("production_mode")) {
        tree->SetBranchAddress("production_mode", &production_mode);
    }

    tree->SetBranchAddress("n_events_generated", &n_events_generated);
    tree->SetBranchAddress(emitter_count_branch, &n_emitter_total);
    tree->SetBranchAddress("n_mcp_total", &n_mcp_total);
    tree->SetBranchAddress("n_mcp_accepted", &n_mcp_accepted);

    std::map<std::pair<double, int>, Accumulator> aggregate;
    std::map<int, int> modes;

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);

        std::pair<double, int> key(mass, pdg);
        Accumulator& entry = aggregate[key];

        entry.n_events_generated += n_events_generated;
        entry.n_emitter_total += n_emitter_total;
        entry.n_mcp_total += n_mcp_total;
        entry.n_mcp_accepted += n_mcp_accepted;

        modes[pdg] = production_mode;
    }

    const std::map<int, ParentInfo> parent_info = make_parent_table();

    std::map<int, std::vector<double> > masses_by_parent;
    std::map<int, std::vector<double> > acceptance_by_parent;
    std::map<int, std::vector<double> > flux_proxy_by_parent;

    double x_min = std::numeric_limits<double>::infinity();
    double x_max = 0.0;
    double acceptance_max = 0.0;
    double flux_min_positive = std::numeric_limits<double>::infinity();
    double flux_max = 0.0;

    for (std::map<std::pair<double, int>, Accumulator>::const_iterator it = aggregate.begin();
         it != aggregate.end();
         ++it) {
        const double m_mcp = it->first.first;
        const int parent_pdg = it->first.second;
        const Accumulator& counts = it->second;

        std::map<int, ParentInfo>::const_iterator info_it = parent_info.find(parent_pdg);
        if (info_it == parent_info.end()) {
            continue;
        }

        const ParentInfo& info = info_it->second;

        const double acceptance =
            counts.n_mcp_total > 0
                ? static_cast<double>(counts.n_mcp_accepted) / counts.n_mcp_total
                : 0.0;

        const double yield_per_event =
            counts.n_events_generated > 0
                ? static_cast<double>(counts.n_emitter_total) / counts.n_events_generated
                : 0.0;

        const double weight =
            epsilon * epsilon * info.br * phase_space(info.type, m_mcp, info.mass);

        const double flux_proxy = pot * yield_per_event * 2.0 * acceptance * weight;

        masses_by_parent[parent_pdg].push_back(m_mcp);
        acceptance_by_parent[parent_pdg].push_back(acceptance);
        flux_proxy_by_parent[parent_pdg].push_back(flux_proxy);

        update_x_range(m_mcp, x_min, x_max);
        acceptance_max = std::max(acceptance_max, acceptance);

        if (flux_proxy > 0.0) {
            flux_min_positive = std::min(flux_min_positive, flux_proxy);
            flux_max = std::max(flux_max, flux_proxy);
        }
    }

    if (!std::isfinite(x_min) || x_max <= 0.0) {
        std::cerr << "No positive MCP masses found; cannot draw log-x plots.\n";
        return;
    }

    expand_log_range(x_min, x_max);

    if (acceptance_max <= 0.0) {
        acceptance_max = 1.0;
    }

    const double acceptance_y_max = std::min(1.0, 1.20 * acceptance_max);

    TCanvas* c_acceptance =
        new TCanvas("c_acceptance_emitters", "Acceptance by emitter", 900, 700);

    c_acceptance->SetLogx();

    TH1* acceptance_frame = c_acceptance->DrawFrame(
        x_min,
        0.0,
        x_max,
        acceptance_y_max
    );

    acceptance_frame->SetTitle(
        "Acceptance vs MCP mass;MCP mass [GeV];accepted / produced MCP"
    );

    TLegend* acceptance_legend = new TLegend(0.62, 0.55, 0.88, 0.88);

    int color = 1;
    for (std::map<int, std::vector<double> >::const_iterator it = masses_by_parent.begin();
         it != masses_by_parent.end();
         ++it) {
        const int parent_pdg = it->first;

        TGraph* graph = make_graph(
            it->second,
            acceptance_by_parent[parent_pdg],
            20,
            color
        );

        graph->Draw("PL SAME");
        acceptance_legend->AddEntry(graph, parent_info.at(parent_pdg).name, "lp");

        ++color;
    }

    acceptance_legend->Draw();
    c_acceptance->Modified();
    c_acceptance->Update();

    if (flux_max <= 0.0 || !std::isfinite(flux_min_positive)) {
        std::cerr << "No positive flux-proxy values found; skipping log-y flux plot.\n";
    } else {
        TCanvas* c_flux =
            new TCanvas("c_flux_proxy_emitters", "Flux proxy by emitter", 900, 700);

        c_flux->SetLogx();
        c_flux->SetLogy();

        double flux_y_min = flux_min_positive;
        double flux_y_max = flux_max;
        expand_log_range(flux_y_min, flux_y_max, 0.10);

        TH1* flux_frame = c_flux->DrawFrame(
            x_min,
            flux_y_min,
            x_max,
            flux_y_max
        );

        flux_frame->SetTitle(
            "Accepted MCP flux proxy;MCP mass [GeV];proxy (placeholder norm)"
        );

        TLegend* flux_legend = new TLegend(0.62, 0.55, 0.88, 0.88);

        color = 1;
        for (std::map<int, std::vector<double> >::const_iterator it = masses_by_parent.begin();
             it != masses_by_parent.end();
             ++it) {
            const int parent_pdg = it->first;

            TGraph* graph = make_positive_y_graph(
                it->second,
                flux_proxy_by_parent[parent_pdg],
                21,
                color
            );

            if (graph) {
                graph->Draw("PL SAME");
                flux_legend->AddEntry(graph, parent_info.at(parent_pdg).name, "lp");
            }

            ++color;
        }

        flux_legend->Draw();
        c_flux->Modified();
        c_flux->Update();
    }

    if (!sum_charmonium_with_light) {
        std::cout
            << "Light-meson and charmonium modes are kept separate by default; "
            << "pass true only for compatible toy sums.\n";
    }
}
