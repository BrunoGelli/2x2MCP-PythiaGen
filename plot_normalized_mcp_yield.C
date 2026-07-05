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
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

namespace {

constexpr int kPseudoscalar = 0;
constexpr int kVector = 1;

constexpr int kLightMesons = 0;
constexpr int kCharmonium = 1;

constexpr double kAlpha = 1.0 / 137.0;

struct ParentPhysics {
    int pdg = 0;
    std::string name;
    int type = kPseudoscalar;
    double mass = 0.0;
    double reference_br = 0.0;
};

struct Row {
    double mass = 0.0;
    int pdg = 0;
    int mode = 0;
    long long events = 0;
    long long total = 0;
    long long accepted = 0;
    double sigma_sum = 0.0;
    double sigma_mean = 0.0;
    int nfiles = 1;
};

struct Point {
    double mass = 0.0;
    double acceptance = 0.0;
    double yield = 0.0;
    double old_proxy = 0.0;
};

struct PlotRange {
    double x_min = 1.0;
    double x_max = 10.0;
    double y_min = 1.0e-30;
    double y_max = 1.0;
};

std::map<int, ParentPhysics> parentTable()
{
    return {
        {111, {111, "#pi^{0}", kPseudoscalar, 0.1349768, 0.98823}},
        {221, {221, "#eta", kPseudoscalar, 0.5478620, 0.39410}},
        {331, {331, "#eta'", kPseudoscalar, 0.9577800, 0.0220}},
        {113, {113, "#rho^{0}", kVector, 0.7752600, 4.72e-5}},
        {223, {223, "#omega", kVector, 0.7826500, 7.36e-5}},
        {333, {333, "#phi", kVector, 1.0194610, 2.973e-4}},
        {443, {443, "J/#psi", kVector, 3.0969000, 5.971e-2}},
    };
}

std::vector<int> parentDrawOrder()
{
    return {111, 221, 331, 113, 223, 333, 443};
}

std::vector<int> graphColors()
{
    return {1, 2, 4, 8, 6, 46, 9};
}

bool isPositiveFinite(double value)
{
    return std::isfinite(value) && value > 0.0;
}

void updateRange(double value, double& min_value, double& max_value)
{
    if (!std::isfinite(value)) {
        return;
    }

    min_value = std::min(min_value, value);
    max_value = std::max(max_value, value);
}

void updatePositiveRange(double value, double& min_value, double& max_value)
{
    if (!isPositiveFinite(value)) {
        return;
    }

    updateRange(value, min_value, max_value);
}

std::pair<double, double> paddedLogRange(double min_value, double max_value)
{
    if (!isPositiveFinite(min_value) || !isPositiveFinite(max_value)) {
        return {1.0e-3, 1.0};
    }

    if (min_value == max_value) {
        return {min_value / 1.5, max_value * 1.5};
    }

    const double log_min = std::log10(min_value);
    const double log_max = std::log10(max_value);
    const double padding = 0.04 * (log_max - log_min);

    return {
        std::pow(10.0, log_min - padding),
        std::pow(10.0, log_max + padding),
    };
}

PlotRange makeLogYRange(double x_min, double x_max, double y_min, double y_max)
{
    auto padded_x = paddedLogRange(x_min, x_max);
    auto padded_y = paddedLogRange(y_min, y_max);

    return {
        padded_x.first,
        padded_x.second,
        padded_y.first,
        padded_y.second,
    };
}

double phaseSpace(const ParentPhysics& parent, double mchi)
{
    if (parent.mass <= 0.0 || 2.0 * mchi >= parent.mass) {
        return 0.0;
    }

    const double r = (mchi * mchi) / (parent.mass * parent.mass);

    if (parent.type == kPseudoscalar) {
        return std::pow(std::max(0.0, 1.0 - 4.0 * r), 3);
    }

    const double beta = std::sqrt(std::max(0.0, 1.0 - 4.0 * r));
    return beta * (1.0 + 2.0 * r);
}

double brExotic(const ParentPhysics& parent, double mchi, double epsilon)
{
    return (
        epsilon *
        epsilon *
        kAlpha *
        parent.reference_br *
        phaseSpace(parent, mchi)
    );
}

bool setLongLongBranch(TTree* tree, const char* name, long long* address)
{
    if (!tree->GetBranch(name)) {
        return false;
    }

    tree->SetBranchAddress(name, address);
    return true;
}

bool setIntBranch(TTree* tree, const char* name, int* address)
{
    if (!tree->GetBranch(name)) {
        return false;
    }

    tree->SetBranchAddress(name, address);
    return true;
}

bool setDoubleBranch(TTree* tree, const char* name, double* address)
{
    if (!tree->GetBranch(name)) {
        return false;
    }

    tree->SetBranchAddress(name, address);
    return true;
}

TTree* getSummaryTree(TFile* file)
{
    TTree* tree = static_cast<TTree*>(file->Get("mcp_summary"));

    if (!tree) {
        tree = static_cast<TTree*>(file->Get("aggregate_summary"));
    }

    return tree;
}

TGraph* makeGraph(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int marker_style,
    int color,
    int line_width = 1
) {
    TGraph* graph = new TGraph(
        static_cast<int>(x.size()),
        x.data(),
        y.data()
    );

    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineWidth(line_width);

    return graph;
}

void drawAcceptanceByEmitter(
    const std::map<int, std::vector<Point>>& by_parent,
    const std::map<int, ParentPhysics>& parents,
    double x_min,
    double x_max,
    double acceptance_max
) {
    auto x_range = paddedLogRange(x_min, x_max);

    const double y_max = (acceptance_max > 0.0)
        ? std::min(1.05, std::max(0.02, 1.25 * acceptance_max))
        : 1.0;

    TCanvas* canvas = new TCanvas(
        "c_acceptance_emitters_norm",
        "Acceptance by emitter",
        900,
        700
    );

    canvas->SetLogx();

    TH1* frame = canvas->DrawFrame(x_range.first, 0.0, x_range.second, y_max);
    frame->SetTitle(
        "Acceptance vs MCP mass;"
        "MCP mass [GeV];"
        "n_{MCP accepted}/n_{MCP total}"
    );

    TLegend* legend = new TLegend(0.62, 0.55, 0.88, 0.88);
    const auto order = parentDrawOrder();
    const auto colors = graphColors();

    for (std::size_t i = 0; i < order.size(); ++i) {
        const int pdg = order[i];
        auto found = by_parent.find(pdg);

        if (found == by_parent.end()) {
            continue;
        }

        std::vector<double> x;
        std::vector<double> y;

        for (const Point& point : found->second) {
            if (!isPositiveFinite(point.mass) || !std::isfinite(point.acceptance)) {
                continue;
            }

            x.push_back(point.mass);
            y.push_back(point.acceptance);
        }

        if (x.empty()) {
            continue;
        }

        TGraph* graph = makeGraph(
            x,
            y,
            20,
            colors[i % colors.size()]
        );

        graph->Draw("PL SAME");
        legend->AddEntry(graph, parents.at(pdg).name.c_str(), "lp");
    }

    legend->Draw();
    canvas->Modified();
    canvas->Update();
}

void drawYieldByEmitter(
    const std::map<int, std::vector<Point>>& by_parent,
    const std::map<int, ParentPhysics>& parents,
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    double pot,
    double epsilon
) {
    const PlotRange range = makeLogYRange(x_min, x_max, y_min, y_max);

    TCanvas* canvas = new TCanvas(
        "c_normalized_yield_emitters",
        "MCP accepted yield by emitter",
        900,
        700
    );

    canvas->SetLogx();
    canvas->SetLogy();

    TH1* frame = canvas->DrawFrame(
        range.x_min,
        range.y_min,
        range.x_max,
        range.y_max
    );

    frame->SetTitle(Form(
        "MCP accepted yield by emitter;"
        "MCP mass [GeV];"
        "Accepted MCPs, N_{POT}=%.2g, #epsilon=%.2g "
        "(placeholder branching normalization)",
        pot,
        epsilon
    ));

    TLegend* legend = new TLegend(0.50, 0.50, 0.88, 0.88);
    const auto order = parentDrawOrder();
    const auto colors = graphColors();

    for (std::size_t i = 0; i < order.size(); ++i) {
        const int pdg = order[i];
        auto found = by_parent.find(pdg);

        if (found == by_parent.end()) {
            continue;
        }

        std::vector<double> x;
        std::vector<double> y;

        for (const Point& point : found->second) {
            if (!isPositiveFinite(point.mass) || !isPositiveFinite(point.yield)) {
                continue;
            }

            x.push_back(point.mass);
            y.push_back(point.yield);
        }

        if (x.empty()) {
            continue;
        }

        TGraph* graph = makeGraph(
            x,
            y,
            21,
            colors[i % colors.size()]
        );

        graph->Draw("PL SAME");
        legend->AddEntry(graph, parents.at(pdg).name.c_str(), "lp");
    }

    legend->Draw();
    canvas->Modified();
    canvas->Update();
}

void drawTotalYield(
    const std::map<double, double>& total_yield_by_mass,
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    double pot,
    double epsilon
) {
    const PlotRange range = makeLogYRange(x_min, x_max, y_min, y_max);

    TCanvas* canvas = new TCanvas(
        "c_normalized_yield_total",
        "Total accepted MCP yield",
        900,
        700
    );

    canvas->SetLogx();
    canvas->SetLogy();

    TH1* frame = canvas->DrawFrame(
        range.x_min,
        range.y_min,
        range.x_max,
        range.y_max
    );

    frame->SetTitle(Form(
        "Total accepted MCP yield;"
        "MCP mass [GeV];"
        "Accepted MCPs, N_{POT}=%.2g, #epsilon=%.2g "
        "(not final sensitivity)",
        pot,
        epsilon
    ));

    std::vector<double> x;
    std::vector<double> y;

    for (const auto& entry : total_yield_by_mass) {
        if (!isPositiveFinite(entry.first) || !isPositiveFinite(entry.second)) {
            continue;
        }

        x.push_back(entry.first);
        y.push_back(entry.second);
    }

    if (!x.empty()) {
        TGraph* graph = makeGraph(x, y, 22, 1, 2);
        graph->Draw("PL SAME");
    }

    canvas->Modified();
    canvas->Update();
}

void drawOldProxyByEmitter(
    const std::map<int, std::vector<Point>>& by_parent,
    const std::map<int, ParentPhysics>& parents,
    double x_min,
    double x_max,
    double y_min,
    double y_max
) {
    const PlotRange range = makeLogYRange(x_min, x_max, y_min, y_max);

    TCanvas* canvas = new TCanvas(
        "c_old_proxy_compare",
        "Old proxy comparison",
        900,
        700
    );

    canvas->SetLogx();
    canvas->SetLogy();

    TH1* frame = canvas->DrawFrame(
        range.x_min,
        range.y_min,
        range.x_max,
        range.y_max
    );

    frame->SetTitle("Old proxy yield;MCP mass [GeV];old proxy");

    TLegend* legend = new TLegend(0.50, 0.50, 0.88, 0.88);
    const auto order = parentDrawOrder();
    const auto colors = graphColors();

    for (std::size_t i = 0; i < order.size(); ++i) {
        const int pdg = order[i];
        auto found = by_parent.find(pdg);

        if (found == by_parent.end()) {
            continue;
        }

        std::vector<double> x;
        std::vector<double> y;

        for (const Point& point : found->second) {
            if (!isPositiveFinite(point.mass) || !isPositiveFinite(point.old_proxy)) {
                continue;
            }

            x.push_back(point.mass);
            y.push_back(point.old_proxy);
        }

        if (x.empty()) {
            continue;
        }

        TGraph* graph = makeGraph(
            x,
            y,
            24,
            colors[i % colors.size()]
        );

        graph->Draw("PL SAME");
        legend->AddEntry(graph, parents.at(pdg).name.c_str(), "lp");
    }

    legend->Draw();
    canvas->Modified();
    canvas->Update();
}

}  // namespace

void plot_normalized_mcp_yield(
    const char* inputFile = "outputs/aggregate_summary_4h.root",
    double pot = 1.5e19,
    double epsilon = 1e-2,
    bool draw_old_proxy = false
) {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* file = TFile::Open(inputFile, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << "\n";
        return;
    }

    TTree* tree = getSummaryTree(file);

    if (!tree) {
        std::cerr
            << "Need mcp_summary or aggregate_summary tree in "
            << inputFile
            << "\n";
        return;
    }

    double mass = 0.0;
    double sigma_sum = 0.0;
    double sigma_mean = 0.0;

    int pdg = 0;
    int mode = kLightMesons;
    int nfiles = 1;

    long long events = 0;
    long long total = 0;
    long long accepted = 0;

    if (!setDoubleBranch(tree, "mcp_mass_GeV", &mass)) {
        setDoubleBranch(tree, "mcp_mass", &mass);
    }

    if (!setIntBranch(tree, "emitter_pdg", &pdg)) {
        setIntBranch(tree, "parent_pdg", &pdg);
    }

    setIntBranch(tree, "production_mode", &mode);
    setLongLongBranch(tree, "n_events_generated", &events);
    setLongLongBranch(tree, "n_mcp_total", &total);
    setLongLongBranch(tree, "n_mcp_accepted", &accepted);
    setDoubleBranch(tree, "sigma_gen_mb_sum", &sigma_sum);
    setDoubleBranch(tree, "sigma_gen_mb_mean", &sigma_mean);
    setIntBranch(tree, "n_input_files", &nfiles);

    std::vector<Row> rows;

    double soft_sigma_sum = 0.0;
    double charm_sigma_sum = 0.0;

    int n_soft = 0;
    int n_charm = 0;

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        sigma_sum = 0.0;
        sigma_mean = 0.0;
        nfiles = 1;

        tree->GetEntry(i);

        const double mean_sigma = (sigma_mean > 0.0)
            ? sigma_mean
            : ((nfiles > 0) ? sigma_sum / nfiles : sigma_sum);

        rows.push_back({
            mass,
            pdg,
            mode,
            events,
            total,
            accepted,
            sigma_sum,
            mean_sigma,
            nfiles,
        });

        if (mean_sigma > 0.0 && mode == kLightMesons) {
            soft_sigma_sum += mean_sigma;
            ++n_soft;
        }

        if (mean_sigma > 0.0 && mode == kCharmonium) {
            charm_sigma_sum += mean_sigma;
            ++n_charm;
        }
    }

    const double sigma_soft = n_soft ? soft_sigma_sum / n_soft : 1.0;
    const double sigma_charm = n_charm ? charm_sigma_sum / n_charm : 0.0;
    const double charm_scale = (sigma_soft > 0.0)
        ? sigma_charm / sigma_soft
        : 0.0;

    std::cout
        << "N_POT=" << pot
        << " epsilon=" << epsilon
        << " sigma_softqcd_mb_mean=" << sigma_soft
        << " sigma_charmonium_mb_mean=" << sigma_charm
        << " charmonium_process_scale=" << charm_scale
        << "\n";

    const auto parents = parentTable();

    std::map<int, std::vector<Point>> by_parent;
    std::map<double, double> total_yield_by_mass;

    double global_x_min = std::numeric_limits<double>::infinity();
    double global_x_max = -std::numeric_limits<double>::infinity();

    double acceptance_y_max = 0.0;

    double emitter_y_min = std::numeric_limits<double>::infinity();
    double emitter_y_max = -std::numeric_limits<double>::infinity();

    double total_y_min = std::numeric_limits<double>::infinity();
    double total_y_max = -std::numeric_limits<double>::infinity();

    double old_proxy_y_min = std::numeric_limits<double>::infinity();
    double old_proxy_y_max = -std::numeric_limits<double>::infinity();

    for (const Row& row : rows) {
        auto parent_it = parents.find(row.pdg);

        if (parent_it == parents.end()) {
            continue;
        }

        const ParentPhysics& parent = parent_it->second;

        const double acceptance = (row.total > 0)
            ? static_cast<double>(row.accepted) / static_cast<double>(row.total)
            : 0.0;

        const double branching_ratio = brExotic(parent, row.mass, epsilon);
        const double process_scale = (row.mode == kCharmonium)
            ? charm_scale
            : 1.0;

        const double accepted_per_event = (row.events > 0)
            ? static_cast<double>(row.accepted) / static_cast<double>(row.events)
            : 0.0;

        const double yield = pot * process_scale * accepted_per_event * branching_ratio;
        const double old_proxy = pot * accepted_per_event * branching_ratio;

        by_parent[row.pdg].push_back({
            row.mass,
            acceptance,
            yield,
            old_proxy,
        });

        total_yield_by_mass[row.mass] += yield;

        updatePositiveRange(row.mass, global_x_min, global_x_max);
        acceptance_y_max = std::max(acceptance_y_max, acceptance);
        updatePositiveRange(yield, emitter_y_min, emitter_y_max);
        updatePositiveRange(old_proxy, old_proxy_y_min, old_proxy_y_max);
    }

    for (auto& entry : by_parent) {
        auto& points = entry.second;
        std::sort(
            points.begin(),
            points.end(),
            [](const Point& left, const Point& right) {
                return left.mass < right.mass;
            }
        );
    }

    for (const auto& entry : total_yield_by_mass) {
        updatePositiveRange(entry.second, total_y_min, total_y_max);
    }

    if (!isPositiveFinite(global_x_min) || !isPositiveFinite(global_x_max)) {
        std::cerr << "No positive MCP masses found for known parent PDGs.\n";
        return;
    }

    drawAcceptanceByEmitter(
        by_parent,
        parents,
        global_x_min,
        global_x_max,
        acceptance_y_max
    );

    if (isPositiveFinite(emitter_y_min) && isPositiveFinite(emitter_y_max)) {
        drawYieldByEmitter(
            by_parent,
            parents,
            global_x_min,
            global_x_max,
            emitter_y_min,
            emitter_y_max,
            pot,
            epsilon
        );
    } else {
        std::cerr << "No positive per-emitter yields to draw on log-y scale.\n";
    }

    if (isPositiveFinite(total_y_min) && isPositiveFinite(total_y_max)) {
        drawTotalYield(
            total_yield_by_mass,
            global_x_min,
            global_x_max,
            total_y_min,
            total_y_max,
            pot,
            epsilon
        );
    } else {
        std::cerr << "No positive total yields to draw on log-y scale.\n";
    }

    if (
        draw_old_proxy &&
        isPositiveFinite(old_proxy_y_min) &&
        isPositiveFinite(old_proxy_y_max)
    ) {
        drawOldProxyByEmitter(
            by_parent,
            parents,
            global_x_min,
            global_x_max,
            old_proxy_y_min,
            old_proxy_y_max
        );
    }

    std::cout
        << "Note: branching constants are placeholder/model inputs; "
        << "do not interpret as final sensitivity.\n";
}
