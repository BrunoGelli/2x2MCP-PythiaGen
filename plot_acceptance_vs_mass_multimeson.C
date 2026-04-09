#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

namespace {

enum ParentType {
    kPseudoscalar = 0,
    kVector = 1,
    kDrellYan = 2
};

struct ParentPhysics {
    int pdg = 0;
    std::string name;
    int type = kPseudoscalar;
    double mass = 0.0;
    double reference_br = 0.0;
};

struct Aggregate {
    long long n_events_generated = 0;
    long long n_parent_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;

    double acceptance() const {
        return (n_mcp_total > 0)
            ? static_cast<double>(n_mcp_accepted) / static_cast<double>(n_mcp_total)
            : 0.0;
    }

    double acceptance_unc() const {
        if (n_mcp_total <= 0) return 0.0;
        const double p = acceptance();
        return std::sqrt(std::max(0.0, p * (1.0 - p) / static_cast<double>(n_mcp_total)));
    }

    double parent_yield_per_event() const {
        return (n_events_generated > 0)
            ? static_cast<double>(n_parent_total) / static_cast<double>(n_events_generated)
            : 0.0;
    }
};

std::map<int, ParentPhysics> parentTable() {
    return {
        {111, {111, "#pi^{0}",  kPseudoscalar, 0.1349770, 0.98823}},
        {221, {221, "#eta",  kPseudoscalar, 0.5478620, 0.39410}},
        {331, {331, "#eta'", kPseudoscalar, 0.9577800, 0.02220}},
        {113, {113, "#rho^{0}", kVector,       0.7752600, 4.72e-5}},
        {223, {223, "#omega",kVector,       0.7826500, 7.36e-5}},
        {333, {333, "#phi",  kVector,       1.0194610, 2.973e-4}},
        {443, {443, "J/#psi", kVector,       3.0969000, 5.971e-2}},
        {0,   {0,   "DY",   kDrellYan,     0.0,       1.0}}
    };
}

double phaseSpacePseudoscalar(double mChi, double mParent) {
    if (mParent <= 0.0 || 2.0 * mChi >= mParent) return 0.0;
    const double x = 4.0 * mChi * mChi / (mParent * mParent);
    return std::pow(std::max(0.0, 1.0 - x), 3);
}

double phaseSpaceVector(double mChi, double mParent) {
    if (mParent <= 0.0 || 2.0 * mChi >= mParent) return 0.0;
    const double x = 4.0 * mChi * mChi / (mParent * mParent);
    const double beta = std::sqrt(std::max(0.0, 1.0 - x));
    return beta * (1.0 + 2.0 * mChi * mChi / (mParent * mParent));
}

// Placeholder hook for DY weighting.
// Today this returns zero unless the tree already contains DY rows and the user
// chooses a nonzero normalization. The point is to keep the pipeline structure ready.
double dyWeightPlaceholder(double mChi, double epsilon, double dyNorm = 0.0) {
    (void)mChi;
    return dyNorm * epsilon * epsilon;
}

double productionWeight(const ParentPhysics& parent, double mChi, double epsilon, double alpha = 1.0 / 137.0, double dyNorm = 0.0) {
    if (parent.type == kPseudoscalar) {
        return epsilon * epsilon * alpha * parent.reference_br * phaseSpacePseudoscalar(mChi, parent.mass);
    }
    if (parent.type == kVector) {
        return epsilon * epsilon * alpha * parent.reference_br * phaseSpaceVector(mChi, parent.mass);
    }
    return dyWeightPlaceholder(mChi, epsilon, dyNorm);
}

} // namespace

void plot_acceptance_vs_mass_multimeson(const char* inputFile = "mcp_acceptance_results.root",
                                        double pot = 1.5e19,
                                        double epsilon = 1,
                                        double dyNorm = 0.0) {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open ROOT file: " << inputFile << "\n";
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("mcp_summary"));
    if (!tree) {
        std::cerr << "Could not find tree 'mcp_summary' in " << inputFile << "\n";
        file->Close();
        return;
    }

    double mcp_mass = 0.0;
    int parent_pdg = 0;
    int parent_type = -1;
    long long n_events_generated = 0;
    long long n_parent_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;

    tree->SetBranchAddress("mcp_mass", &mcp_mass);
    tree->SetBranchAddress("parent_pdg", &parent_pdg);
    tree->SetBranchAddress("parent_type", &parent_type);
    tree->SetBranchAddress("n_events_generated", &n_events_generated);
    tree->SetBranchAddress("n_parent_total", &n_parent_total);
    tree->SetBranchAddress("n_mcp_total", &n_mcp_total);
    tree->SetBranchAddress("n_mcp_accepted", &n_mcp_accepted);

    std::map<double, std::map<int, Aggregate>> agg;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        Aggregate& a = agg[mcp_mass][parent_pdg];
        a.n_events_generated += n_events_generated;
        a.n_parent_total += n_parent_total;
        a.n_mcp_total += n_mcp_total;
        a.n_mcp_accepted += n_mcp_accepted;
        (void)parent_type;
    }

    const std::map<int, ParentPhysics> parentInfo = parentTable();

    std::vector<double> massVec;
    std::vector<double> totalFluxVec;
    std::vector<double> totalAcceptanceVec;
    std::vector<double> totalAcceptanceErrVec;

    std::map<int, std::vector<double>> perParentMass;
    std::map<int, std::vector<double>> perParentAcceptance;
    std::map<int, std::vector<double>> perParentFlux;

    for (const auto& [mass, perParent] : agg) {
        massVec.push_back(mass);
        double weightedAccepted = 0.0;
        double weightedTotal = 0.0;
        double totalFlux = 0.0;

        for (const auto& [pdg, a] : perParent) {
            auto it = parentInfo.find(pdg);
            if (it == parentInfo.end()) continue;
            const ParentPhysics& parent = it->second;

            const double acceptance = a.acceptance();
            const double yieldPerEvent = a.parent_yield_per_event();
            const double weight = productionWeight(parent, mass, epsilon, 1.0 / 137.0, dyNorm);
            const double flux = (pdg == 0)
                ? pot * acceptance * weight
                : pot * yieldPerEvent * 2.0 * acceptance * weight;

            perParentMass[pdg].push_back(mass);
            perParentAcceptance[pdg].push_back(acceptance);
            perParentFlux[pdg].push_back(flux);

            weightedAccepted += static_cast<double>(a.n_mcp_accepted);
            weightedTotal += static_cast<double>(a.n_mcp_total);
            totalFlux += flux;
        }

        const double totalAcceptance = (weightedTotal > 0.0) ? (weightedAccepted / weightedTotal) : 0.0;
        const double totalAcceptanceErr = (weightedTotal > 0.0)
            ? std::sqrt(std::max(0.0, totalAcceptance * (1.0 - totalAcceptance) / weightedTotal))
            : 0.0;

        totalAcceptanceVec.push_back(totalAcceptance);
        totalAcceptanceErrVec.push_back(totalAcceptanceErr);
        totalFluxVec.push_back(totalFlux);
    }

    if (massVec.empty()) {
        std::cerr << "No entries found in mcp_summary.\n";
        file->Close();
        return;
    }

    TCanvas* c1 = new TCanvas("c_acceptance_total", "Total Acceptance vs MCP Mass", 900, 700);
    c1->SetLogx();
    TGraphErrors* gAccTotal = new TGraphErrors(
        static_cast<int>(massVec.size()), massVec.data(), totalAcceptanceVec.data(), nullptr, totalAcceptanceErrVec.data());
    gAccTotal->SetTitle("Total acceptance vs MCP mass;MCP mass [GeV];Acceptance fraction");
    gAccTotal->SetMarkerStyle(20);
    gAccTotal->SetLineWidth(2);
    gAccTotal->Draw("APL");

    TCanvas* c2 = new TCanvas("c_acceptance_parents", "Per-parent Acceptance vs MCP Mass", 900, 700);
    c2->SetLogx();
    TLegend* leg2 = new TLegend(0.60, 0.58, 0.88, 0.88);
    bool firstDraw = true;
    int color = 1;
    for (const auto& [pdg, masses] : perParentMass) {
        auto it = parentInfo.find(pdg);
        if (it == parentInfo.end()) continue;
        TGraph* g = new TGraph(static_cast<int>(masses.size()), masses.data(), perParentAcceptance[pdg].data());
        g->SetLineColor(color);
        g->SetMarkerColor(color);
        g->SetMarkerStyle(20 + (color % 5));
        g->SetLineWidth(2);
        g->SetTitle("Per-parent acceptance vs MCP mass;MCP mass [GeV];Acceptance fraction");
        g->Draw(firstDraw ? "APL" : "PL SAME");
        leg2->AddEntry(g, it->second.name.c_str(), "lp");
        firstDraw = false;
        ++color;
    }
    leg2->Draw();

    TCanvas* c3 = new TCanvas("c_flux_total", "Total MCP Flux vs MCP Mass", 900, 700);
    c3->SetLogx();
    c3->SetLogy();
    TGraph* gFluxTotal = new TGraph(static_cast<int>(massVec.size()), massVec.data(), totalFluxVec.data());
    gFluxTotal->SetTitle("Total MCP flux vs MCP mass;MCP mass [GeV];N_{#chi}/#epsilon^{2}/(1.5#times10^{19} POT)");
    gFluxTotal->SetMarkerStyle(20);
    gFluxTotal->SetLineWidth(2);
    gFluxTotal->Draw("APL");

    TCanvas* c4 = new TCanvas("c_flux_parents", "Per-parent MCP Flux vs MCP Mass", 900, 700);
    c4->SetLogx();
    c4->SetLogy();
    TLegend* leg4 = new TLegend(0.60, 0.58, 0.88, 0.88);
    firstDraw = false;
    color = 2;
    gFluxTotal->Draw("APL");
    leg4->AddEntry(gFluxTotal, "Total Flux", "lp");
    for (const auto& [pdg, masses] : perParentMass) {
        auto it = parentInfo.find(pdg);
        if (it == parentInfo.end()) continue;
        TGraph* g = new TGraph(static_cast<int>(masses.size()), masses.data(), perParentFlux[pdg].data());
        g->SetLineColor(color);
        g->SetMarkerColor(color);
        g->SetMarkerStyle(20);
        g->SetLineWidth(2);
        g->SetTitle("Per-parent MCP flux vs MCP mass;MCP mass [GeV];Flux at detector");
        g->Draw(firstDraw ? "APL" : "PL SAME");
        leg4->AddEntry(g, it->second.name.c_str(), "lp");
        firstDraw = false;
        ++color;
    }


    leg4->Draw();
    std::cout << "Loaded " << tree->GetEntries() << " summary rows from " << inputFile << "\n";
    std::cout << "Using epsilon = " << epsilon << ", POT = " << pot << ", DY normalization = " << dyNorm << "\n";
    std::cout << "Reminder: reference_br and the DY normalization placeholder are model inputs in parentTable()/dyWeightPlaceholder().\n";

    file->Close();

    for (int i = 0; i < massVec.size(); ++i)
    {
        cout << massVec[i] << " " << totalFluxVec[i] << endl;
    }
}
