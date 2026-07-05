#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
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

struct ParentPhysics { int pdg; std::string name; int type; double mass; double reference_br; };
struct Row { double mass=0; int pdg=0, mode=0; long long events=0, total=0, accepted=0; double sigma_sum=0, sigma_mean=0; int nfiles=1; };
struct Point { double mass=0, acceptance=0, yield=0, old_proxy=0; };

std::map<int, ParentPhysics> parentTable() {
    return {
        {111, {111, "#pi^{0}", kPseudoscalar, 0.1349768, 0.98823}},
        {221, {221, "#eta", kPseudoscalar, 0.5478620, 0.39410}},
        {331, {331, "#eta'", kPseudoscalar, 0.9577800, 0.0220}},
        {113, {113, "#rho^{0}", kVector, 0.7752600, 4.72e-5}},
        {223, {223, "#omega", kVector, 0.7826500, 7.36e-5}},
        {333, {333, "#phi", kVector, 1.0194610, 2.973e-4}},
        {443, {443, "J/#psi", kVector, 3.0969000, 5.971e-2}}
    };
}

double phaseSpace(const ParentPhysics& p, double mchi) {
    if (p.mass <= 0.0 || 2.0 * mchi >= p.mass) return 0.0;
    const double r = (mchi * mchi) / (p.mass * p.mass);
    if (p.type == kPseudoscalar) return std::pow(std::max(0.0, 1.0 - 4.0 * r), 3);
    const double beta = std::sqrt(std::max(0.0, 1.0 - 4.0 * r));
    return beta * (1.0 + 2.0 * r);
}

double brExotic(const ParentPhysics& p, double mchi, double epsilon) {
    return epsilon * epsilon * kAlpha * p.reference_br * phaseSpace(p, mchi);
}

double branchD(TTree* t, const char* preferred, const char* fallback, double def=0.0) {
    double v = def;
    if (t->GetBranch(preferred)) t->SetBranchAddress(preferred, &v);
    else if (fallback && t->GetBranch(fallback)) t->SetBranchAddress(fallback, &v);
    return v;
}
}

void plot_normalized_mcp_yield(const char* inputFile="outputs/aggregate_summary_4h.root",
                               double pot=1.5e19,
                               double epsilon=1e-2,
                               bool draw_old_proxy=false) {
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << inputFile << "\n"; return; }
    TTree* t = (TTree*)f->Get("mcp_summary");
    if (!t) { std::cerr << "Need mcp_summary in " << inputFile << "\n"; return; }

    double mass=0, sigma_sum=0, sigma_mean=0;
    int pdg=0, mode=0, nfiles=1;
    long long events=0, total=0, accepted=0;
    t->SetBranchAddress(t->GetBranch("mcp_mass_GeV") ? "mcp_mass_GeV" : "mcp_mass", &mass);
    t->SetBranchAddress(t->GetBranch("emitter_pdg") ? "emitter_pdg" : "parent_pdg", &pdg);
    if (t->GetBranch("production_mode")) t->SetBranchAddress("production_mode", &mode);
    t->SetBranchAddress("n_events_generated", &events);
    t->SetBranchAddress("n_mcp_total", &total);
    t->SetBranchAddress("n_mcp_accepted", &accepted);
    if (t->GetBranch("sigma_gen_mb_sum")) t->SetBranchAddress("sigma_gen_mb_sum", &sigma_sum);
    if (t->GetBranch("sigma_gen_mb_mean")) t->SetBranchAddress("sigma_gen_mb_mean", &sigma_mean);
    if (t->GetBranch("n_input_files")) t->SetBranchAddress("n_input_files", &nfiles);

    std::vector<Row> rows;
    double softSigmaSum=0, charmSigmaSum=0; int nSoft=0, nCharm=0;
    for (Long64_t i=0; i<t->GetEntries(); ++i) {
        sigma_sum=0; sigma_mean=0; nfiles=1; t->GetEntry(i);
        double mean = sigma_mean > 0 ? sigma_mean : (nfiles > 0 ? sigma_sum / nfiles : sigma_sum);
        rows.push_back({mass,pdg,mode,events,total,accepted,sigma_sum,mean,nfiles});
        if (mean > 0 && mode == kLightMesons) { softSigmaSum += mean; ++nSoft; }
        if (mean > 0 && mode == kCharmonium) { charmSigmaSum += mean; ++nCharm; }
    }
    const double sigmaSoft = nSoft ? softSigmaSum / nSoft : 1.0;
    const double sigmaCharm = nCharm ? charmSigmaSum / nCharm : 0.0;
    const double charmScale = (sigmaSoft > 0.0) ? sigmaCharm / sigmaSoft : 0.0;
    std::cout << "N_POT=" << pot << " epsilon=" << epsilon
              << " sigma_softqcd_mb_mean=" << sigmaSoft
              << " sigma_charmonium_mb_mean=" << sigmaCharm
              << " charmonium_process_scale=" << charmScale << "\n";

    auto parents = parentTable();
    std::map<int, std::vector<Point>> byParent;
    std::map<double, double> totalYieldByMass;
    for (const Row& r : rows) {
        if (!parents.count(r.pdg)) continue;
        const ParentPhysics& p = parents[r.pdg];
        const double acc = r.total ? double(r.accepted) / double(r.total) : 0.0;
        const double br = brExotic(p, r.mass, epsilon);
        const double processScale = (r.mode == kCharmonium) ? charmScale : 1.0;
        const double yield = (r.events > 0) ? pot * processScale * (double(r.accepted) / double(r.events)) * br : 0.0;
        const double oldProxy = pot * (r.events ? double(r.accepted) / double(r.events) : 0.0) * br;
        byParent[r.pdg].push_back({r.mass, acc, yield, oldProxy});
        totalYieldByMass[r.mass] += yield;
    }

    TCanvas* c1 = new TCanvas("c_acceptance_emitters_norm", "Acceptance by emitter", 900, 700);
    c1->SetLogx(); TLegend* l1 = new TLegend(0.62,0.55,0.88,0.88); bool first=true; int color=1;
    for (auto& kv : byParent) {
        std::sort(kv.second.begin(), kv.second.end(), [](const Point&a,const Point&b){return a.mass<b.mass;});
        std::vector<double> x,y; for (auto&p:kv.second){x.push_back(p.mass); y.push_back(p.acceptance);} 
        TGraph* g = new TGraph(x.size(), x.data(), y.data()); g->SetMarkerStyle(20); g->SetMarkerColor(color); g->SetLineColor(color++);
        g->SetTitle("Acceptance vs MCP mass;MCP mass [GeV];n_{MCP accepted}/n_{MCP total}"); g->Draw(first?"APL":"PL"); first=false; l1->AddEntry(g, parents[kv.first].name.c_str(), "lp");
    }
    l1->Draw();

    TCanvas* c2 = new TCanvas("c_normalized_yield_emitters", "Accepted MCP yield by emitter", 900, 700);
    c2->SetLogx(); c2->SetLogy(); TLegend* l2 = new TLegend(0.50,0.50,0.88,0.88); first=true; color=1;
    for (auto& kv : byParent) {
        std::vector<double> x,y; for (auto&p:kv.second){x.push_back(p.mass); y.push_back(p.yield);} 
        TGraph* g = new TGraph(x.size(), x.data(), y.data()); g->SetMarkerStyle(21); g->SetMarkerColor(color); g->SetLineColor(color++);
        g->SetTitle(Form("Accepted MCPs at detector;MCP mass [GeV];Accepted MCPs, N_{POT}=%.2g, #epsilon=%.2g (placeholder branching normalization)", pot, epsilon));
        g->Draw(first?"APL":"PL"); first=false; l2->AddEntry(g, parents[kv.first].name.c_str(), "lp");
    }
    l2->Draw();

    TCanvas* c3 = new TCanvas("c_normalized_yield_total", "Total accepted MCP yield", 900, 700);
    c3->SetLogx(); c3->SetLogy(); std::vector<double> xt,yt; for (auto& kv: totalYieldByMass){xt.push_back(kv.first); yt.push_back(kv.second);} 
    TGraph* gt = new TGraph(xt.size(), xt.data(), yt.data()); gt->SetMarkerStyle(22); gt->SetLineWidth(2);
    gt->SetTitle(Form("Total accepted MCPs at detector;MCP mass [GeV];Accepted MCPs, N_{POT}=%.2g, #epsilon=%.2g (not final sensitivity)", pot, epsilon)); gt->Draw("APL");

    if (draw_old_proxy) {
        TCanvas* c4 = new TCanvas("c_old_proxy_compare", "Old proxy comparison", 900, 700);
        c4->SetLogx(); c4->SetLogy(); first=true; color=1; TLegend* l4 = new TLegend(0.50,0.50,0.88,0.88);
        for (auto& kv : byParent) { std::vector<double>x,y; for(auto&p:kv.second){x.push_back(p.mass); y.push_back(p.old_proxy);} TGraph* g=new TGraph(x.size(),x.data(),y.data()); g->SetMarkerColor(color); g->SetLineColor(color++); g->SetMarkerStyle(24); g->SetTitle("Old proxy yield;MCP mass [GeV];old proxy"); g->Draw(first?"APL":"PL"); first=false; l4->AddEntry(g, parents[kv.first].name.c_str(), "lp"); }
        l4->Draw();
    }

    std::cout << "Note: branching constants are placeholder/model inputs; do not interpret as final sensitivity.\n";
}
