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
namespace { struct P{int pdg,type,mode; const char* name; double mass,br;}; std::map<int,P> tab(){return {{111,{111,0,0,"#pi^{0}",0.1349768,0.98823}},{221,{221,0,0,"#eta",0.547862,0.39410}},{331,{331,0,0,"#eta'",0.95778,0.022}},{113,{113,1,0,"#rho^{0}",0.77526,4.72e-5}},{223,{223,1,0,"#omega",0.78265,7.36e-5}},{333,{333,1,0,"#phi",1.019461,2.973e-4}},{443,{443,1,1,"J/#psi",3.0969,5.971e-2}}};} double ps(int type,double m,double mp){ if(mp<=0||2*m>=mp)return 0; double x=4*m*m/(mp*mp); if(type==0)return pow(std::max(0.0,1-x),3); return sqrt(std::max(0.0,1-x)); } struct A{long long ev=0,em=0,tot=0,acc=0;}; }
void plot_acceptance_vs_mass_multimeson(const char* inputFile="outputs/aggregate_summary.root", double pot=1.5e19, double epsilon=1e-2, bool sum_charmonium_with_light=false){
 gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); TFile*f=TFile::Open(inputFile); if(!f||f->IsZombie()){std::cerr<<"Cannot open "<<inputFile<<"\n";return;} TTree*t=(TTree*)f->Get("mcp_summary"); if(!t) t=(TTree*)f->Get("aggregate_summary"); if(!t){std::cerr<<"Need mcp_summary/aggregate_summary tree. CSV-only aggregate can be plotted after conversion or use raw ROOT files.\n";return;}
 double mass=0; int pdg=0,pmode=0; long long ev=0,em=0,tot=0,acc=0; t->SetBranchAddress(t->GetBranch("mcp_mass_GeV")?"mcp_mass_GeV":"mcp_mass",&mass); t->SetBranchAddress(t->GetBranch("emitter_pdg")?"emitter_pdg":"parent_pdg",&pdg); if(t->GetBranch("production_mode")) t->SetBranchAddress("production_mode",&pmode); t->SetBranchAddress("n_events_generated",&ev); t->SetBranchAddress(t->GetBranch("n_emitter_total")?"n_emitter_total":"n_parent_total",&em); t->SetBranchAddress("n_mcp_total",&tot); t->SetBranchAddress("n_mcp_accepted",&acc);
 std::map<std::pair<double,int>,A> agg; std::map<int,int> modes; for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); auto& a=agg[{mass,pdg}]; a.ev+=ev; a.em+=em; a.tot+=tot; a.acc+=acc; modes[pdg]=pmode;} auto info=tab(); std::map<int,std::vector<double>> xs,ys,fs; std::vector<double> xt,yt; for(auto& kv:agg){double m=kv.first.first; int p=kv.first.second; auto a=kv.second; if(!info.count(p))continue; double ac=a.tot?double(a.acc)/a.tot:0; double yld=a.ev?double(a.em)/a.ev:0; double w=epsilon*epsilon*info[p].br*ps(info[p].type,m,info[p].mass); double flux=pot*yld*2*ac*w; xs[p].push_back(m); ys[p].push_back(ac); fs[p].push_back(flux);} TCanvas*c1=new TCanvas("c_acceptance_emitters","Acceptance by emitter",900,700); c1->SetLogx(); TLegend*l=new TLegend(.62,.55,.88,.88); bool first=true; int col=1; for(auto& kv:xs){auto g=new TGraph(kv.second.size(),kv.second.data(),ys[kv.first].data()); g->SetMarkerStyle(20); g->SetMarkerColor(col); g->SetLineColor(col++); g->SetTitle("Acceptance vs MCP mass;MCP mass [GeV];accepted / produced MCP"); g->Draw(first?"APL":"PL"); first=false; l->AddEntry(g,info[kv.first].name,"lp");} l->Draw(); TCanvas*c2=new TCanvas("c_flux_proxy_emitters","Flux proxy by emitter",900,700); c2->SetLogx(); c2->SetLogy(); TLegend*l2=new TLegend(.62,.55,.88,.88); first=true; col=1; for(auto& kv:xs){auto g=new TGraph(kv.second.size(),kv.second.data(),fs[kv.first].data()); g->SetMarkerStyle(21); g->SetMarkerColor(col); g->SetLineColor(col++); g->SetTitle("Accepted MCP flux proxy;MCP mass [GeV];proxy (placeholder norm)"); g->Draw(first?"APL":"PL"); first=false; l2->AddEntry(g,info[kv.first].name,"lp");} l2->Draw(); if(!sum_charmonium_with_light) std::cout<<"Light-meson and charmonium modes are kept separate by default; pass true only for compatible toy sums.\n"; }
