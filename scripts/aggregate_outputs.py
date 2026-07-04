#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path
try:
    import ROOT
except Exception as e:
    raise SystemExit(f"PyROOT is required: {e}")
ap=argparse.ArgumentParser(); ap.add_argument('inputs',nargs='*',default=['outputs/raw']); ap.add_argument('--csv',default='outputs/aggregate_summary.csv'); ap.add_argument('--root',default='outputs/aggregate_summary.root'); ap.add_argument('--merge-spectra',action='store_true')
a=ap.parse_args(); files=[]
for x in a.inputs:
    p=Path(x); files += ([p] if p.is_file() else list(p.rglob('*.root')))
agg={}
for fn in files:
    f=ROOT.TFile.Open(str(fn));
    if not f or f.IsZombie(): continue
    t=f.Get('mcp_summary')
    if not t: continue
    for row in t:
        key=(round(float(row.mcp_mass_GeV),12), int(row.emitter_pdg), int(row.production_mode), int(row.geometry_id))
        d=agg.setdefault(key, dict(n_events_generated=0,n_emitter_total=0,n_mcp_total=0,n_mcp_accepted=0,n_mcp_wrong_mother=0,n_mcp_no_mother=0))
        for k in d: d[k]+=int(getattr(row,k))
    f.Close()
Path(a.csv).parent.mkdir(parents=True,exist_ok=True)
with open(a.csv,'w',newline='') as f:
    fields=['mcp_mass_GeV','emitter_pdg','production_mode','geometry_id','n_events_generated','n_emitter_total','n_mcp_total','n_mcp_accepted','n_mcp_wrong_mother','n_mcp_no_mother','acceptance_fraction','acceptance_uncertainty_binomial']
    w=csv.DictWriter(f,fields); w.writeheader()
    for key,d in sorted(agg.items()):
        n=d['n_mcp_total']; p=d['n_mcp_accepted']/n if n else 0; unc=math.sqrt(max(0,p*(1-p)/n)) if n else 0
        w.writerow(dict(zip(fields[:4],key))|d|{'acceptance_fraction':p,'acceptance_uncertainty_binomial':unc})

# Also write a compact ROOT aggregate tree for the plotting macro.
rout=ROOT.TFile(a.root,'RECREATE')
tout=ROOT.TTree('mcp_summary','Aggregated MCP summary')
from array import array
m=array('d',[0.0]); p=array('i',[0]); pm=array('i',[0]); g=array('i',[0]); ev=array('l',[0]); em=array('l',[0]); tot=array('l',[0]); ac=array('l',[0]); wm=array('l',[0]); nm=array('l',[0]); af=array('d',[0.0]); au=array('d',[0.0])
for name,var,leaf in [('mcp_mass_GeV',m,'mcp_mass_GeV/D'),('emitter_pdg',p,'emitter_pdg/I'),('production_mode',pm,'production_mode/I'),('geometry_id',g,'geometry_id/I'),('n_events_generated',ev,'n_events_generated/L'),('n_emitter_total',em,'n_emitter_total/L'),('n_mcp_total',tot,'n_mcp_total/L'),('n_mcp_accepted',ac,'n_mcp_accepted/L'),('n_mcp_wrong_mother',wm,'n_mcp_wrong_mother/L'),('n_mcp_no_mother',nm,'n_mcp_no_mother/L'),('acceptance_fraction',af,'acceptance_fraction/D'),('acceptance_uncertainty_binomial',au,'acceptance_uncertainty_binomial/D')]: tout.Branch(name,var,leaf)
for key,d in sorted(agg.items()):
    m[0],p[0],pm[0],g[0]=key; ev[0]=d['n_events_generated']; em[0]=d['n_emitter_total']; tot[0]=d['n_mcp_total']; ac[0]=d['n_mcp_accepted']; wm[0]=d['n_mcp_wrong_mother']; nm[0]=d['n_mcp_no_mother']; af[0]=ac[0]/tot[0] if tot[0] else 0; au[0]=math.sqrt(max(0,af[0]*(1-af[0])/tot[0])) if tot[0] else 0; tout.Fill()
tout.Write(); rout.Close()
print(f"Wrote {len(agg)} aggregate rows to {a.csv} and {a.root}")

