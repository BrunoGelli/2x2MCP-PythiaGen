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
count_fields=['n_events_generated','n_emitter_record_entries','n_emitter_decayed_to_mcp','n_mcp_pairs','n_mcp_pairing_anomalies','n_emitter_total','n_mcp_total','n_mcp_accepted','n_mcp_wrong_mother','n_mcp_no_mother']
agg={}
for fn in files:
    f=ROOT.TFile.Open(str(fn))
    if not f or f.IsZombie(): continue
    t=f.Get('mcp_summary')
    if not t: continue
    for row in t:
        key=(round(float(row.mcp_mass_GeV),12), int(row.emitter_pdg), int(row.production_mode), int(row.geometry_id))
        d=agg.setdefault(key, {k:0 for k in count_fields} | {'sigma_gen_mb_sum':0.0,'sigma_err_mb_quadrature':0.0})
        for k in count_fields:
            if hasattr(row,k): d[k]+=int(getattr(row,k))
            elif k=='n_emitter_decayed_to_mcp' and hasattr(row,'n_emitter_total'): d[k]+=int(row.n_emitter_total)
            elif k=='n_emitter_total' and hasattr(row,'n_emitter_decayed_to_mcp'): d[k]+=int(row.n_emitter_decayed_to_mcp)
        if hasattr(row,'sigma_gen_mb'): d['sigma_gen_mb_sum']+=float(row.sigma_gen_mb)
        if hasattr(row,'sigma_err_mb'): d['sigma_err_mb_quadrature']+=float(row.sigma_err_mb)**2
    f.Close()
Path(a.csv).parent.mkdir(parents=True,exist_ok=True)
fields=['mcp_mass_GeV','emitter_pdg','production_mode','geometry_id']+count_fields+['acceptance_fraction','acceptance_uncertainty_binomial','sigma_gen_mb_sum','sigma_err_mb_quadrature','weight_per_event_mb_sum']
with open(a.csv,'w',newline='') as f:
    w=csv.DictWriter(f,fields); w.writeheader()
    for key,d in sorted(agg.items()):
        n=d['n_mcp_total']; p=d['n_mcp_accepted']/n if n else 0; unc=math.sqrt(max(0,p*(1-p)/n)) if n else 0
        d['sigma_err_mb_quadrature']=math.sqrt(d['sigma_err_mb_quadrature'])
        d['weight_per_event_mb_sum']=d['sigma_gen_mb_sum']/d['n_events_generated'] if d['n_events_generated'] else 0
        w.writerow(dict(zip(fields[:4],key))|d|{'acceptance_fraction':p,'acceptance_uncertainty_binomial':unc})
# Compact ROOT tree mirrors CSV. Cross-section fields are carried for diagnostics;
# do not interpret summed charmonium sigma as an absolute normalization without validation.
rout=ROOT.TFile(a.root,'RECREATE'); tout=ROOT.TTree('mcp_summary','Aggregated MCP summary')
from array import array
vals={}
for name in fields:
    if name in ('mcp_mass_GeV','acceptance_fraction','acceptance_uncertainty_binomial','sigma_gen_mb_sum','sigma_err_mb_quadrature','weight_per_event_mb_sum'):
        vals[name]=array('d',[0.0]); leaf=f'{name}/D'
    elif name in ('emitter_pdg','production_mode','geometry_id'):
        vals[name]=array('i',[0]); leaf=f'{name}/I'
    else:
        vals[name]=array('l',[0]); leaf=f'{name}/L'
    tout.Branch(name,vals[name],leaf)
for key,d in sorted(agg.items()):
    n=d['n_mcp_total']; p=d['n_mcp_accepted']/n if n else 0; unc=math.sqrt(max(0,p*(1-p)/n)) if n else 0
    out=dict(zip(fields[:4],key))|d|{'acceptance_fraction':p,'acceptance_uncertainty_binomial':unc}
    for name,val in out.items(): vals[name][0]=val
    tout.Fill()
tout.Write(); rout.Close()
print(f"Wrote {len(agg)} aggregate rows to {a.csv} and {a.root}")
