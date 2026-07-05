#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path
try:
    import ROOT
except Exception as e:
    raise SystemExit(f"PyROOT is required: {e}")

ALPHA = 1.0 / 137.0
PARENTS = {
    111: ('pi0', 0, 0.1349768, 0.98823),
    221: ('eta', 0, 0.5478620, 0.39410),
    331: ('etap', 0, 0.9577800, 0.0220),
    113: ('rho0', 1, 0.7752600, 4.72e-5),
    223: ('omega', 1, 0.7826500, 7.36e-5),
    333: ('phi', 1, 1.0194610, 2.973e-4),
    443: ('jpsi', 1, 3.0969000, 5.971e-2),
}
MODE_NAME = {0: 'light_mesons', 1: 'charmonium', 2: 'dy_reserved'}

def br_exotic(pdg, mchi, eps):
    name, typ, mp, br = PARENTS[pdg]
    if mp <= 0 or 2*mchi >= mp:
        return 0.0
    r = mchi*mchi/(mp*mp)
    if typ == 0:
        ps = max(0.0, 1.0 - 4.0*r)**3
    else:
        beta = math.sqrt(max(0.0, 1.0 - 4.0*r))
        ps = beta * (1.0 + 2.0*r)
    return eps*eps*ALPHA*br*ps

def read_summary(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (round(float(r['mcp_mass_GeV']), 12), int(r['emitter_pdg']), int(r['production_mode']), int(r['geometry_id']))
            rows[key] = r
    soft = [float(r.get('sigma_gen_mb_mean') or 0) for r in rows.values() if int(r['production_mode']) == 0 and float(r.get('sigma_gen_mb_mean') or 0) > 0]
    charm = [float(r.get('sigma_gen_mb_mean') or 0) for r in rows.values() if int(r['production_mode']) == 1 and float(r.get('sigma_gen_mb_mean') or 0) > 0]
    sigma_soft = sum(soft)/len(soft) if soft else 1.0
    sigma_charm = sum(charm)/len(charm) if charm else 0.0
    return rows, sigma_soft, sigma_charm

def root_files(paths):
    out=[]
    for x in paths:
        p=Path(x)
        out += ([p] if p.is_file() else list(p.rglob('*.root')))
    return out

def main():
    ap=argparse.ArgumentParser(description='Export weighted MCP spectra rows for toy detector MC.')
    ap.add_argument('--summary', required=True)
    ap.add_argument('--spectra-dir', action='append', required=True)
    ap.add_argument('--pot', type=float, default=1.5e19)
    ap.add_argument('--epsilon', type=float, default=1e-2)
    ap.add_argument('--output-csv', required=True)
    ap.add_argument('--accepted-only', action='store_true')
    ap.add_argument('--spectra-prescale-default', type=int, default=1)
    args=ap.parse_args()
    summary, sigma_soft, sigma_charm = read_summary(args.summary)
    charm_scale = sigma_charm/sigma_soft if sigma_soft > 0 else 0.0
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    fields=['mcp_mass_GeV','emitter_pdg','emitter_name','production_mode','E_GeV','px_GeV','py_GeV','pz_GeV','p_GeV','theta_x_rad','theta_y_rad','theta_rad','x_at_detector_m','y_at_detector_m','passed_geometry','event_weight','weight_per_POT','N_POT','epsilon','process_scale','br_exotic','spectra_prescale']
    sums={}
    nout=0
    with open(args.output_csv,'w',newline='') as fout:
        w=csv.DictWriter(fout,fields); w.writeheader()
        for fn in root_files(args.spectra_dir):
            f=ROOT.TFile.Open(str(fn))
            if not f or f.IsZombie(): continue
            t=f.Get('mcp_spectra')
            if not t: f.Close(); continue
            for row in t:
                if args.accepted_only and not int(row.passed_geometry):
                    continue
                key=(round(float(row.mcp_mass_GeV),12), int(row.emitter_pdg), int(row.production_mode), int(row.geometry_id))
                if key not in summary:
                    print(f'WARNING: no summary row for spectra key {key} in {fn}')
                    continue
                sr=summary[key]
                nev=int(sr['n_events_generated'])
                if nev <= 0: continue
                pdg=int(row.emitter_pdg); mode=int(row.production_mode)
                process_scale = charm_scale if mode == 1 else 1.0
                br=br_exotic(pdg, float(row.mcp_mass_GeV), args.epsilon)
                prescale=int(sr.get('spectra_prescale_max') or sr.get('spectra_prescale') or args.spectra_prescale_default)
                weight=args.pot*process_scale*br*prescale/nev
                out={
                    'mcp_mass_GeV':float(row.mcp_mass_GeV), 'emitter_pdg':pdg, 'emitter_name':PARENTS.get(pdg,('unknown',))[0],
                    'production_mode':MODE_NAME.get(mode,str(mode)), 'E_GeV':float(row.E_GeV), 'px_GeV':float(row.px_GeV),
                    'py_GeV':float(row.py_GeV), 'pz_GeV':float(row.pz_GeV), 'p_GeV':float(row.p_GeV),
                    'theta_x_rad':float(row.theta_x_rad), 'theta_y_rad':float(row.theta_y_rad), 'theta_rad':float(row.theta_rad),
                    'x_at_detector_m':float(row.x_at_detector_m), 'y_at_detector_m':float(row.y_at_detector_m),
                    'passed_geometry':int(row.passed_geometry), 'event_weight':weight, 'weight_per_POT':weight/args.pot,
                    'N_POT':args.pot, 'epsilon':args.epsilon, 'process_scale':process_scale, 'br_exotic':br,
                    'spectra_prescale':prescale,
                }
                w.writerow(out); nout+=1
                if int(row.passed_geometry): sums[key]=sums.get(key,0.0)+weight
            f.Close()
    print(f'Wrote {nout} weighted spectra rows to {args.output_csv}')
    print(f'sigma_softqcd_mb_mean={sigma_soft} sigma_charmonium_mb_mean={sigma_charm} charm_process_scale={charm_scale}')
    for key, sr in sorted(summary.items()):
        pdg=key[1]; nev=int(sr['n_events_generated'])
        if nev <= 0 or pdg not in PARENTS: continue
        process_scale = charm_scale if key[2] == 1 else 1.0
        expected=args.pot*process_scale*(int(sr['n_mcp_accepted'])/nev)*br_exotic(pdg,key[0],args.epsilon)
        got=sums.get(key,0.0)
        if expected or got:
            rel=(got-expected)/expected if expected else float('inf')
            print(f'validation key={key}: accepted_weight_sum={got:.6g} normalized_yield={expected:.6g} rel_diff={rel:.3g}')

if __name__ == '__main__':
    main()
