#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path
EMITTERS={"pi0":(111,0.1349768,"light_mesons"),"eta":(221,0.5478620,"light_mesons"),"etap":(331,0.9577800,"light_mesons"),"rho0":(113,0.7752600,"light_mesons"),"omega":(223,0.7826500,"light_mesons"),"phi":(333,1.0194610,"light_mesons"),"jpsi":(443,3.0969000,"charmonium")}
def tag(x): return f"{x:.6f}".replace('.','p')
ap=argparse.ArgumentParser()
ap.add_argument('--masses',default='masses.txt'); ap.add_argument('--emitters',default='pi0,eta,etap,rho0,omega,phi,jpsi')
ap.add_argument('--output',default='job_list.csv'); ap.add_argument('--output-dir',default='outputs/raw')
ap.add_argument('--n-events',type=int,default=10000); ap.add_argument('--geometry',default='2x2'); ap.add_argument('--seed-base',type=int,default=1000); ap.add_argument('--shards',type=int,default=1); ap.add_argument('--include-closed',action='store_true')
a=ap.parse_args(); masses=[]
for line in Path(a.masses).read_text().splitlines():
    line=line.strip()
    if line and not line.startswith('#'): masses.append(float(line))
rows=[]; jid=0
for em in [x.strip() for x in a.emitters.split(',') if x.strip()]:
    pdg,mp,mode=EMITTERS[em]
    for m in masses:
        if not a.include_closed and not (mp>2*m): continue
        for shard in range(a.shards):
            out=Path(a.output_dir)/f"emitter_{em}"/f"mass_{tag(m)}"/f"job_{jid:06d}.root"
            rows.append(dict(job_id=jid,emitter_name=em,emitter_pdg=pdg,production_mode=mode,mcp_mass_GeV=m,seed=a.seed_base+jid,n_events=a.n_events,geometry=a.geometry,output_file=str(out)))
            jid+=1
with open(a.output,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['job_id','emitter_name','emitter_pdg','production_mode','mcp_mass_GeV','seed','n_events','geometry','output_file']); w.writeheader(); w.writerows(rows)
print(f"Wrote {len(rows)} jobs to {a.output}")
