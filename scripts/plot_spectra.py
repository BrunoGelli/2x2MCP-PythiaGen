#!/usr/bin/env python3
import argparse
from pathlib import Path
try:
    import ROOT
except Exception as e:
    raise SystemExit(f"PyROOT is required: {e}")

def add_files(chain, paths):
    n=0
    for x in paths:
        p=Path(x)
        files=[p] if p.is_file() else list(p.rglob('*.root'))
        for f in files:
            chain.Add(str(f)); n+=1
    return n

def save(c, outdir, name):
    c.SaveAs(str(Path(outdir)/f'{name}.png'))

def main():
    ap=argparse.ArgumentParser(description='Plot MCP spectra from raw ROOT mcp_spectra trees.')
    ap.add_argument('--spectra-dir', action='append', required=True)
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--accepted-only', action='store_true')
    args=ap.parse_args()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    ch=ROOT.TChain('mcp_spectra')
    nfiles=add_files(ch,args.spectra_dir)
    if ch.GetEntries()==0:
        raise SystemExit(f'No mcp_spectra entries found in {nfiles} files')
    cut='passed_geometry==1' if args.accepted_only else ''
    ROOT.gStyle.SetOptStat(0)
    c=ROOT.TCanvas('c_theta_xy','theta_x vs theta_y',900,750)
    ch.Draw('theta_y_rad:theta_x_rad>>h_theta_xy(160,-0.2,0.2,160,-0.2,0.2)',cut,'COLZ')
    save(c,args.output_dir,'theta_x_vs_theta_y')
    c=ROOT.TCanvas('c_theta','theta distribution',900,750); c.SetLogy()
    ch.Draw('sqrt(theta_x_rad*theta_x_rad+theta_y_rad*theta_y_rad)>>h_theta(160,0,0.3)',cut,'HIST')
    save(c,args.output_dir,'theta_distribution')
    c=ROOT.TCanvas('c_pz_theta','pz vs theta',900,750); c.SetLogz()
    ch.Draw('pz_GeV:sqrt(theta_x_rad*theta_x_rad+theta_y_rad*theta_y_rad)>>h_pz_theta(160,0,0.3,160,0,120)',cut,'COLZ')
    save(c,args.output_dir,'pz_vs_theta')
    c=ROOT.TCanvas('c_energy','E distribution',900,750); c.SetLogy()
    ch.Draw('E_GeV>>h_E(160,0,120)',cut,'HIST')
    save(c,args.output_dir,'energy_distribution')
    c=ROOT.TCanvas('c_xy','x/y at detector',900,750); c.SetLogz()
    ch.Draw('y_at_detector_m:x_at_detector_m>>h_xy(200,-5,5,200,-5,5)',cut,'COLZ')
    boxes=[(-0.65,-0.7,-0.05,0.7),(0.05,-0.7,0.65,0.7)]
    drawn_boxes=[]
    for x1,y1,x2,y2 in boxes:
        b=ROOT.TBox(x1,y1,x2,y2); b.SetFillStyle(0); b.SetLineColor(ROOT.kRed); b.SetLineWidth(3); b.Draw('same'); drawn_boxes.append(b)
    save(c,args.output_dir,'x_y_at_detector')
    c=ROOT.TCanvas('c_xy_acc','accepted x/y hit map',900,750)
    ch.Draw('y_at_detector_m:x_at_detector_m>>h_xy_acc(120,-0.8,0.8,120,-0.8,0.8)','passed_geometry==1','COLZ')
    drawn_boxes_acc=[]
    for x1,y1,x2,y2 in boxes:
        b=ROOT.TBox(x1,y1,x2,y2); b.SetFillStyle(0); b.SetLineColor(ROOT.kRed); b.SetLineWidth(3); b.Draw('same'); drawn_boxes_acc.append(b)
    save(c,args.output_dir,'accepted_x_y_hit_map')
    print(f'Wrote spectra plots to {args.output_dir} from {int(ch.GetEntries())} entries in {nfiles} files')

if __name__ == '__main__':
    main()
