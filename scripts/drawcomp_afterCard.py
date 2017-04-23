# to compare histograms - max samples with different options
# python scripts/drawcomp_afterBDT.py -n -c -w 0 -b 20170116-181049 

import json
import os
from glob import glob
import importlib

# ROOT imports
from ROOT import TChain, TH1F, TFile, vector, TCanvas, gROOT

from Analysis.alp_analysis.samplelists import samlists
from Analysis.alp_analysis.alpSamples import samples
from Analysis.alp_plots.histOpt import hist_opt
import Analysis.alp_plots.UtilsDraw as UtilsDraw

TH1F.AddDirectory(0)
gROOT.SetBatch(True)

# parsing parameters
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-w", "--whichPlots", help="which plots to be produced", type=int, default='-1')
parser.add_argument("-b", "--bdt"       , help="bdt version, equal to input file name", default="")
parser.add_argument("-o", "--oDir"     , help="output directory"        , default="plots_moriond")
parser.add_argument("-r", "--clrebin", help="to rebin (classifier output)"    , type=int, default=-1)
parser.add_argument("--res", dest="plotResidual", help="to plot residuals (2==fit)" , type=int, default=0)
parser.add_argument("-n", "--doNorm"    , help="do not normalize"       , action='store_false')
parser.add_argument("-c", "--customCol" , help="do not use custom colors"       , action='store_false')
parser.add_argument('--ra', help="run benchmarks (to change naming convention to script)", default=-1, type=int)
parser.set_defaults(doNorm=True, customCol=True, plotResidual=False)
args = parser.parse_args()

ra = args.ra
#if ra == -1 :  iDir       = '../hh2bbbb_limit/' #'/lustre/cmswork/hh/alp_afterMVA/'
#else :  
iDir       = '/afs/cern.ch/work/a/acarvalh/codeCMSHHH4b/CMSSW_8_0_25/src/Analysis/hh2bbbb_limit/datacards/' #'/lustre/cmswork/hh/alp_afterMVA/'
fileto = ["SM", "BMbox","BM2","BM3","BM4","BM5","BM6","BM7","BM8","BM9","BM10","BM11","BM12","BM13"]
filename1 = iDir+"hists_VOBM_plots_"
filename1end="_reweighted.root"
filename2 = iDir+"hists_VOBM_plots_ggh_hh_bbbb_"
filename2end="_plain.root"

# exe parameters
histList   = [ 'h_jet0_pt', 'h_jet1_pt', 'h_jet2_pt', 'h_jet3_pt', 'h_jets_ht', 
              'h_jet2_csv','h_jet3_csv',
              'h_jet2_cmva','h_jet3_cmva',
              'h_H0_mass','h_H0_csthst0_a','h_H0_dr', #'h_H0_pt','h_H0_eta','h_H0_csthst1_a','h_H0_deta_a','h_H0_dphi_a',
              'h_H1_mass','h_H1_csthst2_a','h_H1_dr', #'h_H1_pt','h_H1_eta','h_H1_csthst2_a','h_H1_csthst3_a','h_H1_dr','h_H1_deta_a','h_H1_dphi_a',
              'h_H0H1_mass', 'h_H0H1_csthst0_a', 'h_H0H1_dr', #'h_H0H1_pt','h_H0H1_eta','h_H0H1_csthst1_a',,'h_H0H1_deta_a','h_H0H1_dphi_a',
              'h_jet0_eta', 'h_jet1_eta', 'h_jet2_eta', 'h_jet3_eta',
              'h_X_mass', 'h_jets_ht_r',               
             ]
histList2  = ["DiJets[0].mass()-DiJets[1].mass()", "CSV_Jet2-CSV_Jet3", "CMVA_Jet2-CMVA_Jet3",] #2D histos,
intLumi_fb = 1. # plots normalized to this

#h_X_mass_signal
#h_X_mass_signal

for ii in range(0,len(fileto)-1) :
    samples = ['',''] #['sig', 'sig']
    fractions = ['test','']
    regions = ['','']
    legList = ["ggHH4b "+fileto[ii]+" V0","Pangea"]
    colorList = [632, 630]
    dofill = [True,True]
    oname = 'comp_'+fileto[ii]+'_afterCard'
    if not args.customCol: colors = []
    else: colors = colorList
    plotDir = ''
    oDir = args.oDir
    oDir += "/"
    if not os.path.exists(oDir): os.mkdir(oDir)
    oDir += "/"+oname
    if args.doNorm: oDir = oDir+fileto[ii]+"_norm/"
    else: oDir = oDir+"/"
    if not os.path.exists(oDir): os.mkdir(oDir)
    #----------------------------------
    for h in histList:
      hOpt = hist_opt[h]
      hs1 = UtilsDraw.getHistos_bdt(h+"_signal", filename2+fileto[ii]+filename2end, plotDir)
      hs2 = UtilsDraw.getHistos_bdt(h+"_signal", filename1+fileto[ii]+filename1end, plotDir)
      if hs1 and hs2: UtilsDraw.drawH1(hs1, legList[0], hs2, legList[1], hOpt, args.plotResidual, args.doNorm, oDir, colors, dofill, args.clrebin)

