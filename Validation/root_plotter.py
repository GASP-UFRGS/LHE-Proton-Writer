from __future__ import division
from subprocess import call
from math import *
from ROOT import *
import numpy as np

rootinput = np.loadtxt('root_input.txt', unpack=True, dtype=str, delimiter='=')

# CROSS SECTION(S) (pb):
xsec    = eval(rootinput[1][0])

# Plotting options:
JOB     = eval(rootinput[1][1])
PDF     = eval(rootinput[1][2])
scale   = eval(rootinput[1][3]) 
cuts    = eval(rootinput[1][4])
setLog  = eval(rootinput[1][5])
filled  = eval(rootinput[1][6])
stacked = eval(rootinput[1][7])
data    = eval(rootinput[1][8])
LUMI    = eval(rootinput[1][9])

# EVENT SAMPLE INPUT:
Nevt    = eval(rootinput[1][10])
Nmax    = eval(rootinput[1][11])
EVTINPUT= str(int(Nmax/1000))+"k";
SQRTS   = eval(rootinput[1][12])
lumi    = float(rootinput[1][14])


# LABELS:
LABEL = "FULL_inner.final.madgraph";
if cuts: LABEL+=".cuts";
if scale: LABEL+=".scaled";
if setLog: LABEL+=".log";
if filled: LABEL+=".filled";
if stacked: LABEL+=".stacked";
if data: LABEL+=".data";

# IMAGE FORMATS TO BE CREATED:
FILE_TYPES = [LABEL+".png"];

# CREATE INDIVIDUAL DIRS FOR IMAGE TYPES:
for l in range(len(FILE_TYPES)):
	call(["mkdir","-p",FILE_TYPES[l]]);

FILEROOT = TFile.Open("histos"+LABEL+".root", "READ")

keys = FILEROOT.GetListOfKeys()

# ARRAYS FOR EACH TYPE OF DISTRIBUTIONS:
#
# 1D:
protpz       = []
proten       = []
protxi       = []
protpt       = []
proteta      = []
mpp          = []
mupz         = []
muen         = []
mupt         = []
ivmmu        = []
mueta        = []
phopz        = []
phopt        = []
phoen        = []
phoivm       = []
phopsrap2    = [] 
phopsrap1    = []
phoY         = []

# 2D:
DDmppmmumu   = []
DDxipximu    = []

# READING FROM ROOT FILE AND ALLOCATING DATA IN LISTS
for k in keys:
    name = k.GetName().split('_')[1]
    eval(name).append(k.ReadObj())

# THE ARRAYS STORE THE LABELS FOR AXIS AND UNITS:
# 1D
histoslog        = [protpz,proten,protxi,protpt,proteta,mpp,mupz,muen,mupt,ivmmu,mueta,phopz,phopt,phoen,phoivm,phopsrap2,phopsrap1,phoY]
histoslog_label  = ["protpz","proten",'protxi','protpt','proteta','mpp',"mupz","muen","mupt",'ivm_mu','mueta','phopz','phopt','phoen','phoivm','phopsrap2','phopsrap1','phoY']
histoslog_axis   = ["p_{z}(p)","E(p)",'#chi(p)','p_{T}(p)','#eta(p^{+}p^{-})','M(p^{+}p^{-})',"p_{z}(#mu)","E(#mu)","p_{T}(#mu)",'M(#mu^{+}#mu^{-})','#eta(#mu^{+}#mu^{-})','p_{z}(#gamma#gamma)','p_{T}(#gamma#gamma)','E(#gamma#gamma)','M(#gamma#gamma)','#eta(#gamma#gamma)','#eta(#gamma)','Y(#gamma#gamma)']
histoslog_varx   = ["(GeV)","(GeV)",'','(GeV)','','(GeV)',"(GeV)","(GeV)","(GeV)",'(GeV)','','(GeV)','(GeV)','(GeV)','(GeV)','','','']

# 2D
DDlog         = [DDmppmmumu,DDxipximu] 
DDlog_label   = ["2Dmppmmumu",'2Dxipximu'];
DDlog_xaxis   = ["M(p^{+}p^{-})",'#xi(p^{+})']
DDlog_yaxis   = ["M(#mu^{+}#mu^{-})",'#xi(#mu^{+}#mu^{-})']
DDlog_varx    = ["(GeV)",'']
DDlog_vary    = ["(GeV)",'']

# FETCHING NUMBER OF FILES
NUMFILES = len(histoslog[0])

# Starting Drawing step:

# SETTING THE NUMBER OF DIGITS ON AXIS
TGaxis.SetMaxDigits(2)

# Defining the top label in the plots:
plotlabel = TPaveText(0.50,0.91,0.84,0.95,"NDC");
plotlabel.SetTextAlign(33);
plotlabel.SetTextColor(1);
plotlabel.SetFillColor(0);
plotlabel.SetBorderSize(0);
plotlabel.SetTextSize(0.035);
plotlabel.SetTextFont(42);
plotlabel.AddText("MadGraphv5 #bullet #sqrt{s}=13 TeV #bullet "+EVTINPUT+" evt");

# Legend:
leg = TLegend(0.55,0.72,0.75,0.87);
leg.SetTextSize(0.035);
leg.SetFillColor(0);
leg.SetBorderSize(0);

# Setting pads:
gStyle.SetOptStat(0);
gStyle.SetPadTickY(1);
gStyle.SetPadTickX(1);
gStyle.SetOptTitle(0);
gStyle.SetLegendBorderSize(0);

# Canvas
canvas = TCanvas("plots","Plots",0,0,700,860);

for i in range(len(histoslog)):
    globals()["hs_histoslog"+str(i)] = THStack("hs","");

# Starting loop over histograms in the arrays for each set:

# 1: 1D log-scaled plots:
canvas.SetLeftMargin(0.2);
canvas.SetBottomMargin(0.11);
canvas.SetRightMargin(0.18);
if setLog: gPad.SetLogy(1);
else: gPad.SetLogy(0);
legs=0;
for l in range(len(histoslog)):
    for m in range(NUMFILES):
            if lumi:
                histoslog[l][m].Scale(xsec[m]/Nevt*lumi*histoslog[l][m].GetBinWidth(1))
            if scale:
                    histoslog[l][m].Scale(xsec[m]/Nevt*histoslog[l][m].GetBinWidth(1));
            histoslog[l][m].SetLineColor(m+1);
            if (m == 4): histoslog[l][m].SetLineColor(m+2);
            if filled:
                    histoslog[l][m].SetFillColor(m+1);
                    if (m == 4): histoslog[l][m].SetFillColor(m+2);
            histoslog[l][m].SetLineWidth(3);
            histoslog[l][m].SetLineStyle(1);
            globals()["hs_histoslog"+str(l)].Add(histoslog[l][m]);
            leg.AddEntry(histoslog[l][m]," "+PDF[m],"f");
            if data:
                if m == 0:
                        datapoints = histoslog[l][m].Clone();
                        dataonly = histoslog[l][m].Clone();
                else:
                        datapoints.Add(histoslog[l][m]);
                        dataonly.Add(histoslog[l][m]);
                        datapoints.SetFillStyle(0);
                        datapoints.SetLineWidth(0);
                        datapoints.SetLineStyle(0);
                        datapoints.SetMarkerStyle(20);
    if stacked:
            globals()["hs_histoslog"+str(l)].Draw("");
    else:
            globals()["hs_histoslog"+str(l)].Draw("nostack hist");
    if scale:
            globals()["hs_histoslog"+str(l)].GetYaxis().SetTitle("d#sigma/d"+str(histoslog_axis[l])+str(histoslog_varx[l])+" (pb)");
    else:
            globals()["hs_histoslog"+str(l)].GetYaxis().SetTitle("Events");
    if LUMI:
            globals()["hs_histoslog"+str(l)].GetYaxis().SetTitle("Yield")
    globals()["hs_histoslog"+str(l)].GetXaxis().SetTitle(str(histoslog_axis[l])+" "+str(histoslog_varx[l]));
    globals()["hs_histoslog"+str(l)].GetXaxis().SetTitleFont(42);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetTitleFont(42);
    globals()["hs_histoslog"+str(l)].GetXaxis().SetTitleSize(0.05);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetTitleSize(0.05);
    globals()["hs_histoslog"+str(l)].GetXaxis().SetLabelFont(42);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetLabelFont(42);
    globals()["hs_histoslog"+str(l)].GetXaxis().SetTitleOffset(1.);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetTitleOffset(1.6);
    globals()["hs_histoslog"+str(l)].GetXaxis().SetLabelSize(0.04);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetLabelSize(0.04);
    globals()["hs_histoslog"+str(l)].GetXaxis().SetDecimals(True);
    globals()["hs_histoslog"+str(l)].GetYaxis().SetNdivisions(5, optim=kTRUE)
    globals()["hs_histoslog"+str(l)].GetXaxis().SetNdivisions(10, optim=kTRUE)
    if data:
            datapoints.Draw("E2,SAME");
            leg.AddEntry(datapoints,"data","p");
    leg.Draw("SAME");
    plotlabel.Draw("SAME");
    for k in range(len(FILE_TYPES)):
        canvas.Print(FILE_TYPES[k]+"/"+JOB+"_"+EVTINPUT+"evt_"+histoslog_label[l]+"."+FILE_TYPES[k]);
    leg.Clear();
    if data:
        dataonly.SetLineStyle(2);
        dataonly.SetFillColor(0);
        dataonly.SaveAs(FILE_TYPES[k]+"/"+JOB+"_"+EVTINPUT+"evt_"+histoslog_label[l]+".root");
# END loop over plots in log scale

# 2: 2D plot in log-scale:
canvas.SetLeftMargin(0.15);
canvas.SetBottomMargin(0.11);
canvas.SetRightMargin(0.23);
canvas.SetFrameFillColor(887);
gPad.SetLogy(0);
for l in range(len(DDlog)):
    # printar os nomes
    for m in range(NUMFILES):
        if (scale):
            DDlog[l][m].Scale(xsec[m]/Nevt*DDlog[l][m].GetXaxis().GetBinWidth(1));
            DDlog[l][m].Draw("COLZ");
            leg.AddEntry(DDlog[l][m]," "+PDF[m],"f");
        plotlabel.Draw("SAME");
        DDlog[l][m].Draw("COLZ")
        DDlog[l][m].GetXaxis().SetTitleOffset(1.);
        DDlog[l][m].GetYaxis().SetTitleOffset(1.2);
        DDlog[l][m].GetZaxis().SetTitleOffset(1.55);
        DDlog[l][m].GetXaxis().SetTitleFont(42);
        DDlog[l][m].GetYaxis().SetTitleFont(42);
        DDlog[l][m].GetXaxis().SetLabelFont(42);
        DDlog[l][m].GetYaxis().SetLabelFont(42);
        DDlog[l][m].GetXaxis().SetTitleSize(0.05);
        DDlog[l][m].GetYaxis().SetTitleSize(0.05);
        DDlog[l][m].GetZaxis().SetTitleSize(0.05);
        DDlog[l][m].GetXaxis().SetLabelSize(0.04);
        DDlog[l][m].GetYaxis().SetLabelSize(0.04);
        DDlog[l][m].GetZaxis().SetLabelSize(0.04);
        if (scale):
            DDlog[l][m].GetZaxis().SetTitle("d#sigma/d"+str(DDlog_xaxis[l])+"d"+str(DDlog_yaxis[l])+str(DDlog_varz[l]));
        else:
            DDlog[l][m].GetZaxis().SetTitle("Events");
        DDlog[l][m].GetYaxis().SetTitle(DDlog_yaxis[l]+" "+DDlog_vary[l]+"");
        DDlog[l][m].GetXaxis().SetTitle(DDlog_xaxis[l]+" "+DDlog_varx[l]+"");
        DDlog[l][m].GetXaxis().SetDecimals(True);
        for k in range(len(FILE_TYPES)):
            canvas.Print(FILE_TYPES[k]+"/"+JOB+"_"+EVTINPUT+"evt_"+DDlog_label[l]+'_'+PDF[m]+"."+FILE_TYPES[k]);
        leg.Clear();
# END loop over DDD Lego plots in log scale --

#####################################################################
#
# C'ESTI FINI
#
#####################################################################
