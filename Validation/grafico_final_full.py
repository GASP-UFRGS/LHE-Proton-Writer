from __future__ import division
from subprocess import call
from math import *
from ROOT import *

#####################################################################
# GGS (CERN-CMS/UFRGS) ---
# the muons are collected considering the ID codes in the event
# sample produced with SuperCHICv2 in LHE format.
#####################################################################

#####################################################################
# USER INPUT:

# CROSS SECTION(S) (pb):
xsec    = [ 0.186818512E-03, 9.985100e-01, 0.13183148E+02, 1.5393433571E+00]; #FIXME
#xsec = [ 1. , 1. , 1. , 1. , .1 ];

# PDF "_"+LABEL FOR OUTPUT FILES:
JOB     = "histos";
PDF     = [ 'superchic', 'MadGraph', 'FPMC', 'LPAIR']; #FIXME
scale   = False; #bug, use False, carefull with Nevts
cuts    = True;
setLog  = False;
filled  = False;
stacked = False;
data    = False;

# KINEMATICAL CUTS: #FIXME
INVMCUTUPPER = 14000.0; # (NO CUT 9999.0 )
INVMCUTLOWER = 10.0; # (NO CUT 0.0)

PTPAIRCUTUPPER = 9999.0; # (NO CUT 0.0 )
PTPAIRCUTLOWER = 0.0; # (NO CUT 0.0)

ETACUT = 2.5  # NO CUT: INNER TURE, ETACUT = 20
ETAPAIRCUT = 5000.; # (NO CUT 100.)
INNER = True; # (TRUE: -x < y < +x ; FALSE: y < -x AND y > +x)

PTCUTUPPER = 9999.0; # (NO CUT 9999.0 )
PTCUTLOWER = 10.0; # (NO CUT 0.0)

# INPUT FILES:

#processo 3
FILES   = [
"samples/newevrectest.dat", 'samples/newunweighted_events.3.5.lhe', 'samples/fpmc_14tev.lhe', 'samples/Artur_gammagammamumu-lpair_elel_pt10_14tev_20k.lhe']
#FIXME

# EVENT SAMPLE INPUT:
Nevt    = 200000; #FIXME
Nmax    = 10000   # number of max events to obtain from samples
EVTINPUT= str(int(Nevt/1000))+"k";
SQRTS   = 14000         # in GeV

#####################################################################

# LABELS:
STRING	= "";
for m in range(len(PDF)):
	if (PDF[m]==PDF[-1]):
		STRING+=PDF[m]+"_";
	else:
		STRING+=PDF[m]+"-";

LABEL = "FULL_inner.final.madgraph";
if cuts: LABEL+=".cuts";
if scale: LABEL+=".scaled";
if setLog: LABEL+=".log";
if filled: LABEL+=".filled";
if stacked: LABEL+=".stacked";
if data: LABEL+=".data";

# IMAGE FORMATS TO BE CREATED:
FILE_TYPES = [LABEL+".png"];
print ("*****");
print ("Os arquivos gravados em %s" % (FILE_TYPES[0]));
print ("*****");
# SAVING HISTOS INTO ROOT FILE:
FILEROOT = TFile("histos"+LABEL+".root","RECREATE");

# CREATE INDIVIDUAL DIRS FOR IMAGE TYPES:
#print(len(FILE_TYPES))
for l in range(len(FILE_TYPES)):
	call(["mkdir","-p",FILE_TYPES[l]]);

#####################################################################

# ARRAYS FOR EACH TYPE OF DISTRIBUTIONS:
#
# 1D:
protpz          = [];
proten          = [];
protxi          = []
protpt          = []
proteta         = []
mpp             = []
mupz            = []
muen            = []
mupt            = []
ivm_mu          = []
mueta           = []
phopz           = []
phopt           = []
phoen           = []

'''# 3D:
DDDpt1pt2	= [];
DDDphi1phi2	= [];
DDDptsumphi	= [];
DDDpt1ptsum	= [];
DDDpt2ptsum	= [];
DDDmllptsum	= [];
DDDetaptsum	= [];
DDDetatheta	= [];
DDDetacost	= [];
DDDmllcost	= [];
DDDth1th2	= [];'''

# 2D:
DDmppmmumu      = []
DDxipximu     = []

# SETTING THE NUMBER OF DIGITS ON AXIS
TGaxis.SetMaxDigits(2)

# SORTING THE DISTRIBUTIONS WITHIN THE SETS:
# THE ARRAYS STORE THE LABELS FOR AXIS AND UNITS:

# 1D
histoslog        = [protpz,proten,protxi,protpt,proteta,mpp,mupz,muen,mupt,ivm_mu,mueta,phopz,phopt,phoen];
histoslog_label  = ["protpz","proten",'protxi','protpt','proteta','mpp',"mupz","muen","mupt",'ivm_mu','mueta','phopz','phopt','phoen'];
histoslog_axis   = ["p_{z}(p)","E(p)",'#chi(p)','p_{T}(p)','#eta(p^{+}p^{-})','M(p^{+}p^{-})',"p_{z}(#mu)","E(#mu)","p_{T}(#mu)",'M(#mu^{+}#mu^{-})','#eta(#mu^{+}#mu^{-})','p_{z}(#alpha)','p_{T}(#alpha)','E(#alpha)'];
histoslog_varx   = ["(GeV)","(GeV)",'','(GeV)','','(GeV)',"(GeV)","(GeV)","(GeV)",'(GeV)','','(GeV)','(GeV)','(GeV)'];

# 2D
DDlog         = [DDmppmmumu,DDxipximu] 
DDlog_label   = ["2Dmppmmumu",'2Dxipximu'];
DDlog_xaxis   = ["M(p^{+}p^{-})",'#xi(p^{+})']
DDlog_yaxis   = ["M(#mu^{+}#mu^{-})",'#xi(#mu^{+}#mu^{-})']
DDlog_varx    = ["(GeV)",'']
DDlog_vary    = ["(GeV)",'']

'''# 3D
legoslog         = [DDDpt1pt2,DDDphi1phi2,DDDptsumphi,DDDpt1ptsum,DDDpt2ptsum,DDDmllptsum,DDDetaptsum,DDDetatheta,DDDetacost,DDDmllcost,DDDth1th2];
legoslog_label   = ["3Dpt1pt2","3Dphi1phi2","3Dptsumphi","3Dpt1ptsum","3Dpt2ptsum","3Dmllptsum","3Detaptsum","3Detatheta","3Detacost","3Dmllcost","3Dth1th2"];
legoslog_xaxis   = ["p_{T}(x^{+})","#phi(x^{+})","p_{T}(x^{+}x^{-})","p_{T}(x^{#pm})","p_{T}(x^{#pm})","M(x^{+}x^{-})","#eta(x^{+}x^{-})","#eta(x^{+}x^{-})","#eta(x^{+}x^{-})","M(x^{+}x^{-})","#theta_{1}"];
legoslog_yaxis   = ["p_{T}(x^{-})","#phi(x^{-})","#phi(x^{#pm})","p_{T}(x^{+}x^{-})","p_{T}(x^{+}x^{-})","p_{T}(x^{+}x^{-})","p_{T}(x^{+}x^{-})","#theta(x^{+}x^{-})","cos#theta(x^{+}x^{-})","cos#theta(x^{+}x^{-})","#theta_{2}"];
legoslog_varz    = ["(nb/GeV^{2})","(nb)","(nb/GeV*deg)","(nb/GeV^{2})","(nb/GeV^{2})","(nb/GeV^{2})","(nb/GeV^{2})","(nb/deg)","(nb)","(nb/GeV)","(nb)"];
legoslog_varx    = ["(GeV)","(deg)","(GeV)","(GeV)","(GeV)","(GeV)","","","","(GeV)","(deg)"];
legoslog_vary    = ["(GeV)","(deg)","","(GeV)","(GeV)","(GeV)","(GeV)","(deg)","","","(deg)"];'''

# STARTING THE LOOP OVER FILES:
for i in range(len(FILES)):
    f = open(FILES[i],'r');
    print (f"Opening file {i}: {FILES[i]}");

    # SORTING THE DISTRIBUTIONS IN THE ARRAYS FOR EACH FILE:
    # EACH ARRAYS IS FORMATTED LIKE: array[] = [plots_file1, plots_file2, plots_file3, ...

    # 1D
    protpz.append(TH1D("1D_protpz"+"_"+PDF[i]       , "", 50,4300., 7200.))
    proten.append(TH1D("1D_proten"+"_"+PDF[i]       , "", 50,4300., 7200.))
    protxi.append(TH1D("1D_protxi"+"_"+PDF[i]       , "", 60,-0.003,0.03))
    protpt.append(TH1D("1D_protpt"+"_"+PDF[i]       , "", 50,-0.1, 1.))
    proteta.append(TH1D("1D_proteta"+"_"+PDF[i]       , "", 50,-20., 20.))
    mpp.append(TH1D("1D_mpp"+"_"+PDF[i]       , "", 60,0., 120.))
    mupz.append(TH1D("1D_mupz"+"_"+PDF[i]       , "", 50,-2500.,2500.))
    muen.append(TH1D("1D_muen"+"_"+PDF[i]       , "", 50,-100., 900.))
    mupt.append(TH1D("1D_mupt"+"_"+PDF[i]       , "", 50,-5., 40.0))
    ivm_mu.append(TH1D("1D_ivm_mu"+"_"+PDF[i]       , "", 60,0., 120.0))
    mueta.append(TH1D("1D_mueta"+"_"+PDF[i]       , "", 50,-15., 15.))
    phopz.append(TH1D("1D_phopz"+"_"+PDF[i]       , "", 50,4973., 4980.))
    phopt.append(TH1D("1D_phopt"+"_"+PDF[i]       , "", 50,0., 500.))
    phoen.append(TH1D("1D_phoen"+"_"+PDF[i]       , "", 50,-100., 2000.))
    #mopz.append(TH1D("1D_mupz"+"_"+PDF[i]       , "", 50,-2500.,2500.))
    #moen.append(TH1D("1D_muen"+"_"+PDF[i]       , "", 50,-100., 900.))
    #mopt.append(TH1D("1D_mupt"+"_"+PDF[i]       , "", 50,-5., 40.0))

    # 2D
    DDmppmmumu.append(TH2D('2D_mpp_mmumu_'+PDF[i]       , '', 50, 0., 1400., 50, 0., 1400.))
    DDxipximu.append(TH2D('2D_xip_ximu_'+PDF[i]     , '', 50, 0., 1., 50, 0., 1.))

    # 3D
    #DDDpt1pt2.append(TH2D("3D_pt1_pt2_"+PDF[i]      , "", 50,  0.,  70., 50, 0.,  70.));
    #DDDphi1phi2.append(TH2D("3D_phi1_phi2_"+PDF[i]  , "", 45,  0., 180., 45, 0., 180.));
    #DDDptsumphi.append(TH2D("3D_ptsum_phi_"+PDF[i]	, "", 50,  0., 160., 45, 0., 180.));
    #DDDpt1ptsum.append(TH2D("3D_pt1_ptsum_"+PDF[i]	, "", 50,  0.,  80., 50, 0., 120.));
    #DDDpt2ptsum.append(TH2D("3D_pt2_ptsum_"+PDF[i]	, "", 50,  0.,  80., 50, 0., 120.));
    #DDDmllptsum.append(TH2D("3D_mll_ptsum_"+PDF[i]	, "", 50,  0., 140., 50, 0., 120.));
    #DDDetaptsum.append(TH2D("3D_eta_ptsum_"+PDF[i]	, "", 50,-15.,  15., 50, 0., 100.));
    #DDDetatheta.append(TH2D("3D_eta_theta_"+PDF[i]  , "", 50,-10.,  10., 45, 0., 180.));
    #DDDetacost.append(TH2D("3D_eta_cost_"+PDF[i]    , "", 50,-15.,  15., 50,-1.,   1.));
    #DDDmllcost.append(TH2D("3D_mll_cost_"+PDF[i]    , "", 50,  0., 300., 50,-1.,   1.));
    #DDDth1th2.append(TH2D("3D_th1_th2_"+PDF[i]      , "", 45,  0., 180., 45, 0., 180.));


    # LOOP OVER LINES IN LHE SAMPLE:

    # RESET EVENT COUNTING:
    event  = 0;
    evPASS = 0;
    # START LOOP: <<<<<<<<<<<<<<<<!!! REMMEMBER TO COUNT CORRECTLY HERE
    if (i == 0):
        for j in range(87): # skip first 434 lines to avoid MG5 comments
            f.readline()
    elif (i == 1):
        for j in range(384): # skip first 431 lines to avoid MG5 comments
            f.readline()
    elif (i == 2):
        for j in range(7): # skip first 430 lines to avoid MG5 comments
            f.readline()
    elif (i == 3):
        for j in range(8): # skip first 440 lines to avoid MG5 comments
            f.readline();
    for line in f:
        # SKIP BLANK LINES:
        line = line.strip();
        if not line: continue;
        # STORE LINES INTO ARRAY:
        coll = line.split();
        # READ EVENT CONTENT:
        if coll[0] == "<event>":
            event += 1;
            # SET A SCREEN OUTPUT FOR CONTROL:
            if Nevt < 10000: evtsplit = 1000;
            else: evtsplit = 10000;
            perct = event / Nevt * 100.;
            if event%evtsplit==0: print ("Event %i [%.2f%%]" % (event,perct));
            elif event>Nevt: break;
        # 4-VECTORS FOR DECAY PRODUCTS:
        elif coll[0] == '22' and eval(coll[8]) > 0:
            dp = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dp.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '22' and eval(coll[8]) < 0:
            dm = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dm.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '13' and coll[1] == '1':
            dmu = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dmu.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '-13' and coll[1] == '1':
            damu = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            damu.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '90':
            dmo = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dmo.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '2212' and coll[1] == '1' and eval(coll[8]) > 0 and abs(eval(coll[8])) < (SQRTS/2):
            dpp = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dpp.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '2212' and coll[1] == '1' and eval(coll[8]) < 0 and abs(eval(coll[8])) < (SQRTS/2):
            dpm = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dpm.SetPxPyPzE(px,py,pz,en);
        # CLOSE EVENT AND FILL HISTOGRAMS:
        elif coll[0] == "</event>":
            # KINEMATICS OF DECAY PRODUCTS:
            if ( cuts and INNER                 # TRY EACH ONE
                and (dmu+damu).M() >= INVMCUTLOWER
                and (dmu+damu).M() <= INVMCUTUPPER
                and (dmu+damu).Pt() >= PTPAIRCUTLOWER
        	and (dmu+damu).Pt() <= PTPAIRCUTUPPER
                #and abs((damu+dmu).Eta()) <= ETAPAIRCUT
                and dmu.Pt() >= PTCUTLOWER
                and damu.Pt() >= PTCUTLOWER
                and dp.Pt() <= PTCUTUPPER
                and dm.Pt() <= PTCUTUPPER
                and abs(damu.Eta()) <= ETACUT
                and abs(dmu.Eta()) <= ETACUT):
                # 1D:
                #-------------------------Medidas dos prótons
                protpz[i].Fill(dpp.Pz());
                protpz[i].Fill(dpm.Pz());
                proten[i].Fill(dpp.E())
                proten[i].Fill(dpm.E())
                protxi[i].Fill(1-(dpp.Pz()/(SQRTS/2)))
                protxi[i].Fill(1-(dpm.Pz()/(-(SQRTS/2))))
                mpp[i].Fill(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS/2)))))*SQRTS)
                protpt[i].Fill(dpp.Pt())
                protpt[i].Fill(dpm.Pt())
                proteta[i].Fill(dpp.Eta())
                proteta[i].Fill(dpm.Eta())
                #-------------------------Medidas dos Múons
                mupz[i].Fill(dmu.Pz())
                muen[i].Fill(dmu.E())
                muen[i].Fill(damu.E())
                mupt[i].Fill(dmu.Pt())
                mupt[i].Fill(damu.Pt())
                ivm_mu[i].Fill((dmu+damu).M())
                mueta[i].Fill(dmu.Eta())
                mueta[i].Fill(damu.Eta())
                #-------------------------Medidas dos fótons
                phopt[i].Fill(dp.Pt())
                phopt[i].Fill(dm.Pt())
                phopz[i].Fill(dp.Pz())
                phopz[i].Fill(dm.Pz())
                phoen[i].Fill(dm.E())
                phoen[i].Fill(dp.E())
                #-------------------------Medidas do monopolo
                #mopz[i].Fill(dmo.Pz())
                #moen[i].Fill(dmo.E())
                #mopt[i].Fill(dmo.Pt())

                # 2D:
                DDmppmmumu[i].Fill(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS//2)))))*SQRTS, (dmu+damu).M())
                DDxipximu[i].Fill(1-(dpp.Pz()/(SQRTS/2)), (1/SQRTS)*(dmu.Pt()*exp(dmu.Eta())+damu.Pt()*exp(damu.Eta())))

                evPASS += 1;
            elif ( cuts and not INNER
                and (dp+dm).M() >= INVMCUTLOWER
                and (dp+dm).M() <= INVMCUTUPPER
                and (dp+dm).Pt() >= PTPAIRCUTLOWER
                and (dp+dm).Pt() <= PTPAIRCUTUPPER
                and abs((dp+dm).Eta()) >= ETAPAIRCUT
                and dp.Pt() >= PTCUTLOWER
                and dm.Pt() >= PTCUTLOWER
                and dp.Pt() <= PTCUTUPPER
                and dm.Pt() <= PTCUTUPPER):
                # 1D:
                invm_decay[i].Fill((dp+dm).M());# !
                pt_decay[i].Fill(dp.Pt());
                pt_decay[i].Fill(dm.Pt());
                ptsum_decay[i].Fill((dp+dm).Pt());
                eta_decay[i].Fill((dp).Eta());
                eta_decay[i].Fill((dm).Eta());
                phi_decay[i].Fill(dp.Phi());
                phi_decay[i].Fill(dm.Phi());
                E_decay[i].Fill(dp.E());
                E_decay[i].Fill(dm.E());
                dpt_decay[i].Fill(abs(dp.Pt()-dm.Pt()));
                dphi[i].Fill(abs(dp.DeltaPhi(dm))*180./3.141592);
                dphi_zoom[i].Fill(abs(dp.DeltaPhi(dm))*180./3.141592);
                acop_zoom[i].Fill((1. - abs(dp.DeltaPhi(dm))/3.141592));
                acop[i].Fill((1. - abs(dp.DeltaPhi(dm))/3.141592));
                # 3D:
                DDDetaptsum[i].Fill((dp+dm).Eta(),(dp+dm).Pt());
                DDDmllcost[i].Fill((dp+dm).M(),(dp+dm).CosTheta());
                DDDpt1pt2[i].Fill(dp.Pt(),dm.Pt());
                DDDphi1phi2[i].Fill(dp.Phi()*180./3.141592,dm.Phi()*180./3.141592);
                DDDpt1ptsum[i].Fill(dp.Pt(),(dp+dm).Pt());
                DDDpt2ptsum[i].Fill(dm.Pt(),(dp+dm).Pt());
                DDDmllptsum[i].Fill((dp+dm).M(),(dp+dm).Pt());
                DDDptsumphi[i].Fill((dp+dm).Pt(),dp.Phi()*180./3.141592);
                DDDptsumphi[i].Fill((dp+dm).Pt(),dm.Phi()*180./3.141592);
                DDDetatheta[i].Fill((dp+dm).Eta(),(dp+dm).Theta()*180./3.141592);
                DDDetacost[i].Fill((dp+dm).Eta(),(dp+dm).CosTheta());
                DDDth1th2[i].Fill(dp.Theta()*180./3.141592,dm.Theta()*180./3.141592);
                # 2D:
                DDpt1pt2[i].Fill(dp.Pt(),dm.Pt());
                DDphi1phi2[i].Fill(dp.Phi()*180./3.141592,dm.Phi()*180./3.141592);
                DDpt1ptsum[i].Fill(dp.Pt(),(dp+dm).Pt());
                DDpt2ptsum[i].Fill(dm.Pt(),(dp+dm).Pt());
                DDmllptsum[i].Fill((dp+dm).M(),(dp+dm).Pt());
                DDptsumphi[i].Fill((dp+dm).Pt(),dp.Phi()*180./3.141592);
                DDptsumphi[i].Fill((dp+dm).Pt(),dm.Phi()*180./3.141592);
                DDetatheta[i].Fill((dp+dm).Eta(),(dp+dm).Theta()*180./3.141592);
                DDetacost[i].Fill((dp+dm).Eta(),(dp+dm).CosTheta());
                DDmllcost[i].Fill((dp+dm).M(),(dp+dm).CosTheta());
                DDth1th2[i].Fill(dp.Theta()*180./3.141592,dm.Theta()*180./3.141592);
                evPASS += 1;
            elif not cuts:
                # 1D:
                #-------------------------Medidas dos prótons
                protpz[i].Fill(dpp.Pz());
                protpz[i].Fill(dpm.Pz());
                proten[i].Fill(dpp.E())
                proten[i].Fill(dpm.E())
                protxi[i].Fill(1-(dpp.Pz()/(SQRTS/2)))
                protxi[i].Fill(1-(dpm.Pz()/(-(SQRTS/2))))
                mpp[i].Fill(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS/2)))))*SQRTS)
                protpt[i].Fill(dpp.Pt())
                protpt[i].Fill(dpm.Pt())
                proteta[i].Fill(dpp.Eta())
                proteta[i].Fill(dpm.Eta())
                #-------------------------Medidas dos Múons
                mupz[i].Fill(dmu.Pz())
                muen[i].Fill(dmu.E())
                muen[i].Fill(damu.E())
                mupt[i].Fill(dmu.Pt())
                mupt[i].Fill(damu.Pt())
                ivm_mu[i].Fill((dmu+damu).M())
                mueta[i].Fill(dmu.Eta())
                mueta[i].Fill(damu.Eta())
                #-------------------------Medidas dos fótons
                phopt[i].Fill(dp.Pt())
                phopt[i].Fill(dm.Pt())
                phopz[i].Fill(dp.Pz())
                phopz[i].Fill(dm.Pz())
                phoen[i].Fill(dm.E())
                phoen[i].Fill(dp.E())
                #-------------------------Medidas do monopolo
                #mopz[i].Fill(dmo.Pz())
                #moen[i].Fill(dmo.E())
                #mopt[i].Fill(dmo.Pt())

                # 2D:
                DDmppmmumu[i].Fill(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS//2)))))*SQRTS, (dmu+damu).M())
                DDxipximu[i].Fill(1-(dpp.Pz()/(SQRTS/2)), (1/SQRTS)*(dmu.Pt()*exp(dmu.Eta())+damu.Pt()*exp(damu.Eta())))

                '''# 3D:
                DDDetaptsum[i].Fill((dp+dm).Eta(),(dp+dm).Pt());
                DDDmllcost[i].Fill((dp+dm).M(),(dp+dm).CosTheta());
                DDDpt1pt2[i].Fill(dp.Pt(),dm.Pt());
                DDDphi1phi2[i].Fill(dp.Phi()*180./3.141592,dm.Phi()*180./3.141592);
                DDDpt1ptsum[i].Fill(dp.Pt(),(dp+dm).Pt());
                DDDpt2ptsum[i].Fill(dm.Pt(),(dp+dm).Pt());
                DDDmllptsum[i].Fill((dp+dm).M(),(dp+dm).Pt());
                DDDptsumphi[i].Fill((dp+dm).Pt(),dp.Phi()*180./3.141592);
                DDDptsumphi[i].Fill((dp+dm).Pt(),dm.Phi()*180./3.141592);
                DDDetatheta[i].Fill((dp+dm).Eta(),(dp+dm).Theta()*180./3.141592);
                DDDetacost[i].Fill((dp+dm).Eta(),(dp+dm).CosTheta());
                DDDth1th2[i].Fill(dp.Theta()*180./3.141592,dm.Theta()*180./3.141592);'''

        # End of loop over lines
        if evPASS >= Nmax: break
    if cuts: print ("Events passing acceptance: %i/%i" % (evPASS,event));
        #print ("Integral of %s: %.6f nb" % (PDF[i],evPASS*xsec[i]/event));
# End of loop over files

#############################################################
#
#-----------------------KS TEST------------------------------
#
#############################################################

with open(f'ks-test.txt', 'w') as f:
    f.write('>>>>>>>KOLMOGOROV-SMIRNOV TEST<<<<<<<\n\n')
    f.write(f'Invariant-Mass protons\n')
    for i in range(4):
        for j in range(i+1, 4):
            ks = mpp[i].KolmogorovTest(mpp[j])
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')
    f.write(f'\nInvariant-Mass leptons\n')
    for i in range(4):
        for j in range(i+1, 4):
            ks = ivm_mu[i].KolmogorovTest(ivm_mu[j])
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')
    f.write(f'\n\u03A7 of protons\n')
    for i in range(4):
        for j in range(i+1, 4):
            ks = protxi[i].KolmogorovTest(protxi[j])
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')

#############################################################
#
#-----------------------END OF KS TEST-----------------------
#
#############################################################


# Starting Drawing step:

# Defining the top label in the plots:
plotlabel = TPaveText(0.50,0.91,0.84,0.95,"NDC");
plotlabel.SetTextAlign(33);
plotlabel.SetTextColor(1);
plotlabel.SetFillColor(0);
plotlabel.SetBorderSize(0);
plotlabel.SetTextSize(0.035);
plotlabel.SetTextFont(42);
plotlabel.AddText("MadGraphv5 #bullet #sqrt{s}=14 TeV #bullet "+EVTINPUT+" evt");

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
    for m in range(len(FILES)):
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
    for m in range(len(FILES)):
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
FILEROOT.Write();

#####################################################################
#
# C'ESTI FINI
#
#####################################################################
