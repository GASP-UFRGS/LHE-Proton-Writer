from __future__ import division
from subprocess import call
from math import *
from ROOT import *
from array import array
import numpy as np

#####################################################################
# GGS (CERN-CMS/UFRGS) ---
# the muons are collected considering the ID codes in the event
# sample produced with SuperCHICv2 in LHE format.
#####################################################################

#####################################################################
# USER INPUT:

rootinput = np.loadtxt('root_input.txt', unpack=True, dtype=str, delimiter='=')

# CROSS SECTION(S) (pb):
xsec    = eval(rootinput[1][0]) #FIXME

# PDF "_"+LABEL FOR OUTPUT FILES:
JOB     = eval(rootinput[1][1])
PDF     = eval(rootinput[1][2]) #FIXME
scale   = eval(rootinput[1][3]) #bug, use False, carefull with Nevts
cuts    = eval(rootinput[1][4])
setLog  = eval(rootinput[1][5])
filled  = eval(rootinput[1][6])
stacked = eval(rootinput[1][7])
data    = eval(rootinput[1][8])

# PARTICLES OF INTEREST
IDS     = eval(rootinput[1][15])


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
FILES   = eval(rootinput[1][13])#FIXME


# EVENT SAMPLE INPUT:
Nevt    = eval(rootinput[1][10]) #FIXME
Nmax    = eval(rootinput[1][11]) # number of max events to obtain from samples
EVTINPUT= str(int(Nevt/1000))+"k";
SQRTS   = eval(rootinput[1][12]) # in GeV

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
protpz       = []
proten       = []
protxi       = []
protpt       = []
proteta      = []
mpp          = []
mupz         = []
muen         = []
mupt         = []
ivmmu       = []
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


# SORTING THE DISTRIBUTIONS WITHIN THE SETS:

# 1D
histoslog  = [protpz,proten,protxi,protpt,proteta,mpp,mupz,muen,mupt,ivmmu,mueta,phopz,phopt,phoen,phoivm,phopsrap2,phopsrap1,phoY]

# 2D
DDlog      = [DDmppmmumu,DDxipximu] 

#------------ Lists for KS test ----------------------

KS_ivm_pp = list([] for i in range(len(FILES)))
KS_ivm_mu = list([] for i in range(len(FILES)))
KS_protxi = list([] for i in range(len(FILES)))

#-----------------------------------------------------

# STARTING THE LOOP OVER FILES:
for i in range(len(FILES)):
    f = open(FILES[i],'r');
    print (f"Opening file {i}: {FILES[i]}");

    # SORTING THE DISTRIBUTIONS IN THE ARRAYS FOR EACH FILE:
    # EACH ARRAYS IS FORMATTED LIKE: array[] = [plots_file1, plots_file2, plots_file3, ...

    # 1D
    protpz.append(TH1D("1D_protpz"+"_"+PDF[i]       , "", 50,4300., 7200.))
    proten.append(TH1D("1D_proten"+"_"+PDF[i]       , "", 50,4300., 7200.))
    protxi.append(TH1D("1D_protxi"+"_"+PDF[i]       , "", 50,-0.003,0.03))
    protpt.append(TH1D("1D_protpt"+"_"+PDF[i]       , "", 50,-0.1, 1.))
    proteta.append(TH1D("1D_proteta"+"_"+PDF[i]       , "", 50,-20., 20.))
    mpp.append(TH1D("1D_mpp"+"_"+PDF[i]       , "", 50,0., 120.))
    mupz.append(TH1D("1D_mupz"+"_"+PDF[i]       , "", 50,-2500.,2500.))
    muen.append(TH1D("1D_muen"+"_"+PDF[i]       , "", 50,-100., 900.))
    mupt.append(TH1D("1D_mupt"+"_"+PDF[i]       , "", 50,-5., 40.0))
    ivmmu.append(TH1D("1D_ivmmu"+"_"+PDF[i]       , "", 50,0., 120.0))
    mueta.append(TH1D("1D_mueta"+"_"+PDF[i]       , "", 50,-15., 15.))
    phopz.append(TH1D("1D_phopz"+"_"+PDF[i]       , "", 50,-10000., 10000.))
    phopt.append(TH1D("1D_phopt"+"_"+PDF[i]       , "", 50,0., 8000.))
    phoen.append(TH1D("1D_phoen"+"_"+PDF[i]       , "", 50,-100., 5000.))
    phoivm.append(TH1D("1D_phoivm"+"_"+PDF[i]     , "", 50, -1000, 10000))
    phopsrap2.append(TH1D('1D_phopsrap2'+'_'+PDF[i] , '', 50, -1, 1))
    phopsrap1.append(TH1D('1D_phopsrap1'+'_'+PDF[i] , '', 50, -10, 10))
    phoY.append(TH1D('1D_phoY'+'_'+PDF[i]         , '', 50, -10000, 10000))
    #mopz.append(TH1D("1D_mupz"+"_"+PDF[i]       , "", 50,-2500.,2500.))
    #moen.append(TH1D("1D_muen"+"_"+PDF[i]       , "", 50,-100., 900.))
    #mopt.append(TH1D("1D_mupt"+"_"+PDF[i]       , "", 50,-5., 40.0))

    # 2D
    DDmppmmumu.append(TH2D('2D_DDmppmmumu_'+PDF[i]       , '', 50, 0., 1400., 50, 0., 1400.))
    DDxipximu.append(TH2D('2D_DDxipximu_'+PDF[i]     , '', 50, 0., 1., 50, 0., 1.))

    # LOOP OVER LINES IN LHE SAMPLE:

    # Flags for differentiating particles
    first_photon = False
    first_muon    = False
    first_proton = False

    # RESET EVENT COUNTING:
    event  = 0;
    evPASS = 0;
    # START LOOP: <<<<<<<<<<<<<<<<!!! REMMEMBER TO COUNT CORRECTLY HERE
    if (i == 0):
        for j in range(87): # skip first 434 lines to avoid comments
            f.readline()
    elif (i == 1):
        for j in range(401): # skip first 431 lines to avoid comments
            f.readline()
    elif (i == 2):
        for j in range(7): # skip first 430 lines to avoid comments
            f.readline()
    elif (i == 3):
        for j in range(8): # skip first 440 lines to avoid comments
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
        elif coll[0] == '22' and coll[1] == '1' and first_photon == False:
            dp = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dp.SetPxPyPzE(px,py,pz,en);
            first_photon = True
        elif coll[0] == '22' and coll[1] == '1' and first_photon:
            dm = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dm.SetPxPyPzE(px,py,pz,en);
            first_photon = False 
        elif coll[0] == '13' and coll[1] == '1' and first_muon == False:
            dmu = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dmu.SetPxPyPzE(px,py,pz,en);
            first_muon = True
        elif coll[0] == '-13' and coll[1] == '1' and first_muon:
            damu = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            damu.SetPxPyPzE(px,py,pz,en);
            first_muon = False
        elif coll[0] == '90':
            dmo = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dmo.SetPxPyPzE(px,py,pz,en);
        elif coll[0] == '2212' and coll[1] == '1' and first_proton == False and abs(eval(coll[8])) < (SQRTS/2):
            dpp = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dpp.SetPxPyPzE(px,py,pz,en)
            first_proton = True
        elif coll[0] == '2212' and coll[1] == '1' and first_proton and abs(eval(coll[8])) < (SQRTS/2):
            dpm = TLorentzVector();
            px = float(coll[6]);
            py = float(coll[7]);
            pz = float(coll[8]);
            en = float(coll[9]);
            dpm.SetPxPyPzE(px,py,pz,en)
            first_proton = False
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
                #and dp.Pt() <= PTCUTUPPER
                #and dm.Pt() <= PTCUTUPPER
                and abs(damu.Eta()) <= ETACUT
                and abs(dmu.Eta()) <= ETACUT):
                # 1D:
                #-------------------------Proton mesurements
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
                #-------------------------Múon mesurements
                if 13 in IDS or -13 in IDS:
                    mupz[i].Fill(dmu.Pz())
                    muen[i].Fill(dmu.E())
                    muen[i].Fill(damu.E())
                    mupt[i].Fill(dmu.Pt())
                    mupt[i].Fill(damu.Pt())
                    ivmmu[i].Fill((dmu+damu).M())
                    mueta[i].Fill(dmu.Eta())
                    mueta[i].Fill(damu.Eta())
                #-------------------------Photon mesurements
                if 22 in IDS:
                    phopt[i].Fill(dp.Pt())
                    phopt[i].Fill(dm.Pt())
                    phopz[i].Fill(dp.Pz())
                    phopz[i].Fill(dm.Pz())
                    phoen[i].Fill(dm.E())
                    phoen[i].Fill(dp.E())
                #-------------------------Monopole mesurements
                if 90 in IDS:
                    mopz[i].Fill(dmo.Pz())
                    moen[i].Fill(dmo.E())
                    mopt[i].Fill(dmo.Pt())
  
                #-------------------------Mesurements for the KS test
                if 13 in IDS and -13 in IDS:
                    KS_ivm_pp[i].append(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS/2)))))*SQRTS)
                    KS_protxi[i].append(1-(dpp.Pz()/(SQRTS/2)))
                    KS_protxi[i].append(1-(dpm.Pz()/(-(SQRTS/2))))
                    KS_ivm_mu[i].append((dmu+damu).M())
  
                # 2D:
                if 13 in IDS and -13 in IDS:
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
                if 13 in IDS and -13 in IDS:
                    mupz[i].Fill(dmu.Pz())
                    muen[i].Fill(dmu.E())
                    muen[i].Fill(damu.E())
                    mupt[i].Fill(dmu.Pt())
                    mupt[i].Fill(damu.Pt())
                    ivmmu[i].Fill((dmu+damu).M())
                    mueta[i].Fill(dmu.Eta())
                    mueta[i].Fill(damu.Eta())
                #-------------------------Medidas dos fótons
                if 22 in IDS:
                    phopt[i].Fill(dp.Pt())
                    phopt[i].Fill(dm.Pt())
                    phopz[i].Fill(dp.Pz())
                    phopz[i].Fill(dm.Pz())
                    phoen[i].Fill(dp.E())
                    phoen[i].Fill(dm.E())
                    phoivm[i].Fill((dp+dm).M())
                    phopsrap2[i].Fill((dp+dm).Eta())
                    phopsrap1[i].Fill(dp.Eta())
                    phopsrap1[i].Fill(dm.Eta())               
                    phoY[i].Fill(dp.Y())
                    phoY[i].Fill(dm.Y())
                #-------------------------Medidas do monopolo
                if 90 in IDS:
                    mopz[i].Fill(dmo.Pz())
                    moen[i].Fill(dmo.E())
                    mopt[i].Fill(dmo.Pt())
  
                # 2D:
                if 13 in IDS and -13 in IDS:
                    DDmppmmumu[i].Fill(sqrt((1-(dpp.Pz()/(SQRTS/2)))*(1-(dpm.Pz()/(-(SQRTS//2)))))*SQRTS, (dmu+damu).M())
                    DDxipximu[i].Fill(1-(dpp.Pz()/(SQRTS/2)), (1/SQRTS)*(dmu.Pt()*exp(dmu.Eta())+damu.Pt()*exp(damu.Eta())))

                evPASS += 1
        # End of loop over lines
        if evPASS >= Nmax: break   
    #print(phoivm[i].Integral())
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
            ks = TMath.KolmogorovTest(len(KS_ivm_pp[i]), array('d', sorted(KS_ivm_pp[i])), len(KS_ivm_pp[j]), array('d', sorted(KS_ivm_pp[j])), 'D') 
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')
    f.write(f'\nInvariant-Mass leptons\n')
    for i in range(4):
        for j in range(i+1, 4):
            ks = TMath.KolmogorovTest(len(KS_ivm_mu[i]), array('d', sorted(KS_ivm_mu[i])), len(KS_ivm_mu[j]), array('d', sorted(KS_ivm_mu[j])), 'D')
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')
    f.write(f'\n\u03A7 of protons\n')
    for i in range(4):
        for j in range(i+1, 4):
            ks = TMath.KolmogorovTest(len(KS_protxi[i]), array('d', sorted(KS_protxi[i])), len(KS_protxi[j]), array('d', sorted(KS_protxi[j])), 'D')
            f.write(f'{PDF[i]:>10} X {PDF[j]:<10}: {ks}\n')

#############################################################
#
#-----------------------END OF KS TEST-----------------------
#
#############################################################
FILEROOT.Write()

#####################################################################
#
# C'ESTI FINI
#
#####################################################################
