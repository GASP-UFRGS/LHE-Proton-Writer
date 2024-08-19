import sys
import numpy as np
import random
from ROOT import TLorentzVector 
'''
Writes the protons into an LHE file from  MadGraph or Superchic 
'''

def Draw():
    # introduces pileup protons
    # the fractional momentum loss (xi) of the protons is assumed to follow a f(xi) = 1/xi distribution.
    xi_min = 0.015      # minimum fractional momentum loss of protons detected by the CT-PPS
    xi_max = 0.2        # maximum fractional momentum loss of protons detected by the CT-PPS
    mu = 1              # mean for the gaussian blur
    sigma = 0.02        # standard deviation for the gaussian blur, taken as 2%.
    PUexpected = 50     # expected value of pileup events. 
    
    p_s = 0.0115695     # calculated probability of finding a pileup proton in one of the PPS arms.
    
    xi1_list = list()
    xi2_list = list()
    puEvents = np.random.poisson(PUexpected)    # the number of pileup events is assumed to follow a poisson distribution.
    for k in range(puEvents):
        tag = np.random.uniform(0,1,size=None)
    
        if tag <= p_s:
            xi1_list.append( xi_min*pow( xi_max/xi_min, np.random.uniform(0,1,size=None) )*random.gauss(mu, sigma) )

        elif tag > p_s and tag <= 2*p_s:
            xi2_list.append( xi_min*pow( xi_max/xi_min, np.random.uniform(0,1,size=None) )*random.gauss(mu, sigma) )

    return [xi1_list, xi2_list]
    
# Receives file path, generator of origin and particle IDs as arguments
if len(sys.argv) < 4:
    print('Missing arguments')
    print('Syntax: python3 LHE-Proton-Writer.py <path of .lhe file> <generator of origin> <IDs>')
    sys.exit()

path = sys.argv[1]
generator = sys.argv[2].lower()
ID = sys.argv[3:]
new = 'new_'+path

if not generator == 'madgraph' or generator == 'superchic':
    print('<<'+generator+'>> generator unsupported. Only "madgraph" and "superchic" are supported. Exiting.')
    sys.exit()

# Flag for printing the Invariant Mass of protons and leptons
printivm = False

# Option to include pileup protons:
pileup = True

# Setting flags according to chosen generator
if(generator == 'madgraph'):
    flag0 = '<event>\n'
    flag1 = '</event>\n'
    header = '<init>\n'
    end = '</LesHouchesEvents>\n'
if(generator == 'superchic'):
    flag0 = ' <event>\n'
    flag1 = ' </event>\n'
    header = ' <init>\n'
    end = ' </LesHouchesEvents>\n'

# Create list of lines from the LHE file and if needed finding the beginning of events 
with open(path, 'r+') as f:
    l = f.readline()
    while l != header:
        l = f.readline()
    beamenergy = int(float(f.readline().split()[2]))

pzini = beamenergy		# protons inicial pz 
eini = beamenergy		# protons inicial energy

# Proton rest mass: 0.9382720882E+00    # from: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-phys-constants.pdf
m0 = 0.9382720882


event = []
lvector = TLorentzVector()
with open(path, 'r+') as f, open(new, 'w') as new:
    l = f.readline()
    while l != end:
        l = f.readline()
        if l == flag0:
            event.append(l)
            while l != flag1:
                l = f.readline()
                event.append(l)
            for i in range(len(event)):
                if event[i].split()[0] in ID and event[i].split()[1] == '-1':
                    # Updating number of particles on header
                    neweventheader = event[1].split()
                    if generator == 'superchic':
                        neweventheader[0] = ' '+str(int(neweventheader[0])+1)
                        event.pop(1)
                        event.insert(1, '    '.join(neweventheader)+'\n')
                    if generator == 'madgraph':
                        neweventheader[0] = ' '+str(int(neweventheader[0])+1)+'     '
                        event.pop(1)
                        event.insert(1, ' '.join(neweventheader)+'\n')
                    # Creating the 4-momentum vector
                    line = event[i].split()
                    px = f'{-eval(line[6]):.9e}' if -eval(line[6]) < 0 else f'+{-eval(line[6]):.9e}'
                    py = f'{-eval(line[7]):.9e}' if -eval(line[7]) < 0 else f'+{-eval(line[7]):.9e}'
                    pzf = eval(line[8])
                    sign = (pzf/abs(pzf))
                    pzp = pzini*sign - pzf
                    ef = eval(line[9])
                    ep = eini - ef
                    if generator == 'superchic':
                        if sign > 0: event.insert(i+1, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} +{pzp:.9e}  {ep:.9e}  {m0:.9e} 0. 1.\n')
                        else: event.insert(i+1, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} {pzp:.9e}  {ep:.9e}  {m0:.9e} 0. -1.\n')
                    if generator == 'madgraph':
                        if sign > 0: event.insert(i+1, ' '*5+f'2212{" "*2}1    0    0    0    0 +{0:.10e} +{0:.10e} +{pzp:.10e} {ep:.10e} {m0:.10e} {0:.4e} {1:.4e}\n')
                        else: event.insert(i+1, ' '*5+f'2212{" "*2}1    0    0    0    0 -{0:.10e} -{0:.10e} {pzp:.10e} {ep:.10e} {m0:.10e} {0:.4e} {-1:.4e}\n')
                    PU_protons = Draw() # Draws pileup protons
                    # the 6th number is the xi, while px = py = m0 = 0, and spin/helicity = 9 for identification as a pileup proton. 
                    if pileup:
                        for xi in PU_protons[0]: # loop over 'left PPS arm xis'
                            if generator == 'superchic':
                                event.insert(i+1, ' '*13+f'2212{" "*8}1    0    0    0    {xi:.9e} {0:.9e} {0:.9e} +{(1-xi)*pzini:.9e}  {(1-xi)*pzini:.9e}  {0:.9e} 0. 9.\n')
                            if generator == 'madgraph':
                                event.insert(i+1,' '*5+f'2212{" "*2}1    0    0    0    {xi:.10e} +{0:.5e} +{0:.5e} +{(1-xi)*pzini:.10e} {(1-xi)*pzini:.10e} {0:.5e} {0:.4e} {9:.4e}\n')
                        for xi in PU_protons[1]: # loop over 'right PPS arm xis'
                            if generator == 'superchic':
                                event.insert(i+1,' '*13+f'2212{" "*8}1    0    0    0    {xi:.9e} {0:.9e} {0:.9e} -{(1-xi)*pzini:.9e}  {(1-xi)*pzini:.9e}  {0:.9e} 0. 9.\n')
                            if generator == 'madgraph':
                                event.insert(i+1,' '*5+f'2212{" "*2}1    0    0    0    {xi:.10e} +{0:.5e} +{0:.5e} -{(1-xi)*pzini:.10e} {(1-xi)*pzini:.10e} {0:.5e} {0:.4e} {9:.4e}\n')

            for nl in event:
                new.write(nl)
                if(printivm):
                   line = nl.split()
                   if(len(line) > 1 and line[0] == '2212'):
                       lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                       print(f'M(p) = {lvector.M()}')
                   if(len(line) > 1 and line[0] == '13'):
                       lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                       print(f'M(\u03BC+) = {lvector.M()}')
                   if(len(line) > 1 and line[0] == '-13'):
                       lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                       print(f'M(\u03BC) = {lvector.M()}\n*-*-*-*-*-*-*-*-*-*-*')
            event.clear()
        else:
            new.write(l)
