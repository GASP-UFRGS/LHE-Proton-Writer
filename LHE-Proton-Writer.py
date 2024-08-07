import sys
from ROOT import TLorentzVector 
'''
Writes the protons into an LHE file from  MadGraph or Superchic 
'''
# Receives particle ID, file name and generator of origin as arguments
if len(sys.argv) < 4:
    print('Missing argument, correct syntax as follows:')
    print('python3 LHE-Proton-Writer.py <path of .lhe file> <generator of origin> <IDs>')
    sys.exit()

path = sys.argv[1]
generator = sys.argv[2].lower()
ID = sys.argv[3:]
new = 'new'+path

# Flag for printing the Invariant Mass of protons and leptons
printivm = False

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
