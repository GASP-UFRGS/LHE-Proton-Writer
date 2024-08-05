import sys
from ROOT import TLorentzVector 
'''
Writes the protons into an LHE file from  MadGraph or Superchic 
'''
# Receives particle ID, file name and generator of origin as arguments
if len(sys.argv) != 4:
    print('Missing arguments')
    print('Syntax: python3 LHE-Proton-Writer.py <id> <input> <generator>')
    sys.exit()

ID = str(sys.argv[1])
path = sys.argv[2]
generator = sys.argv[3].lower()
if not generator == 'madgraph' or generator == 'superchic':
    print('<<'+generator+'>> generator unsupported. Exiting.')
    sys.exit()
new = 'new_'+path

# Flag for printing the Invariant Mass of protons and leptons
printivm = False

# Fill lines:
def func_fill(_gen):
    if _gen == 'superchic':
        if sign > 0: lines.insert(i+j+2, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} +{pzp:.9e}  {ep:.9e}  {m0:.9e} 0. 1.\n')
        else: lines.insert(i+j+2, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} {pzp:.9e}  {ep:.9e}  {m0:.9e} 0. -1.\n')
    elif _gen == 'madgraph':
        if sign > 0: lines.insert(i+j+2, ' '*5+f'2212{" "*2}1    0    0    0    0 +{0:.10e} +{0:.10e} +{pzp:.10e} {ep:.10e} {m0:.10e} {0:.4e} {1:.4e}\n')
        else: lines.insert(i+j+2, ' '*5+f'2212{" "*2}1    0    0    0    0 -{0:.10e} -{0:.10e} {pzp:.10e} {ep:.10e} {m0:.10e} {0:.4e} {-1:.4e}\n')

# Setting flags according to chosen generator
if(generator == 'madgraph'):
    flag0 = '<event>\n'
    flag1 = '</event>\n'
    header = '</header>\n'
if(generator == 'superchic'):
    flag0 = ' <event>\n'
    flag1 = ' </event>\n'
    header = ' </header>\n'

# Create list of lines from the LHE file and if needed finding the beginning of events 
with open(path, 'r+') as f:
    lines = f.readlines()
    if(generator == 'madgraph'):
        lines.index(flag0)
        beamenergy = int(float(lines[lines.index(header)+2].split()[2]))
    if(generator == 'superchic'):
        lines.index(flag0)
        beamenergy = int(float(lines[lines.index(header)+2].split()[2]))

pzini = beamenergy		# protons inicial pz 
eini = beamenergy		# protons inicial energy

# Proton rest mass: 0.9382720882E+00    # from: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-phys-constants.pdf
m0 = 0.9382720882

lvector = TLorentzVector()
with open(new, 'w') as new:
    i = 0
    while True:
        if lines[i] == flag0:
            j = 0
            while lines[i+j] != flag1:
                if lines[i+j].split()[0] == ID and lines[i+j].split()[1] == '-1':
                    line = lines[i+j].split()
                    px = f'{-eval(lines[i+j].split()[6]):.9e}' if -eval(line[6]) < 0 else f'+{-eval(lines[i+j].split()[6]):.9e}'
                    py = f'{-eval(lines[i+j].split()[7]):.9e}' if -eval(line[7]) < 0 else f'+{-eval(lines[i+j].split()[7]):.9e}'
                    print('oi')
                    pzf = eval(line[8])
                    sign = (pzf/abs(pzf))
                    pzp = pzini*sign - pzf
                    ef = eval(line[9])
                    ep = eini - ef
                    func_fill(generator)
                j += 1
            neweventheader = lines[i+1].split()
            if generator == 'superchic':
                neweventheader[0] = ' '+str(int(neweventheader[0])+2)
                lines.pop(i+1)
                lines.insert(i+1, '    '.join(neweventheader)+'\n')
            if generator == 'madgraph':
                neweventheader[0] = ' '+str(int(neweventheader[0])+2)+'     '
                lines.pop(i+1)
                lines.insert(i+1, ' '.join(neweventheader)+'\n')
        i += 1
        if i == len(lines): break
    for line in lines:
        if(printivm):
            line = line.strip().split()
            if(len(line) > 1 and line[0] == '2212'):
                lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                print(f'M(p) = {lvector.M()}')
            if(len(line) > 1 and line[0] == '13'):
                lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                print(f'M(\u03BC+) = {lvector.M()}')
            if(len(line) > 1 and line[0] == '-13'):
                lvector.SetPxPyPzE(float(line[6]), float(line[7]), float(line[8]), float(line[9]))
                print(f'M(\u03BC) = {lvector.M()}\n*-*-*-*-*-*-*-*-*-*-*')
        new.write(line)


