import sys
from ROOT import *
'''
Writes the protons into an LHE file from Superchic 
'''
# Ele recebe o arquivo como argumento depois do python3
if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    # Caso não receba o argumento usa esse como padrão
    path = 'evrectest.dat'
new = 'amostras/new'+path

# Aqui serve pra encontrar o inicio dos eventos no .lhe, caso seja relevante, e inicializar a lista com as linhas 
with open('amostras/'+path, 'r+') as f:
    linhas = f.readlines()
    #print(linhas)
    inicio = linhas.index(' <event>\n')	

# Aqui serve pra inserir os dados
pzini = 7000		# pz inicial dos protons
eini = 7000		# energia inicial dos protons

# Massa de repouso do proton: 0.9382720882E+00

lvector = TLorentzVector()
with open(new, 'w') as new:
    i = 0
    while True:
        if linhas[i] == ' <event>\n':
            j = 0
            while linhas[i+j] != ' </event>\n':
                if linhas[i+j].split()[0] == '22':
                    linha = linhas[i+j].split()
                    px = f'{-eval(linhas[i+j].split()[6]):.9e}' if -eval(linha[6]) < 0 else f'+{-eval(linhas[i+j].split()[6]):.9e}'
                    py = f'{-eval(linhas[i+j].split()[7]):.9e}' if -eval(linha[7]) < 0 else f'+{-eval(linhas[i+j].split()[7]):.9e}'
                    pzf = eval(linha[8])
                    sign = (pzf/abs(pzf))
                    pzp = pzini*sign - pzf
                    ef = eval(linha[9])
                    ep = eini - ef
                    if sign > 0: linhas.insert(i+j+2, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} +{pzp:.9e}  {ep:.9e}  0.938272088E+00 0. 1.\n')
                    else: linhas.insert(i+j+2, ' '*13+f'2212{" "*8}1    0    0    0    0 {px} {py} {pzp:.9e}  {ep:.9e}  0.938272088E+00 0. -1.\n')
                j += 1
        i += 1
        if i == len(linhas): break
    for line in linhas:
        linha = line.strip().split()
        if(len(linha) > 1 and linha[0] == '2212'):
            lvector.SetPxPyPzE(float(linha[6]), float(linha[7]), float(linha[8]), float(linha[9]))
            print(f'M(p) = {lvector.M()}')
        if(len(linha) > 1 and linha[0] == '13'):
            lvector.SetPxPyPzE(float(linha[6]), float(linha[7]), float(linha[8]), float(linha[9]))
            print(f'M(\u03BC+) = {lvector.M()}')
        if(len(linha) > 1 and linha[0] == '-13'):
            lvector.SetPxPyPzE(float(linha[6]), float(linha[7]), float(linha[8]), float(linha[9]))
            print(f'M(\u03BC) = {lvector.M()}\n*-*-*-*-*-*-*-*-*-*-*')
        new.write(line[1:])

