"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0
"""

import cantera as ct

# Parâmetros de Simulação
# REAÇÃO BALANCEADA
# Reagentes: CH4 + CO + H2 + 3.1 (O2 + 3.76 N2) + H2O + CO2 -> 
# Produtos: 3CO2 + 5.826 N2 + 0.1 O2 + 4 H2O

p = ct.one_atm  # Pressão [Pa]
Tin = (25 + 273.15)  # temperatura inicial do gás [K]
#reactants = 'CO:.2215, H2: 0.0546, H2O:0.0988, CH4:0.0102, N2:0.4673, CO2:0.1476'  # Composição da pré mistura
reactants = "CO:.3713, H2:.0915, H2O:.1656, CH4:.0171, CO2:.2474, O2:.0967, N2:0.4161"
width = 0.03  # tamanho (m)
loglevel = 1  # amount of diagnostic output (0 to 8) ??

# Solução gasosa para análise
gas = ct.Solution('gri30.yaml') # verificar mecanismo
gas.TPY = Tin, p, reactants # condições iniciais da mistura 

# Objeto chama do Cantera
f = ct.FreeFlame(gas, width=width) # passa o gás e determina o tamanho da chama 
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12) # verificar parâmetros 
f.show_solution() # exibe a solução

# Solve with mixture-averaged transport model (verificar na api )
f.transport_model = 'Mix' 
f.solve(loglevel=loglevel, auto=True) # chama o solver 

# Salva o cálculo em formato xml
f.save('adiabatic_flame.xml', 
    'mix',
    'solution with mixture-averaged transport'
)

f.show_solution()
print('mixture-averaged flame temperature = {0:7f} K'.format(f.T[0]))

# Resolve a chama utilizando algoritmo de trasnporte multicomponete.
f.transport_model = 'Multi'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show_solution()
print('multicomponent flametemperature = {0:7f} K'.format(f.T[0]))
f.save(
    'adiabatic_flame.xml', 
    'multi',
    'solution with multicomponent transport'
)

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('adiabatic_flame.csv', quiet=False)