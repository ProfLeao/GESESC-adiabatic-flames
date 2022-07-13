import pandas as pd
from CoolProp.CoolProp import PropsSI 
import matplotlib.pyplot as plt
import numpy as np

from gasmistura_pkg import mm_aparente_mistura
from gasmistura_pkg import corr_vazao_normal
from gasmistura_pkg import vaz_combustao
from gasmistura_pkg import temp_adiabatica




df_composicao = pd.read_json("composicao.json", orient='split')

# Intervalo de vazões a serem simuladas em Nm³/h
# Considere o limite superior aberto. 
vaz_nm3 = [14e3, 17e3]
# Discretização do intervalo
del_vaz_nm3 = 500 # Nm³/h
# Inicialização do vetor de vazões para simulações. 
vaz_nm3 = np.arange(vaz_nm3[0], vaz_nm3[1], del_vaz_nm3)

# Temperatura de entrada dos gases de auto-forno
temp_gas_af = 95 # °C
temp_gas_afSI = temp_gas_af + 273.15 # K



mma, r_mistura, z_mistura = mm_aparente_mistura(
    compostos = list(df_composicao.columns),
    mmolar = df_composicao.loc["Massa Molar - (g/mol)"].to_numpy(),
    fracao = (df_composicao.loc["% em Massa"]).to_numpy()/100,
    t_fracao = "fmol"
)

extrato = {
    "MM Aparente": [mma],
    "Rm": [r_mistura],
    "Zm": [z_mistura]
}
extrato = pd.DataFrame(extrato)


vaz_massica, vaz_molar = corr_vazao_normal(
    vaz_norm = 4.54, 
    fator_z = z_mistura,
    mmolar_ap = mma,
    temp_crr = temp_gas_afSI,
)

df_vazoes = pd.DataFrame(
    [

        df_composicao.loc["% em Massa"] / 100 * vaz_massica,
        (df_composicao.loc["% em Massa"] / 100) * (mma /  df_composicao.loc["Massa Molar - (g/mol)"]) * vaz_molar
        
        
    ], 
    index=["vazao mássica individual", "vazao molar individual"]
)


df_vazoes, df_vaz_comb = vaz_combustao(df_vazoes.iloc[1,:].to_frame().transpose())

vgc = df_vaz_comb.iloc[0,[0,2]].sum() # CO2 e N2
o2livre = 0.1908 * vgc
print(
    "Volume dos gases de combustão," 
    f"em base seca: ",
    round(vgc,4),
    "mols/s."
)
print(
    "O2 livre: ",
    round(o2livre,4),
    "mols/s."
)
ar_teo = df_vazoes["O2"][0] # ar teórico

df_vazoes["N2"] += o2livre*3.76
df_vazoes["O2"] += o2livre

ex_ar = abs(1 - df_vazoes["O2"][0]/ar_teo)
ex_ar *= 100

df_entalpias = pd.DataFrame(
    {
        "CO": [-110530, None],
        "CH4": [-74850, None],
        "O2": [0, 10809],
        "N2": [0, 10763],
        "CO2": [-393520, 12148],
        "H2O": [-241820, 12331]
    },
    index = ["hf", "hT"]
)

temp_K, raz_equ, dataframe_b = temp_adiabatica(
    df_reagentes= df_vazoes,
    df_produtos= df_vaz_comb,
    ar_teorico= ar_teo,
    metodo="moran",
    nmax_iter=1e3
    
)