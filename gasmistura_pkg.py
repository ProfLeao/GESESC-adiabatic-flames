import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI 
def mm_aparente_mistura(
    compostos, mmolar, fracao, t_fracao, CNTP=True, TP=(273.153, 101325)
):
    """
        Calcula a massa molecular aparente, a constante de gás e o fator de 
        compressibilidade aparente de uma mistura de gases reais, 
        utilizando a lei de Amagat.

        Recebe uma lista contendo as strings dos gases componentes da mistura, 
        uma série de frações molares e uma string contendo o tipo da descrição 
        estequimétrica:
        · "fmass" - para fração mássica;
        · "fmol" - para fração molar.

        Por padrão a função trabalha nas condições normais de temperatura e 
        pressão 273,153 K e 101325 N/m². Esta condição pode ser modificada 
        alterando o parâmetro padrão CNTP para False e declarando, por meio de 
        uma par ordenado (list, tuple ...) os novos valores de temperatura e 
        pressão no sistema internacional de unidades.

        Retorna a massa molecular aparente da mistura `mma`, a constante 
        de gás da mistura `const_r_m` e o fator de compressibilidade da mistura
        `z_m`.
    """

    if t_fracao.lower() == "fmass":
        # Massa Molar Aparente
        fmass = fracao
        
        mma = 1 / np.sum(
            [
                fmass / PropsSI("MOLARMASS", 'T', TP[0], 'P', TP[1], i)\
                for i in compostos
            ]
        )
    elif t_fracao.lower() == "fmol":
        fmol = fracao
        
        mma = np.sum(
            [
                fmol * PropsSI("MOLARMASS", 'T', TP[0], 'P', TP[1], i)\
                for i in compostos
            ]
        )
    else:
        raise AttributeError(
            "Erro de parametrização ao tentar calcular a massa molecular ",
            "aparente da mistura na função 'mm_aparente_mistura'.\n",
            "Verifique a documentação!"
        )
    
    # Constante de gás da mistura
    const_r_m = 8.314462618/mma


    # Fator de compressibilidade Z da mistura
    if fmol.size:
        z_m = np.sum(
            [
                fmol * PropsSI('Z', 'T', TP[0], 'P', TP[1], i)\
                for i in compostos
            ]
        )
    return mma, const_r_m, z_m

def corr_vazao_normal(
    vaz_norm, 
    fator_z,  
    mmolar_ap,
    temp_crr,
    const_rm = 8.314462618,
    press_pad = 101325
):
    """
        Corrige vazões normais no SI para valores de vazão em T-p, atualmente, 
        usando o fator de compressibilidade como correção para a lei geral 
        dos gases ideias. Outras equações de estado podem ser implementadas 
        futuramente. 

        Parâmetros:
        · voz_norm - a vazão normal a ser corrigida em unidades SI;
        · factor_z - o fator de compressibilidade da mistura de gases ou do 
        gás;
        · mmolar_ap - a massa molar aparente dos componentes da mistura gasosa ou d
        o gás;
        · temp_crr - temperatura absoluta no estado de correção;
        · const_rm (opcional) - a constante de gás para a correção, se não for 
        for informada o valor para os gases ideiais é empregado 
        (8,314 462 618 J/mol.K);

        Retornos:
        Retorna um par de valores:
        1. A vazão mássica em kg/s;
        2. A vazão molar em mol/s.
    """

    vaz_mol = (press_pad * vaz_norm)/(fator_z * const_rm * temp_crr)
    vaz_mass = vaz_mol * mmolar_ap

    return vaz_mass, vaz_mol

def df_vazoes(dataframe):
    dict_vals = {}
    dict_vals["CO2 comb"] = dataframe.iloc[1,0] + dataframe.iloc[1,3]
    dict_vals["H2O"] = dataframe.iloc[1,1] + dataframe.iloc[1,2] + dataframe.iloc[1,3]
    dict_vals["Ar"] = ()
