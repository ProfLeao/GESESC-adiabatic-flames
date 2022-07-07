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


def vaz_combustao(dataframe):
    """
        Realiza o balaço de massa da combustão dos gases de auto-forno.

        Parâmetros: 
        · um dataframe de uma linha com os dados de vazão dos gases do 
        autoforno. Na versão atual a composição deve ser limitada à:
        CO, H2, H2O, CH4, N2, CO2 e Ar. O Ar pode ser substituído por O2 e N2
        atmosféricos na forma n(O2 + 3.72 N2) quando conveniente, ou mesmo 
        omitido, neste último caso, a quantidade de ar é calculada 
        automaticamente.

        Retornos:
        · o dataframe dos gases de combustão com a demanda de ar-teórico caso 
        este não conste no dataframe original passado como argumento
        · um segundo dataframe com os gases de combustão já balanceados em
        massa. 
    """
    compostos = [val.lower() for val in dataframe.columns]
    dfcomb = pd.DataFrame()

    dfcomb.loc["vazao molar individual","CO2"] =\
        dataframe.loc["vazao molar individual","CO"] +\
        dataframe.loc["vazao molar individual","CH4"]

    dfcomb.loc["vazao molar individual","H2O"] =\
        dataframe.loc["vazao molar individual","H2"] +\
        dataframe.loc["vazao molar individual","H2O"] +\
        2 * dataframe.loc["vazao molar individual","CH4"]
    
    if "ar" in compostos:
        dfcomb.loc["vazao molar individual","N2"] =\
            dataframe.loc["vazao molar individual","N2"] +\
            3.72 * dataframe.loc["vazao molar individual","Ar"]

    elif "o2" in compostos:
        dfcomb.loc["vazao molar individual","N2"] =\
            dataframe.loc["vazao molar individual","N2"] +\
            3.72 * dataframe.loc["vazao molar individual","O2"]

    else:
        dataframe.loc["vazao molar individual","O2"] =\
            dataframe.loc["vazao molar individual","CO"]/2 +\
            dataframe.loc["vazao molar individual","H2"]/2 +\
            dataframe.loc["vazao molar individual","CH4"]*2

        dfcomb.loc["vazao molar individual","N2"] =\
            dataframe.loc["vazao molar individual","N2"] +\
            3.72 * dataframe.loc["vazao molar individual","O2"]

        dataframe.loc["vazao molar individual","N2"] +=\
            3.72 * dataframe.loc["vazao molar individual","O2"]
    
    return dataframe, dfcomb

def temp_adiabatica(
    df_reagentes, 
    df_produtos, 
    metodo,
    guess = 1000, 
    salva_dados = True,
    nmax_iter = 1e5):
    """
        Determina a temperatura adiabática da chama para um conjunto de 
        reagentes e produtos. Além disso, exibe no terminal um log do estado
        final alcançado pelo algoritmo característico do método utilizado. 
        Se o parâmetro opcional salva_dados for passado como True, um arquivo 
        .csv é criado no diretório de execução com os dados evolutivos do modelo
        para cada iteração realizada, assim como um arquivo de texto plano 
        .txt com o log exibido no terminal. 

        Parâmetros:
        · df_reagentes - DataFrame com colunas nomeadas segundo a fórmula 
        química dos compostos presentes no combustível e uma linha contendo 
        as vazões molares ou frações desses compostos.
        
        · df_prdutos - DataFrame com colunas nomeadas segundo a fórmula 
        química dos compostos presentes nos gases de combustão e uma linha 
        contendo as vazões molares ou frações desses compostos.

        · metodo - Strig com o método utilizado para a determinação da
        temperatura adiabática. 

            i) "moran" - método iterativo brute force proposto por moran e 
            shapiro.
            
            ii) "bisect" - método da bissecção
            
            iii)  "brentq"  - método clássico de Brent veja (Brent, R. P., 
            Algorithms for Minimization Without Derivatives. Englewood Cliffs, 
            NJ: Prentice-Hall, 1973. Ch. 3-4.)
            
            iv) "brenth" - método de Brent com extrapolação hiperbólica veja 
            (Bus, J. C. P., Dekker, T. J., “Two Efficient Algorithms with 
            Guaranteed Convergence for Finding a Zero of a Function”, ACM 
            Transactions on Mathematical Software, Vol. 1, Issue 4, Dec. 1975, 
            pp. 330-345. Section 3: “Algorithm M”. DOI:10.1145/355656.355659)

            v) "ridder" - método de Ridder

            vi) "toms748" - algorítimo 748 de Alefeld, Potro and Shi veja 
            (APS1995 Alefeld, G. E. and Potra, F. A. and Shi, Yixun, Algorithm 
            748: Enclosing Zeros of Continuous Functions, ACM Trans. Math. 
            Softw. Volume 221(1995) doi = {10.1145/210089.210111})

            vii) "newton" - método de Newton. 

            viii) "secant" - método das secantes.

            ix) "halley" - método de Halley veja (Boyd, J. P. (2013). "Finding 
            the Zeros of a Univariate Equation: Proxy Rootfinders, Chebyshev 
            Interpolation, and the Companion Matrix". SIAM Review. 55 (2): 
            375–396. doi:10.1137/110838297)
        
        · [guess] opcional - lista com sugestão de temperatura(s) para início 
        dos algoritmos. Se apenas um valor for passado, as buscas são realizadas 
        em torno do valor. Se dois valores forem passados, então entre os 
        os valores, inclusivamente.

        · [salva_dados] opcional - booleano; quando True, os dados da busca 
        são salvos no diretório de execução. A saber, um CSV com os valores 
        assumidos pelo algoritmo em cada iteração e um log com o estado final 
        do algoritmo. 

        · [nmax_iter] opcional - inteiro com o número máximo de iterações 
        permitidas. 

        Retorno:
        · float - temperatura adiabática da chama em Kelvins.         

        Para mais informações consultar:
        https://rb.gy/fdcsqf
    """
    
    # Obtem entalpias de formação do NIST
    try:
        url = "https://data.nist.gov/od/ds/mds2-2124/NBS_Tables%20Library.xlsx"
        df_form_enthalpies = pd.read_excel(url)
        df_form_enthalpies = df_form_enthalpies.iloc[:,[0,7]]
        df_form_enthalpies.columns = ["compostos", "entform"]
    except:
        print(
            "Impossível conectar ao Data NIST para extração das entalpias de",
            " formação."
        )

    dados = pd.DataFrame()
    try:
        if metodo.lower() == "moran":

            reagentes = df_reagentes.columns
            produtos = df_produtos.columns

            # consulta ao banco de dados DF o somatório das entalpuias de 
            # formação dos reagente 
            reag_ent_form =   df_form_enthalpies.query(
                f"compostos in {reagentes}"
            )["entform"].map(float).sum()

            # Iterações de busca:

            min_del_temp = 1 # K - menor temperatura a ser incrementada na busca
            del_temp = guess/2   # K - incremento inicial da temperatura
            temp = guess
            n_iter = 0
            while del_temp > min_del_temp:
                n_iter += 1
                prods_ental = np.sum(
                    [
                        float(df_form_enthalpies.loc[c,entform]) +\
                        df_produtos.loc["vazao molar individual",c] *\
                            df_produtos.PropsSI(
                                "HMOLAR", 
                                'T', 
                                temp, 
                                'P', 
                                101325, 
                                c
                            ) -\
                            df_produtos.PropsSI(
                                "HMOLAR",
                                'T', 
                                298, 
                                'P', 
                                101325, 
                                c
                            )
                            for c in produtos
                        ]
                )

                # Armazena os dados apenas se salva_dados == True
                if salva_dados:
                    df = pd.DataFrame(
                        [[n_iter, temp, prods_ental]], 
                        columns=["iteracoes", "temperaturaK", "sumentalpia"]
                    )
                    dados = pd.concat([dados, df])
                
                result = round(prods_ental, 2)
                if round(prods_ental, 2) == round(reag_ent_form, 2):
                    if salva_dados:
                        dados.to_csv("dados.csv", sep=",")
                    return result
                elif del_temp<=1:
                    if salva_dados:
                        dados.to_csv("dados.csv", sep=",")
                    return result
                elif prods_ental > reag_ent_form:
                    temp -= del_temp
                    del_temp *= .9
                elif prods_ental < reag_ent_form:
                    temp += del_temp
                    del_temp *= .9
                elif n_iter > nmax_iter:
                    if salva_dados:
                        dados.to_csv("dados.csv", sep=",")
                    raise Exception("Número máximo de iterações excedido.")
                else:
                    raise Exception("Erro desconhecido.")
    
    except:
        return None





