#!/usr/bin/env python
# coding: utf-8

# ![IFMG.jpg](attachment:IFMG.jpg)

# <h1 style="text-align: center">Combustão dos Gases do Alto Forno na Fornalha da Caldeira</h1>
# 
# <h1 style="text-align: center">Indústria Siderúrgica</h1>
# 
# <p style="font-weight:bold; text-align: center">Mestrando: Moises Teles Fernandes</p>
# <p style="font-weight:bold; text-align: center">Orientador: Prof. Dr. Ricado C. Carpio</p>
# <p style="font-weight:bold; text-align: center">Co-orientador: Prof. Dr. Rogério José da Silva</p>
# <p style="font-weight:bold; text-align: center">Co-autor: Prof. Dr. Reginaldo Gonçalves Leão Junior</p>
# 
# ### Cuidados técnicos
# Antes de executar esse *notebook* certifique-se de estar em um ambiente cujas as seguintes dependências estejam satiasfeitas.
# 
# 1. Cantera
# 2. CoolProp
# 3. Matplotlib
# 4. Pandas
# 6. Champy
# 5. jupyter_contrib_nbextensions
# 
# **IMPORTANTE:** A compatibilidade deste notebook não foi testada em outros editores de notebooks com JupyterLab ou NBViewer.
# 
# Para maximizar a compatibilidade das funcionaldiades, recomendamos o uso do gestor de pacotes Anaconda. As especificações de instalação estão descritas abaixo.
# 
# #### Instalação do Cantera
# No *anaconda prompt* em ambiente `base` use:
# ```
# conda create --name ct-env --channel cantera cantera ipython matplotlib jupyter
# ```
# 
# após a instalação, feche o terminal e abra um novo e então digite a instrução abaixo para ativar o ambiente Cantera. 
# ```
# conda activate ct-env
# ```
# #### Instalação do CoolProp
# O CoolProp deve estar instalado no ambiente `ct-env` para que este notebook funcione corretamente. Caso não esteja, utilize o comando abaixo:
# ```
# conda install -c conda-forge coolprop 
# ```
# mais detalhes técnicos da instalação podem ser obtidos em (http://www.coolprop.org/).
# 
# #### Instalação do Pandas.
# A biblioteca Pandas pode ser instalada em seu ambiente `ct-env` por meio do comando:
# ```
# conda install pandas
# ```
# caso não esteja disponível nos canais padrão de seu ambiente utilize o comando abaixo para forçar a seleção do canal `conda-forge`:
# ```
# conda install -c conda-forge pandas 
# ```
# 
# #### Instalação do openpyxl
# Módulo adicional para a manimulação de planilhas construídas com Excel. Deve estar disponível em seu ambiente `ct-env`. Para instalá-la use:
# 
# ```
#  conda install -c anaconda openpyxl 
# ```
# 
# #### Instalação e configuração de complementos para o jupyter-notebook
# Utilize a seguinte sequência de comandos em seu ambiente `ct-env` para instalar e configurar os complementos necessários do Jupyter Notebook.
# ```
# conda install -c conda-forge jupyter_contrib_nbextensions
# jupyter contrib nbextension install --user
# jupyter nbextension enable python-markdown/main
# ```

# **(Cogeração em uma Indústria Siderúrgica).** Realizar um balanço térmico para um sistema de cogeração em uma indústria siderúgica Não integrada. (Utilizando o programa COOLPROP para o cálculo das entalpias)
# e simular a temperatura adiabática de chama como função da vazão de gases do autoforno via CANTERA (GRIMECH3.0).
# 

# ![gases.jpg](attachment:gases.jpg)

# In[1]:


import pandas as pd
from CoolProp.CoolProp import PropsSI 
import matplotlib.pyplot as plt
import numpy as np


# <h1 style="text-align: center">DADOS DE ENTRADA</h1>
# 

# ### Tabela de Composição dos Gases do Alto Forno (AF):

# In[2]:


df_composicao = pd.read_json("composicao.json", orient='split')
#display(df_composicao)


# In[3]:


df_composicao.loc[["Massa Molar - (g/mol)"]]


# ### Parametrização Física da Simulação:

# In[4]:


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


# ### **Dados Referentes ao combustível que entra na Caldeira:** 
# Tomando a *massa molecular aparente* da mistura como na forma dada pelas Equações *1.a* ou *1.b*;
# 
# $$
# \large M_{ma} =\frac{1}{ \#displaystyle \sum_{i=1}^{k}{\frac{f_{mas_i}}{M_i}}} \quad (1.a)
# $$
# 
# $$
# \large M_{ma} = \#displaystyle \sum_{i=1}^{k}{f_{mol_i} M_i}  \quad (1.b)
# $$
# 
# 
# onde $M_{ma}$ é a massa molecular aparente,$M_i$, $f_{mol_i}$ e $f_{mas_i}$ são, respectivamente, a massa molar, a fração molar e a fração mássica do iésimo composto na mistura, a constante de gás $R_m$ pode ser obtida conforme a Equação *2*. 
# 
# $$
# \large R_m = \frac{R_u}{M_{ma}} \quad (2)
# $$
# 
# Além disso, adimitindo que a relação entre $f_{mas_i}$ e a $f_{mol_i}$ se dá por meio da Equação *3*. 
# 
# $$
# \large f_{mas_i} = f_{mol_i} \frac{M_i}{M_{ma}} \quad (3)
# $$
# 
# Considerando que a lei de Amagat como uma boa aproxiamação para a descrição da mistura, então o fator de compressibilidade da mistura pode ser obtido pela Equação *4*, onde $f_{mol_i}$ é a fração molar de cada um dos componentes da mistura que podem ser obtidos em `df_composicao`. 
# 
# $$
# \large Z_m = \#displaystyle \sum^{k}_{i=1} {f_{mol_i} Z_i} \quad (4)
# $$
# 

# ### Correção de Valores e Parametrização Vetorial
# O módulo `gasmistura` contém funções que implementam o formulismo usado para a obtenção
# das características termodinâmicas da mistura de gases.

# In[5]:


from gasmistura_pkg import mm_aparente_mistura
help(mm_aparente_mistura)


# In[6]:


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
#display(extrato)


# Considerando a lei geral dos gases ideias corrigida pelo fator de compressibilidade como uma boa aproximação para o comportamento da mistura e derivando ambos os lados da igualdade obtem-se a forma da vazão molar da mistura em uma atmosfera conforme as Equações 5.a e 5.b.
# $$
# \large \frac{d(p V)}{dt} = \frac{(Z_m n R T)}{dt}
# $$
# 
# $$
# \large p \dot V = Z_m R T \dot n 
# $$
# 
# $$
# \large \dot n = \frac{p \dot V}{Z_m R T} \quad U[mol.s^{-1}] \quad (5.a) 
# $$
# 
# $$
# \large \dot m = \frac{p \dot V}{Z_m R T} M_{ma} \quad U[kg.s^{-1}] \quad (5.b) 
# $$
# 
# O mólduo `gasmistura` possui uma função chamada `corr_vazao_normal` que automatiza o processo de correção das vazões normais utilizando o formulismo acima. 

# In[7]:


from gasmistura_pkg import corr_vazao_normal
help(corr_vazao_normal)


# In[8]:


vaz_massica, vaz_molar = corr_vazao_normal(
    vaz_norm = 4.54, 
    fator_z = z_mistura,
    mmolar_ap = mma,
    temp_crr = temp_gas_afSI,
)

print(f"À temperatura de {temp_gas_afSI} K vazão mássica do gás é de {vaz_massica:.3f} kg/s e a vazão molar de {vaz_molar:.3f} mol/s.")


# <h1 style="text-align: center">EQUAÇÃO ESTEQUIOMÉTRICA</h1>
# 
# <h1 style="text-align: center">Equação de Balanço de Massa - EBM</h1>
# 
#  

# ### Escrevendo a equação estequiométrica (100\% de ar teórico).
# 
# 
# \begin{equation*}
# \large CO + H_2 + H_2O + CH_4 + N_2 + CO_2 (AF) + (a)\big(O_2 + 3,76 N_2\big)   \longrightarrow  ...
# \end{equation*}
# 
# \begin{equation*}
#  \large  ... \longrightarrow (X) CO_2 (Combustível) + CO_2 (AF) + (Y)H_2O + (Z) N_2   
# \end{equation*}

# As vazões molares e másicas individuais podem ser determinadas simplesmente pela derivada do número de mols $n_i$ e da massa $m_i$ no tempo conforme as equações 6 e 7, onde é possível constatar que o produto das frações (mássicas e molares) pelas vazões mássica e molar do gás podem determinar as grandezas de interesse para o balançao de massa. 
# 
# 
# $$
# \large \frac{d n_i}{dt} = \dot n_i = f_{mol_i} \dot n  \quad(6)
# $$
# 
# $$
# \large \frac{d m_i}{dt} = \dot m_i = f_{mas_i} \dot m  \quad(7)
# $$
# 
# Neste caso, como a composição do gás é dada em frações mássicas em `df_composicao` a Equação 3 pode ser usada para corrigir a Equação 6 na forma dada pela Equação 8. 
# 
# $$
# \large \dot n_i = f_{mass_i} \frac{M_{ma}}{M_i} \dot n \quad(8)
# $$
# 
# Com as Equações 7 e 8 definimos a dataframe `df_vazoes` que contém esses valores individuais.

# In[9]:


df_vazoes = pd.DataFrame(
    [

        df_composicao.loc["% em Massa"] / 100 * vaz_massica,
        (df_composicao.loc["% em Massa"] / 100) * (mma /  df_composicao.loc["Massa Molar - (g/mol)"]) * vaz_molar
        
        
    ], 
    index=["vazao mássica individual", "vazao molar individual"]
)
#df_vazoes


# ### Equação estequiométrica com seus coeficientes:
# 
# ### {{round(df_vazoes.iloc[1,0],4)}} $CO$ +  {{round(df_vazoes.iloc[1,1],4)}} $H_2$ + {{round(df_vazoes.iloc[1,2],4)}} $H_2 O$ + {{round(df_vazoes.iloc[1,3],4)}} $CH_4$ +  {{round(df_vazoes.iloc[1,4],4)}} $N_2$ +  {{round(df_vazoes.iloc[1,5],4)}}  $CO_2$  + a . $(O_2 + 3,76 N_2) \rightarrow$
# ### x . $CO_2$(Comb) + {{round(df_vazoes.iloc[1,5],4)}} $CO_2$ (AF) + y . $H_2O$ + z . $N_2$  
# 
# 
# Equação cuja a solução pode ser generalizada para:
# 
# $\large x = n° CO + n° CH_4$
# 
# $\large y = n° H_2 + n° H_2 O + 2 \cdot n° CH_4$
# 
# $\large a = \frac{n° CO}{2} + 2 \cdot n° CH_4 + \frac{n° H_2}{2}$
# 
# $\large z = n° N_2 + 3.76 \cdot a$
# 
# Um *DataFrame* para a solução dessas equações é devolvido como segundo retorno da função `vaz_combustao` que recebe como argumento `df_vazoes`, além disso a função devolve como primeiro valor a quantidade de ar (vazão) teórica para a combustão completa dos gases. 
# 
# A importação da função e sua documentação é mostrada abaixo. 

# In[10]:


from gasmistura_pkg import vaz_combustao
help(vaz_combustao)


# Deve-se notar que a documentação exclarece que o `dataframe` passado como argumento deve ser `um dataframe de uma linha com os dados de vazão dos gases do autoforno`.
# 
# A extração desses valores pode ser feita como:

# In[11]:


df_vazoes.iloc[1,:].to_frame().transpose()


# In[12]:


df_vazoes, df_vaz_comb = vaz_combustao(df_vazoes.iloc[1,:].to_frame().transpose())
#display(df_vazoes)
#display(df_vaz_comb)


# No código anterior o `dataframe` `df_vazoes`é extendido para a inserção do $O_2$ atmosférico teórico e a modificação do montante de $N_2$. Além disso, a variável `df_vaz_comb` contém o dataframe com a composição dos gases de combustão. 

# ### Equação estequiométrica completa com os respetivos coeficientes fica:
# 
# 
# ### {{round(df_vazoes.iloc[0,0],4)}} $CO$ +  {{round(df_vazoes.iloc[0,1],4)}} $H_2$ + {{round(df_vazoes.iloc[0,2],4)}} $H_2 O$ + {{round(df_vazoes.iloc[0,3],4)}} $CH_4$ +  {{round(df_vazoes.iloc[0,4],4)}} $N_2^{Ar+AF}$ +  {{round(df_vazoes.iloc[0,5],4)}}  $CO_2$  + {{round(df_vazoes.iloc[0,6],4)}}$O_2 \rightarrow$
# 
# ### {{round(df_vaz_comb.iloc[0,0],4)}}  $CO_2$(Comb) + {{round(df_vazoes.iloc[0,5],4)}} $CO_2$ (AF) + {{round(df_vaz_comb.iloc[0,1],4)}} $H_2O$ + {{round(df_vaz_comb.iloc[0,2],4)}} $N_2$  
# 

# <h1 style="text-align: center">CÁLCULO DE $O_2$ LIVRE (EM BASE SECA)</h1>

# ### BALANÇO DE MASSA E ENERGIA:  
# 
# 
# \begin{equation}
# \large \% O_2 livre =  \frac{O_2 livre}{V_{gc} + 4,76 O_2 livre} \quad (9)
# \end{equation}
# 
# ### Dados:
# 
# \begin{equation}
# \large \% O_2 livre =  10\% 
# \end{equation}
# 
# \begin{equation}
# \large O_2 livre =  0,1908 V_{gc} \quad (10)
# \end{equation}
# 
# onde:
# $V_{gc}$  é volume dos gases de combustão em base seca
# 

# ### Cálculo do Volume de Gases de Combustão ($V_{gc}$)
# Pode ser determinado pelo somatório da única linha de `df_vaz_comb` e o $O_2{livre}$ pela Equação 10.

# In[13]:


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
ar_teo


# Considerando a necessidade adicional de $O_2$ estabelecida na Equação 10, as quantidades estequimétricas correspondentes ao ar presentes em `df_vazoes` precisam ser acrescidas desses novos valores. 
# 
# Aqui atualiza-se as frações para contemplar a nova demanda de $O_2$, no entanto, antes, preserva-se a quantidade de ar teórico para a determinação posterior do excesso de ar.

# In[14]:


df_vazoes["N2"] += o2livre*3.76
df_vazoes["O2"] += o2livre
#display(df_vazoes)
ar_teo


# ### Valores de $O_2$ e $N_2$ para a equação de reação considerando $O_2 livre$
# 
#  ### $O_2$  = $O_2 livre$ + $O_2 {estequimétrico}$ = {{round(float(df_vazoes["O2"]),4)}}  mol/s
# 
# 
# ### $N_2$  = $3,76 \times O_2 livre + N_2 {estequiométrico}$ =  {{round(float(df_vazoes["N2"]),4)}}  mol/s
# 

# <h1 style="text-align: center">Equação de reação considerando - O2 livre</h1>
# 

# ### Equação de reação:
# 
# 
# ### {{round(float(df_vazoes["CO"]),4)}} $CO$ + {{round(float(df_vazoes["H2"]),4)}}  $H_2$ + {{round(float(df_vazoes["H2O"]),4)}}  $H_2O$ + {{round(float(df_vazoes["CH4"]),4)}}  $CH_4$ + {{round(float(df_vazoes["N2"]),4)}}  $N_2^{Ar + AF}$ + {{round(float(df_vazoes["CO2"]),4)}}  $CO_2$ + {{round(float(df_vazoes["O2"]),4)}} $O_2 \longrightarrow$ 
# 
# ### {{round(float(df_vaz_comb["CO2"]),4)}}  $CO_2^{comb + AF}$ + {{round(float(df_vaz_comb["H2O"]),4)}} $H_2 O$ + {{round(float(df_vaz_comb["N2"]),4)}}$N_2$ 

# <h1 style="text-align: center">Cálculo do Excesso de Ar</h1>
# 

# ### Equação para calcular o excesso de Ar:
# 
# 
# \begin{equation*}
# \large \% execesso-de-ar =   \frac{Ar_{real}}{Ar_{teórico}} 
# \end{equation*}

# In[15]:


ex_ar = abs(1 - df_vazoes["O2"][0]/ar_teo)
ex_ar *= 100
ex_ar


# ### Excesso de Ar será:
# 
# 
# ### \%_{execesso-de-ar} =   {{round(ex_ar,3)}} % 
# 
# 

# <h1 style="text-align: center">Análise Energético - Equação de Balanço de Energia - EBE</h1>
# <h1 style="text-align: center">Temperatura Adiabática da Chama dentro da Fornalha</h1>
# 
# Considerando uma situação ideal na qual a perca energética durante a combustão se dê exclusivamente através da energia dos produtos da combustão, pode-se considerar como verdadeira a Equação 11.

# ### EBE:
# 
# \begin{equation*}
# \large   \sum_{reagentes}  {\dot n_{i}[ \bar{h}^{o}_f + \bar{h}_{T} - \bar{h}_{Tref}]_{r} }  = \sum_{produtos}  {\dot n_{i}[ \bar{h}^{o}_f + \bar{h}_{T} - \bar{h}_{Tref}]_{p} } \quad (11)
# \end{equation*}
# 
# Além disto, se tomarmos como nula a variação da entalpia padrão ($\bar{h}_{T} - \bar{h}_{Tref}$) dos reagentes, o balanço energético assume a forma dada pela Equação 12 onde $\bar{h}$ é a função de entalpia para a temperatura absoluta $T$.
# 
# \begin{equation*}
# \large \sum_{produtos}  {\dot n_{i}[ \bar{h}^{o}_f + \bar{h}(T)_{pad} - \bar{h}(T)_{ref}]_{p} } - \sum_{reagentes}  {\dot n_{i} \bar{h}^{o}_f} = 0\quad (12)
# \end{equation*}
# 
# Nesta forma a determinação resume-se a um problema de determinação das raízes da Equação 12 em um determinado intervalo. 

# ### Os valores de calor de formação $\bar{h}^{o}_f$ e valores da entalpia a temperatura de referência $\bar{h}_{Tref}$. Tabela A-23 - Moran.
# 
# 
# |substância||$\bar{h}^{o}_f$ (kJ/kmol)||$\bar{h}(T){ref}$ (kJ/kmol)|
# |:-------:||:----:||:----:|
# |$CO$||-110530||-|
# |$CH_4$||-74850||-|
# |$O_2$||0||10809|
# |$N_2$||0||10763|
# |$CO_2$||-393520||12148|
# |$H_2O$||-241820||12331|
# 
# A implementação de seus valores como um *dataframe* (`df_entalpias`) é mostrado abaixo. Note que para manutenção do padrão de referênciamento dos dados, a forma de `df_entalpias` e a transposta da tabela **A-23**.

# In[16]:


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
#display(df_entalpias)


# ### O método iterativo de Moran e Shapiro:
# 
# 
# O método iterativo proposto por Moran e Shapiro (ano?) consiste em uma heurística simples de busca das soluções por método *brute force*. Tecnicamente ele itera de forma sucessiva verificando a igualdade entre as parcelas da Equação 12 até encontrar uma aproximação razoável da solução.
# 
# A função `temp_adiabatica` do módulo `gasmistura_pkg` disponibiliza esse recurso pemo método `moran`. Veja documentação abaixo. 

# In[17]:


from gasmistura_pkg import temp_adiabatica
help(temp_adiabatica)


# Na célula abaixo, passo reagentes é passado os dataframes de reagentes e produtos e o ar teórico determinado anteriormente, além disso o método `moran` é selecionado um um total de no máximo 1000 iterações é selecionado.

# In[18]:


temp_K, raz_equ = temp_adiabatica(
    df_reagentes= df_vazoes,
    df_produtos= df_vaz_comb,
    ar_teorico= ar_teo,
    metodo="moran",
    nmax_iter=1e3
    
)


# ### Substituindo os valores da tabela, aplicando a Eq.2 a cada elemento da Eq. 1:
# 
# 
# \begin{equation*}
# 0,0446\times (-110530) + 0,0036\times (-74850) =  0,0482\times[-393520 - \#displaystyle{\left ( \bar{h}_{CO_2}(T) - 12148\right)}] + ...            
# \end{equation*}
# 
# \begin{equation*}
# ...+ 0,1921\times [-241820 - \#displaystyle{\left ( \bar{h}_{H_2O}(T) - 12331\right)}] + 2,0264 \times[0 + \#displaystyle{\left ( \bar{h}_{N_2}(T) - 10763\right)}] + 0,2326 \times[0 + \#displaystyle{\left ( \bar{h}_{O_2}(T) - 10809\right)}]         
# \end{equation*}

# In[ ]:


reagente = 0.0446*(hfCO) + 0.1540*(hfCH4) 
produtos = 0.0482*(hfCO2 - TrCO2) + 0.1921*(hfH2O - TrH2O) + 2.0264*(-TrN2) + 0.2326*(-TrO2)
RP = reagente - produtos
print("O valor dos reagentes = ","{0:.2f}".format(reagente),'kJ/kmol')
print ("O valor dos produtos = %.3f kJ/kmol" %produtos)
print ("O valor dos reagentes menos produtos  será = %.3f kJ/kmol" %RP)


# ### Podendo ser também:
# 
# \begin{equation*}
# \large  0,0482 \times \bar{h}_{CO_2}(T) +0,1921 \times \bar{h}_{H_2O}(T) + 2,0264 \times \bar{h}_{N_2}(T) + 0,2326 \times \bar{h}_{O_2}(T)= 76243,383 (kJ/kmol) 
# \end{equation*}

# # Iteração 1: Assumindo T = 1000 K
# 
# ### Substituindo os valores do lado esquerdo para 1000 K.  E verificando os valores da Tabela A-23 (Moran), obtemos:

# \begin{equation*}
# \large  0,0482 \times \bar{h}_{CO_2}(1000) +0,1921 \times \bar{h}_{H_2O}(1000) + 2,0264 \times \bar{h}_{N_2}(1000) + 0,2326 \times \bar{h}_{O_2}(1000)= 76243,383 
# \end{equation*}
# 
# $$\;$$
# $$\;$$
# 
# \begin{equation*}
# \large  0,0482 \times \bar{h}_{CO_2}(42769) + 0,1921 \times \bar{h}_{H_2O}(35882) + 2,0264 \times \bar{h}_{N_2}(30129) + ... 
# \end{equation*}
# 
# \begin{equation*}
# \large  ... + 0,2326 \times \bar{h}_{O_2}(31389)= 76243,383 
# \end{equation*}
# 
# $$\;$$
# 
# \begin{equation*}
# \large  77308,89 > 76243,383
# \end{equation*}

# In[ ]:


LE1 = 0.0482*42769 + 0.1921*35882 + 2.0264*30129 + 0.2326*31389
print("O valor do lado Esquerdo = ","{0:.2f}".format(LE1),'kJ/kmol')


# # Iteração 2: Assumindo T = 800 K
# 
# ### Substituindo os valores do lado esquerdo para 800 K.  E verificando os valores da Tabela A-23 (Moran), obtemos:

# \begin{equation*}
# \large  0,0482 \times \bar{h}_{CO_2}(800) +0,1921 \times \bar{h}_{H_2O}(800) + 2,0264 \times \bar{h}_{N_2}(800) + 0,2326 \times \bar{h}_{O_2}(800)= 76243,383 
# \end{equation*}
# 
# $$\;$$
# $$\;$$
# 
# \begin{equation*}
# \large  0,0482 \times \bar{h}_{CO_2}(32179) + 0,1921 \times \bar{h}_{H_2O}(27896) + 2,0264 \times \bar{h}_{N_2}(23714) + ... 
# \end{equation*}
# 
# \begin{equation*}
# \large  ... + 0,2326 \times \bar{h}_{O_2}(24523)= 76243,383 
# \end{equation*}
# $$\;$$
# \begin{equation*}
# \large  60667,95 < 76243,383
# \end{equation*}

# In[ ]:


LE2 = 0.0482*32179 + 0.1921*27896 + 2.0264*23714 + 0.2326*24523
print("O valor do lado Esquerdo = ","{0:.2f}".format(LE2),'kJ/kmol')


# # Interpolando entre T = 1000 K e  T = 800 K: 
# 
# 
# |$h$||T(K)|
# |:--:||:--:|
# |60667,95||800|
# |76243,383||$T$|
# |77308,89||1000|

# # Obtemos a temperatura Adiabática de Chama.
# 
# #                                             T = 987,19 K
