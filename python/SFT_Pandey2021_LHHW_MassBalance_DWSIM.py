# =========================================================
# SCRIPT PYTHON IMPLEMENTADO NO DWSIM PARA A SÍNTESE DE
# FISCHER-TROPSCH (SFT) EM REATOR DE LEITO FIXO CATALÍTICO
# =========================================================
#
# Finalidade:
# Este script modela a etapa de SFT em um reator PFR
# (Plug Flow Reactor) pseudo-homogêneo, usando como variável
# independente a massa acumulada de catalisador (W).
#
# O modelo considera:
# 1) consumo de CO pela reação de Fischer-Tropsch;
# 2) reação de deslocamento gás-água (WGS);
# 3) distribuição de produtos por Anderson-Schulz-Flory (ASF);
# 4) produção explícita de C5–C19;
# 5) agrupamento da fração C20+ em um lump C24.
#
# O DWSIM resolve as propriedades termodinâmicas do stream,
# enquanto este script calcula as taxas cinéticas e atualiza
# a composição do fluxo de saída.
# =========================================================

import math
import System
from System import Array, Double

# =========================================================
# 1. LEITURA DOS STREAMS CONECTADOS AO BLOCO PYTHON
# =========================================================
# ims1 = corrente de entrada (SYNGAS_FEED)
# oms1 = corrente de saída (FT_OUT)
feed = ims1
prod = oms1

# =========================================================
# 2. ÍNDICES DOS COMPONENTES NO DWSIM
# =========================================================
# Estes índices correspondem à ordem interna dos componentes
# cadastrados no flowsheet do DWSIM.
#
# Espécies reagentes e inorgânicas:
iCO  = 0   # Carbon monoxide
iH2  = 1   # Hydrogen
iH2O = 2   # Water
iCO2 = 3   # Carbon dioxide

# Hidrocarbonetos na faixa gasolina sintética:
iC5  = 11  # N-pentane
iC6  = 12  # N-hexane
iC7  = 13  # N-heptane
iC8  = 14  # N-octane
iC9  = 15  # N-nonane
iC10 = 16  # N-decane
iC11 = 17  # N-undecane

iC12 = 18  # N-dodecane
iC13 = 19  # N-tridecane
iC14 = 20  # N-tetradecane
iC15 = 21  # N-pentadecane
iC16 = 22  # N-hexadecane
iC17 = 23  # N-hexadecane
iC18 = 24  # N-hexadecane
iC19 = 25  # N-hexadecane

# Lump pesado representando C20+: 
iC24 = 27  # N-icosano

# Lista dos hidrocarbonetos explicitamente contabilizados
# (C5–C19):
hc_idx = [
    (5, iC5), (6, iC6), (7, iC7),
    (8, iC8), (9, iC9), (10, iC10), (11, iC11),
    (12, iC12), (13, iC13), (14, iC14), (15, iC15),
    (16, iC16), (17, iC17), (18, iC18), (19, iC19)  ]

# Definição do lump pesado: 24
heavy_n, heavy_idx = (24, iC24)

# =========================================================
# 3. PARÂMETROS DO REATOR
# =========================================================
# Wcat_total = massa total de catalisador no leito (kg).
# Esse valor foi variado nas simulações paramétricas.
#
# Nsteps = número de passos de discretização do PFR.
# O reator é dividido numericamente em Nsteps pequenos
# incrementos de massa catalítica.
#
# dW = incremento de massa catalítica em cada passo.
Wcat_total = 15.0   # kg de catalisador (ajustável conforme o caso)
Nsteps = 200        # discretização numérica do PFR
dW = Wcat_total / float(Nsteps)

# =========================================================
# 4. CONSTANTES E FUNÇÕES CINÉTICAS
# =========================================================
# O modelo usa:
# - constante universal dos gases R
# - equação de Arrhenius deslocada ("shifted Arrhenius")
# - taxa FT do tipo LHHW
# - correlação para alpha (probabilidade de crescimento)
# - taxa de WGS
#

R = 8.314  # J/mol/K

def arrhenius_shifted(kref, Ea_Jmol, T):
    """
    Calcula a constante cinética por uma forma de Arrhenius
    deslocada em relação a 483 K.

    k = kref * exp[-Ea/R * (1/T - 1/483)]

    Parâmetros:
    kref     = constante de referência
    Ea_Jmol  = energia de ativação (J/mol)
    T        = temperatura (K)
    """
    return kref * math.exp(-(Ea_Jmol / R) * (1.0 / T - 1.0 / 483.0))

# ---------------------------------------------------------
# 4.1. Taxa de Fischer-Tropsch (forma LHHW)
# ---------------------------------------------------------
# Constantes adotadas no script:
kref = 7.05       # constante de referência
Ea_FT = 92000.0   # energia de ativação (J/mol)
a0 = 12.0         # termo de adsorção associado ao CO
b  = 1.10         # termo associado a H2
f0 = 1.25         # termo associado a H2O
kpH2O = 0.1       # fator corretivo envolvendo água

def rFT_COcons(pCO, pH2, pH2O, T):
    """
    Taxa de consumo de CO pela reação de Fischer-Tropsch.

    A forma funcional é do tipo LHHW:
    - numerador: termo cinético principal
    - denominador: competição por sítios ativos

    Entrada:
    pCO, pH2, pH2O = pressões parciais em MPa
    T              = temperatura em K

    Saída:
    taxa na base do paper (depois convertida para mol/(kg*s))
    """
    k = arrhenius_shifted(kref, Ea_FT, T)
    num = k * pCO * math.sqrt(max(pH2, 0.0)) * (1.0 + kpH2O * pH2O)
    den = (1.0 + a0 * pCO + b * math.sqrt(max(pH2, 0.0)) + f0 * pH2O) ** 2
    val = num / den
    return max(0.0, val)

# ---------------------------------------------------------
# 4.2. Probabilidade de crescimento de cadeia (alpha)
# ---------------------------------------------------------
# Constantes adotadas:
kalpha_ref = 0.118
Ea_alpha = 4770.0
zexp = -0.170
yexp = -0.095
 
def alpha_growth(pCO, pH2O, T):
    """
    Calcula alpha, parâmetro de distribuição ASF.

    Alpha representa a probabilidade de crescimento de cadeia.
    Quanto maior alpha, maior a tendência de formar compostos
    de maior número de carbonos.

    O valor é limitado numericamente entre 0,05 e 0,95
    para evitar instabilidades.
    """
    kalphaT = arrhenius_shifted(kalpha_ref, Ea_alpha, T)
    pH2Oeff = max(pH2O, 1e-12)
    denom = (pCO ** zexp) * (pH2Oeff ** yexp)
    a = 1.0 / (1.0 + kalphaT * (1.0 / denom))
    return min(max(a, 0.05), 0.95)

# ---------------------------------------------------------
# 4.3. Reação de deslocamento gás-água (WGS)
# ---------------------------------------------------------
kWGS = 119.1

def rWGS(pCO, pH2O, pCO2, pH2, T):
    """
    Taxa da reação WGS:
    CO + H2O <-> CO2 + H2

    Inclui dependência de equilíbrio por Keq. J. M. MOE 1962
    """
    kco2 = kWGS * math.exp(-47400.0 / (R * T))
    Keq  = math.exp(4557.8 / T - 4.33)
    return kco2 * (pCO * pH2O - (pCO2 * pH2) / Keq)

# =========================================================
# 5. LEITURA DO STREAM DE ENTRADA
# =========================================================
# A corrente de saída é inicialmente copiada da entrada,
# para preservar T, P e demais propriedades antes da atualização.
prod.Clear()
prod.Assign(feed)

# Temperatura e pressão da corrente de entrada:
T = feed.GetTemperature()          # K
P_MPa = feed.GetPressure() / 1e6   # MPa

# z = composição global da corrente de entrada
# Ftot = vazão molar total da corrente (mol/s)
z = list(feed.GetOverallComposition())
Ftot = feed.GetMolarFlow()

# Se não houver fluxo, apenas calcula a saída e encerra.
if Ftot <= 0.0:
    prod.Calculate()

else:
    # =====================================================
    # 6. INICIALIZAÇÃO DAS VAZÕES MOLARES POR COMPONENTE
    # =====================================================
    # Fi guarda a vazão molar de cada componente (mol/s).
    Fi = [Ftot * zi for zi in z]

    # Conversão de unidades:
    # Assumiu-se que a taxa original está em kmol/(kgcat*h).
    # Portanto, converte-se para mol/(kgcat*s).
    conv = 1000.0 / 3600.0

    # =====================================================
    # 7. INTEGRAÇÃO DO PFR AO LONGO DA MASSA DE CATALISADOR
    # =====================================================
    # O reator é tratado como 200 pequenos reatores diferenciais
    # em série, cada um contendo uma massa dW de catalisador.
    for _ in range(Nsteps):

        # Vazão total instantânea em cada passo
        Ftot_now = sum(Fi)
        if Ftot_now <= 1e-18:
            break

        # -------------------------------------------------
        # 7.1. Cálculo das frações molares atuais
        # -------------------------------------------------
        yCO  = Fi[iCO]  / Ftot_now
        yH2  = Fi[iH2]  / Ftot_now
        yH2O = Fi[iH2O] / Ftot_now
        yCO2 = Fi[iCO2] / Ftot_now

        # -------------------------------------------------
        # 7.2. Cálculo das pressões parciais (MPa)
        # -------------------------------------------------
        pCO  = yCO  * P_MPa
        pH2  = yH2  * P_MPa
        pH2O = yH2O * P_MPa
        pCO2 = yCO2 * P_MPa

        # -------------------------------------------------
        # 7.3. Taxas reacionais
        # -------------------------------------------------
        rft  = rFT_COcons(pCO, pH2, pH2O, T) * conv
        rwgs = rWGS(pCO, pH2O, pCO2, pH2, T) * conv

        # -------------------------------------------------
        # 7.4. Distribuição de produtos por ASF
        # -------------------------------------------------
        a = alpha_growth(pCO, pH2O, T)

        # Pesos ASF para C5–C19
        w = {}
        for n, idx in hc_idx:
            w[n] = (1.0 - a) * (a ** (n - 1))

        # Soma total de carbono da distribuição ASF ideal
        carbon_total = 1.0 / (1.0 - a)

        # Carbono explicitamente alocado 
        carbon_hc = 0.0
        for n, idx in hc_idx:
            carbon_hc += n * w[n]

        # Carbono remanescente é agrupado na cauda C20+
        carbon_tail = max(0.0, carbon_total - carbon_hc)

        # -------------------------------------------------
        # 7.5. Fechamento de carbono
        # -------------------------------------------------
        # O consumo de CO na FT é redistribuído entre C5–C19
        # e o lump pesado C16.
        Rn = {}
        for n, idx in hc_idx:
            Rn[n] = (w[n] / carbon_total) * rft

        # Lump pesado C24 representando a fração C20+
        R24 = (carbon_tail / carbon_total) * rft / float(heavy_n)

        # -------------------------------------------------
        # 7.6. Consumo de H2 pela estequiometria de parafinas
        # -------------------------------------------------
        # n CO + (2n+1) H2 -> CnH(2n+2) + n H2O
        rH2_FT = 0.0
        for n, idx in hc_idx:
            rH2_FT += (2 * n + 1) * Rn[n]
        rH2_FT += (2 * heavy_n + 1) * R24

        # -------------------------------------------------
        # 7.7. Atualizações diferenciais em cada passo dW
        # -------------------------------------------------
        # Contribuições da FT
        dCO_FT  = -rft * dW
        dH2_FT  = -rH2_FT * dW
        dH2O_FT = +rft * dW

        # Contribuições da WGS
        dCO_WGS  = -rwgs * dW
        dH2O_WGS = -rwgs * dW
        dCO2_WGS = +rwgs * dW
        dH2_WGS  = +rwgs * dW

        # -------------------------------------------------
        # 7.8. Atualização dos reagentes/inorgânicos
        # -------------------------------------------------
        Fi[iCO]  = max(Fi[iCO]  + dCO_FT  + dCO_WGS,  0.0)
        Fi[iH2]  = max(Fi[iH2]  + dH2_FT  + dH2_WGS,  0.0)
        Fi[iH2O] = max(Fi[iH2O] + dH2O_FT + dH2O_WGS, 0.0)
        Fi[iCO2] = max(Fi[iCO2] + dCO2_WGS,           0.0)

        # -------------------------------------------------
        # 7.9. Atualização dos produtos C5–C19
        # -------------------------------------------------
        for n, idx in hc_idx:
            Fi[idx] = max(Fi[idx] + Rn[n] * dW, 0.0)

        # Atualização do lump pesado C24
        Fi[heavy_idx] = max(Fi[heavy_idx] + R24 * dW, 0.0)

    # =====================================================
    # 8. CÁLCULO DAS MÉTRICAS DE DESEMPENHO
    # =====================================================
    # Produção de hidrocarbonetos individuais:
    Fhc = 0.0
    for n, idx in hc_idx:
        Fhc += Fi[idx]

    # Produção de lump pesado:
    Fheavy = Fi[heavy_idx]

    # Conversão de CO:
    FCO_in  = Ftot * z[iCO]
    FCO_out = Fi[iCO]
    XCO = 0.0
    if FCO_in > 1e-12:
        XCO = (FCO_in - FCO_out) / FCO_in

    # Exibe resultados resumidos no log do DWSIM
    Flowsheet.ShowMessage(
        "RESULTADOS: Hidrocarbonetos(C5-C19)= %.6g mol/s | C24= %.6g mol/s | X_CO= %.3f"
        % (Fhc, Fheavy, XCO),
        0
    )

    # =====================================================
    # 9. ESCRITA DA CORRENTE DE SAÍDA
    # =====================================================
    Ftot_out = sum(Fi)

    if Ftot_out <= 0.0:
        prod.Calculate()
    else:
        # Nova composição molar global
        znew = [fi / Ftot_out for fi in Fi]

        # Atualiza vazão molar total e composição no DWSIM
        prod.SetMolarFlow(Ftot_out)
        prod.SetOverallComposition(Array[Double](znew))
        prod.Calculate()
