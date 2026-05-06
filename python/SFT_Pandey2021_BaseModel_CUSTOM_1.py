# ==============================================================
# APÊNDICE X – SCRIPT PYTHON IMPLEMENTADO NO DWSIM PARA A
# MODELAGEM DA SÍNTESE DE FISCHER–TROPSCH EM REATOR PFR
# PSEUDO-HOMOGÊNEO FORMULADO EM FUNÇÃO DA MASSA DE CATALISADOR
# ==============================================================
#
# Descrição geral:
# Este script foi implementado em uma unidade personalizada do tipo
# Python Script Unit Operation no DWSIM, com a finalidade de modelar
# a etapa de síntese de Fischer–Tropsch (SFT) em regime permanente.
#
# O reator é tratado como um PFR (Plug Flow Reactor) pseudo-homogêneo,
# discretizado ao longo da massa de catalisador, W, e não ao longo do
# comprimento geométrico do reator.
#
# O algoritmo realiza:
# 1) leitura da corrente de entrada fornecida pelo DWSIM;
# 2) cálculo das taxas da reação principal de Fischer–Tropsch (FT);
# 3) cálculo da reação de deslocamento gás–água (WGS);
# 4) cálculo do parâmetro alfa da distribuição Anderson–Schulz–Flory;
# 5) distribuição dos produtos entre C5–C19 e um lump pesado C24;
# 6) atualização das vazões molares em 200 incrementos de massa
#    catalítica por método de Euler explícito;
# 7) escrita da corrente final de saída no DWSIM.
#
# O DWSIM é responsável pelo fechamento termodinâmico da corrente de
# saída, enquanto o presente script é responsável pelo balanço molar
# reacional e pela cinética.
# ==============================================================

import clr
clr.AddReference("System")

# --------------------------------------------------------------
# 1. DEFINIÇÃO DOS STREAMS DE ENTRADA E SAÍDA NO DWSIM
# --------------------------------------------------------------
# ims1 = corrente material de entrada conectada à unidade Python
# oms1 = corrente material de saída conectada à unidade Python
feed = ims1
prod = oms1

# --------------------------------------------------------------
# 2. ÍNDICES DOS COMPONENTES NA ORDEM INTERNA DO DWSIM
# --------------------------------------------------------------
# Esses índices devem corresponder exatamente à ordem dos compostos
# cadastrados no fluxograma do DWSIM.

# Componentes gasosos reativos
iCO  = 0
iH2  = 1
iH2O = 2
iCO2 = 3

# Hidrocarbonetos explicitamente representados
iC5  = 11
iC6  = 12
iC7  = 13
iC8  = 14
iC9  = 15
iC10 = 16
iC11 = 17
iC12 = 18
iC13 = 20
iC14 = 21
iC15 = 22
iC16 = 19
iC17 = 23
iC18 = 24
iC19 = 25

# Lumps pesados
iC20 = 27
iC24 = 26

# Lista de componentes tratados explicitamente pela distribuição ASF
# Cada tupla contém:
# (número de átomos de carbono, índice do composto)
C5_C19_idx = [
    (5,  iC5),
    (6,  iC6),
    (7,  iC7),
    (8,  iC8),
    (9,  iC9),
    (10, iC10),
    (11, iC11),
    (12, iC12),
    (13, iC13),
    (14, iC14),
    (15, iC15),
    (16, iC16),
    (17, iC17),
    (18, iC18),
    (19, iC19)
]

# Fração pesada representada por lump C24
heavy_n = 24
heavy_idx = iC24

# --------------------------------------------------------------
# 3. PARÂMETROS DO REATOR
# --------------------------------------------------------------
# Wcat_total = massa total de catalisador (kg)
# Nsteps     = número de incrementos da discretização
# dW         = incremento diferencial de massa catalítica
Wcat_total = 5.0
Nsteps = 200
dW = Wcat_total / float(Nsteps)

# --------------------------------------------------------------
# 4. PARÂMETROS CINÉTICOS E CONSTANTES
# --------------------------------------------------------------
E = 2.718281828459045  # base da exponencial
R = 8.314              # constante universal dos gases [J/mol/K]

# 4.1. Parâmetros da reação principal de Fischer–Tropsch (FT)
# kref   = constante cinética de referência
# Ea_FT  = energia de ativação [J/mol]
# a0, b, f0 = coeficientes do denominador LHHW
# kpH2O  = termo do efeito positivo da água no numerador
kref = 7.05
Ea_FT = 92000.0
a0 = 12.0
b = 1.10
f0 = 1.25
kpH2O = 0.1

# 4.2. Parâmetros do modelo de alfa (ASF)
# kalpha_ref = parâmetro de referência da correlação de alfa
# Ea_alpha   = energia de ativação da correlação de alfa [J/mol]
# zexp e yexp = expoentes positivos da forma adotada de Pandey et al.
kalpha_ref = 0.118
Ea_alpha = 4770.0
zexp = 0.170
yexp = 0.095

# 4.3. Parâmetro da reação de deslocamento gás–água (WGS)
kWGS = 119.1

# --------------------------------------------------------------
# 5. INÍCIO DO PROCESSAMENTO DA UNIDADE
# --------------------------------------------------------------
# A corrente de saída é inicialmente limpa e em seguida recebe uma
# cópia da corrente de entrada, de modo a preservar pressão,
# temperatura e demais propriedades antes da atualização composicional.
prod.Clear()
prod.Assign(feed)

# Leitura das condições globais da corrente de entrada
T = float(feed.GetTemperature())          # temperatura [K]
P_MPa = float(feed.GetPressure()) / 1.0e6 # pressão total convertida de Pa para MPa

# Composição molar global da alimentação
z = list(feed.GetOverallComposition())

# Vazão molar total da alimentação [mol/s]
Ftot = float(feed.GetMolarFlow())

# Se não houver fluxo, o DWSIM apenas recalcula a corrente de saída
if Ftot <= 0.0:
    prod.Calculate()

else:
    # ----------------------------------------------------------
    # 6. INICIALIZAÇÃO DAS VAZÕES MOLARES INDIVIDUAIS
    # ----------------------------------------------------------
    # Fi[i] = vazão molar da espécie i [mol/s]
    Fi = []
    for zi in z:
        Fi.append(Ftot * zi)

    # ----------------------------------------------------------
    # 7. CONVERSÃO DE UNIDADES DAS TAXAS
    # ----------------------------------------------------------
    # As constantes cinéticas foram tratadas na base:
    # kmol/(kgcat·h)
    # A integração do reator é feita em:
    # mol/(kgcat·s)
    #
    # Logo:
    # 1 kmol = 1000 mol
    # 1 h    = 3600 s
    conv = 1000.0 / 3600.0

    # ----------------------------------------------------------
    # 8. LOOP DE INTEGRAÇÃO NUMÉRICA DO PFR
    # ----------------------------------------------------------
    # O reator é discretizado em 200 incrementos de massa de
    # catalisador e resolvido por Euler explícito.
    for step in range(Nsteps):

        # Vazão molar total local no passo atual
        Ftot_now = sum(Fi)

        # Critério de interrupção numérica
        if Ftot_now <= 1.0e-18:
            break

        # ------------------------------------------------------
        # 8.1. Frações molares locais
        # ------------------------------------------------------
        yCO  = Fi[iCO]  / Ftot_now
        yH2  = Fi[iH2]  / Ftot_now
        yH2O = Fi[iH2O] / Ftot_now
        yCO2 = Fi[iCO2] / Ftot_now

        # ------------------------------------------------------
        # 8.2. Pressões parciais locais [MPa]
        # ------------------------------------------------------
        pCO  = yCO  * P_MPa
        pH2  = yH2  * P_MPa
        pH2O = yH2O * P_MPa
        pCO2 = yCO2 * P_MPa

        # ------------------------------------------------------
        # 8.3. Proteções numéricas
        # ------------------------------------------------------
        # Evitam divisão por zero e problemas com potências.
        pCOeff = pCO
        if pCOeff < 1.0e-12:
            pCOeff = 1.0e-12

        pH2eff = pH2
        if pH2eff < 0.0:
            pH2eff = 0.0

        pH2Oeff = pH2O
        if pH2Oeff < 1.0e-12:
            pH2Oeff = 1.0e-12

        # ------------------------------------------------------
        # 8.4. Taxa da reação principal de Fischer–Tropsch
        # ------------------------------------------------------
        # kFT(T) = kref * exp[-Ea_FT/R * (1/T - 1/483)]
        #
        # rFT =
        # [kFT * pCO * sqrt(pH2) * (1 + kpH2O*pH2O)] /
        # [1 + a0*pCO + b*sqrt(pH2) + f0*pH2O]^2
        kFT = kref * (E ** (-(Ea_FT / R) * (1.0 / T - 1.0 / 483.0)))

        sqrt_pH2 = pH2eff ** 0.5

        num = kFT * pCOeff * sqrt_pH2 * (1.0 + kpH2O * pH2Oeff)
        den = 1.0 + a0 * pCOeff + b * sqrt_pH2 + f0 * pH2Oeff
        den = den * den

        rft = 0.0
        if den > 0.0:
            rft = num / den
        if rft < 0.0:
            rft = 0.0

        # Conversão para mol/(kgcat·s)
        rft = rft * conv

        # ------------------------------------------------------
        # 8.5. Cálculo do parâmetro alfa da distribuição ASF
        # ------------------------------------------------------
        # kalpha(T) = kalpha_ref * exp[-Ea_alpha/R * (1/T - 1/483)]
        #
        # alpha = 1 / [1 + kalpha(T) * 1/(pCO^z * pH2O^y)]
        kalphaT = kalpha_ref * (E ** (-(Ea_alpha / R) * (1.0 / T - 1.0 / 483.0)))
        denom_alpha = (pCOeff ** zexp) * (pH2Oeff ** yexp)
        a = 1.0 / (1.0 + kalphaT * (1.0 / denom_alpha))

        # Limitação numérica de alfa
        if a < 0.05:
            a = 0.05
        if a > 0.95:
            a = 0.95

        # ------------------------------------------------------
        # 8.6. Distribuição de produtos por ASF
        # ------------------------------------------------------
        # w_n = (1 - a) * a^(n - 1)
        w = {}
        for item in C5_C19_idx:
            n = item[0]
            w[n] = (1.0 - a) * (a ** (n - 1))

        # Carbono total ideal da distribuição
        carbon_total = 1.0 / (1.0 - a)

        # Carbono explicitamente alocado à faixa C5–C19
        carbon_C5_C19 = 0.0
        for item in C5_C19_idx:
            n = item[0]
            carbon_C5_C19 = carbon_C5_C19 + n * w[n]

        # Carbono remanescente agrupado no lump pesado C24
        carbon_tail = carbon_total - carbon_C5_C19
        if carbon_tail < 0.0:
            carbon_tail = 0.0

        # Taxas molares de formação dos hidrocarbonetos explícitos
        Rn = {}
        for item in C5_C19_idx:
            n = item[0]
            Rn[n] = (w[n] / carbon_total) * rft

        # Taxa molar do lump pesado
        R24 = (carbon_tail / carbon_total) * rft / float(heavy_n)

        # ------------------------------------------------------
        # 8.7. Consumo de H2 associado à FT
        # ------------------------------------------------------
        # Para CnH2n+2:
        # nCO + (2n+1)H2 -> CnH2n+2 + nH2O
        rH2_FT = 0.0
        for item in C5_C19_idx:
            n = item[0]
            rH2_FT = rH2_FT + (2 * n + 1) * Rn[n]
        rH2_FT = rH2_FT + (2 * heavy_n + 1) * R24

        # ------------------------------------------------------
        # 8.8. Limitador da FT
        # ------------------------------------------------------
        # Evita consumo maior do que o disponível de CO e H2.
        scale_FT = 1.0

        need_CO = rft * dW
        if need_CO > 1.0e-18:
            scale_CO = Fi[iCO] / need_CO
            if scale_CO < scale_FT:
                scale_FT = scale_CO

        need_H2 = rH2_FT * dW
        if need_H2 > 1.0e-18:
            scale_H2 = Fi[iH2] / need_H2
            if scale_H2 < scale_FT:
                scale_FT = scale_H2

        if scale_FT < 0.0:
            scale_FT = 0.0
        if scale_FT > 1.0:
            scale_FT = 1.0

        if scale_FT < 1.0:
            rft = rft * scale_FT
            rH2_FT = rH2_FT * scale_FT

            for item in C5_C19_idx:
                n = item[0]
                Rn[n] = Rn[n] * scale_FT

            R24 = R24 * scale_FT

        # ------------------------------------------------------
        # 8.9. Taxa da reação de deslocamento gás–água (WGS)
        # ------------------------------------------------------
        # rWGS = kco2 * [pCO*pH2O - (pCO2*pH2)/Keq]
        kco2 = kWGS * (E ** (-47400.0 / (R * T)))
        Keq = E ** (4557.8 / T - 4.33)
        rwgs = kco2 * (pCO * pH2O - (pCO2 * pH2) / Keq)

        # Conversão para mol/(kgcat·s)
        rwgs = rwgs * conv

        # ------------------------------------------------------
        # 8.10. Limitador da WGS
        # ------------------------------------------------------
        if rwgs > 0.0:
            # Sentido direto:
            # CO + H2O -> CO2 + H2
            scale_WGS = 1.0

            need_CO_wgs = rwgs * dW
            if need_CO_wgs > 1.0e-18:
                scale1 = Fi[iCO] / need_CO_wgs
                if scale1 < scale_WGS:
                    scale_WGS = scale1

            need_H2O_wgs = rwgs * dW
            if need_H2O_wgs > 1.0e-18:
                scale2 = Fi[iH2O] / need_H2O_wgs
                if scale2 < scale_WGS:
                    scale_WGS = scale2

            if scale_WGS < 0.0:
                scale_WGS = 0.0
            if scale_WGS > 1.0:
                scale_WGS = 1.0

            rwgs = rwgs * scale_WGS

        else:
            # Sentido reverso:
            # CO2 + H2 -> CO + H2O
            scale_WGS = 1.0

            need_CO2_rev = (-rwgs) * dW
            if need_CO2_rev > 1.0e-18:
                scale1 = Fi[iCO2] / need_CO2_rev
                if scale1 < scale_WGS:
                    scale_WGS = scale1

            need_H2_rev = (-rwgs) * dW
            if need_H2_rev > 1.0e-18:
                scale2 = Fi[iH2] / need_H2_rev
                if scale2 < scale_WGS:
                    scale_WGS = scale2

            if scale_WGS < 0.0:
                scale_WGS = 0.0
            if scale_WGS > 1.0:
                scale_WGS = 1.0

            rwgs = rwgs * scale_WGS

        # ------------------------------------------------------
        # 8.11. Atualização das vazões molares
        # ------------------------------------------------------
        # Contribuições da FT
        dCO_FT = -rft * dW
        dH2_FT = -rH2_FT * dW
        dH2O_FT = rft * dW

        # Contribuições da WGS
        dCO_WGS = -rwgs * dW
        dH2O_WGS = -rwgs * dW
        dCO2_WGS = rwgs * dW
        dH2_WGS = rwgs * dW

        # Atualização das espécies gasosas
        Fi[iCO]  = Fi[iCO]  + dCO_FT  + dCO_WGS
        Fi[iH2]  = Fi[iH2]  + dH2_FT  + dH2_WGS
        Fi[iH2O] = Fi[iH2O] + dH2O_FT + dH2O_WGS
        Fi[iCO2] = Fi[iCO2] + dCO2_WGS

        # Impede valores negativos
        if Fi[iCO] < 0.0:
            Fi[iCO] = 0.0
        if Fi[iH2] < 0.0:
            Fi[iH2] = 0.0
        if Fi[iH2O] < 0.0:
            Fi[iH2O] = 0.0
        if Fi[iCO2] < 0.0:
            Fi[iCO2] = 0.0

        # Atualização dos hidrocarbonetos explícitos C5–C19
        for item in C5_C19_idx:
            n = item[0]
            idx = item[1]
            Fi[idx] = Fi[idx] + Rn[n] * dW
            if Fi[idx] < 0.0:
                Fi[idx] = 0.0

        # Atualização do lump pesado C24
        Fi[heavy_idx] = Fi[heavy_idx] + R24 * dW
        if Fi[heavy_idx] < 0.0:
            Fi[heavy_idx] = 0.0

    # ----------------------------------------------------------
    # 9. ESCRITA DA CORRENTE DE SAÍDA NO DWSIM
    # ----------------------------------------------------------
    # Ao final da integração, o script calcula a nova vazão molar
    # total e escreve espécie por espécie na corrente de saída.
    Ftot_out = sum(Fi)

    if Ftot_out <= 0.0:
        prod.Calculate()
    else:
        prod.SetMolarFlow(Ftot_out)

        for i in range(len(Fi)):
            prod.SetOverallCompoundMolarFlow(i, Fi[i])

        # Recalcula as propriedades termodinâmicas do efluente
        prod.Calculate()
