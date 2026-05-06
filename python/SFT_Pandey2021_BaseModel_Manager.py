import clr
clr.AddReference("System")

feed = ims1
prod = oms1

# =========================================================
# ÍNDICES DOS COMPONENTES - CORRIGIDOS PARA A SUA ORDEM REAL
# =========================================================
iCO  = 0
iH2  = 1
iH2O = 2
iCO2 = 3

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
iC20 = 27
iC24 = 26

C5_C19_idx = [
    (5, iC5),
    (6, iC6),
    (7, iC7),
    (8, iC8),
    (9, iC9),
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

heavy_n = 24
heavy_idx = iC24

# =========================================================
# PARÂMETROS DO REATOR
# =========================================================
Wcat_total = 5
Nsteps = 200
dW = Wcat_total / float(Nsteps)

# =========================================================
# PARÂMETROS CINÉTICOS
# =========================================================
E = 2.718281828459045
R = 8.314

# FT
kref = 7.05
Ea_FT = 92000.0
a0 = 12.0
b = 1.10
f0 = 1.25
kpH2O = 0.1

# alpha
kalpha_ref = 0.118
Ea_alpha = 4770.0
zexp = 0.170
yexp = 0.095

# WGS
kWGS = 119.1

# =========================================================
# INÍCIO
# =========================================================
prod.Clear()
prod.Assign(feed)

T = float(feed.GetTemperature())
P_MPa = float(feed.GetPressure()) / 1.0e6

z = list(feed.GetOverallComposition())
Ftot = float(feed.GetMolarFlow())

if Ftot <= 0.0:
    prod.Calculate()
else:
    Fi = []
    for zi in z:
        Fi.append(Ftot * zi)

    # paper em kmol/(kgcat*h) -> mol/(kgcat*s)
    conv = 1000.0 / 3600.0

    for step in range(Nsteps):
        Ftot_now = sum(Fi)
        if Ftot_now <= 1.0e-18:
            break

        # Frações molares atuais
        yCO = Fi[iCO] / Ftot_now
        yH2 = Fi[iH2] / Ftot_now
        yH2O = Fi[iH2O] / Ftot_now
        yCO2 = Fi[iCO2] / Ftot_now

        # Pressões parciais (MPa)
        pCO = yCO * P_MPa
        pH2 = yH2 * P_MPa
        pH2O = yH2O * P_MPa
        pCO2 = yCO2 * P_MPa

        # Proteções numéricas
        pCOeff = pCO
        if pCOeff < 1.0e-12:
            pCOeff = 1.0e-12

        pH2eff = pH2
        if pH2eff < 0.0:
            pH2eff = 0.0

        pH2Oeff = pH2O
        if pH2Oeff < 1.0e-12:
            pH2Oeff = 1.0e-12

        # =====================================================
        # FT
        # =====================================================
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

        rft = rft * conv

        # =====================================================
        # alpha / ASF
        # =====================================================
        kalphaT = kalpha_ref * (E ** (-(Ea_alpha / R) * (1.0 / T - 1.0 / 483.0)))
        denom_alpha = (pCOeff ** zexp) * (pH2Oeff ** yexp)
        a = 1.0 / (1.0 + kalphaT * (1.0 / denom_alpha))

        if a < 0.05:
            a = 0.05
        if a > 0.95:
            a = 0.95

        w = {}
        for item in C5_C19_idx:
            n = item[0]
            w[n] = (1.0 - a) * (a ** (n - 1))

        carbon_total = 1.0 / (1.0 - a)

        carbon_C5_C19 = 0.0
        for item in C5_C19_idx:
            n = item[0]
            carbon_C5_C19 = carbon_C5_C19 + n * w[n]

        carbon_tail = carbon_total - carbon_C5_C19
        if carbon_tail < 0.0:
            carbon_tail = 0.0

        Rn = {}
        for item in C5_C19_idx:
            n = item[0]
            Rn[n] = (w[n] / carbon_total) * rft

        R24 = (carbon_tail / carbon_total) * rft / float(heavy_n)

        # H2 consumido na FT
        rH2_FT = 0.0
        for item in C5_C19_idx:
            n = item[0]
            rH2_FT = rH2_FT + (2 * n + 1) * Rn[n]
        rH2_FT = rH2_FT + (2 * heavy_n + 1) * R24

        # -----------------------------------------------------
        # LIMITADOR FT: não consumir mais CO ou H2 do que existe
        # -----------------------------------------------------
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

        # =====================================================
        # WGS
        # =====================================================
        kco2 = kWGS * (E ** (-47400.0 / (R * T)))
        Keq = E ** (4557.8 / T - 4.33)
        rwgs = kco2 * (pCO * pH2O - (pCO2 * pH2) / Keq)
        rwgs = rwgs * conv

        # -----------------------------------------------------
        # LIMITADOR WGS: respeitar disponibilidade de reagentes
        # -----------------------------------------------------
        if rwgs > 0.0:
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
            # sentido reverso: CO2 + H2 -> CO + H2O
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

        # =====================================================
        # ATUALIZAÇÕES
        # =====================================================
        dCO_FT = -rft * dW
        dH2_FT = -rH2_FT * dW
        dH2O_FT = rft * dW

        dCO_WGS = -rwgs * dW
        dH2O_WGS = -rwgs * dW
        dCO2_WGS = rwgs * dW
        dH2_WGS = rwgs * dW

        Fi[iCO] = Fi[iCO] + dCO_FT + dCO_WGS
        Fi[iH2] = Fi[iH2] + dH2_FT + dH2_WGS
        Fi[iH2O] = Fi[iH2O] + dH2O_FT + dH2O_WGS
        Fi[iCO2] = Fi[iCO2] + dCO2_WGS

        if Fi[iCO] < 0.0:
            Fi[iCO] = 0.0
        if Fi[iH2] < 0.0:
            Fi[iH2] = 0.0
        if Fi[iH2O] < 0.0:
            Fi[iH2O] = 0.0
        if Fi[iCO2] < 0.0:
            Fi[iCO2] = 0.0

        for item in C5_C19_idx:
            n = item[0]
            idx = item[1]
            Fi[idx] = Fi[idx] + Rn[n] * dW
            if Fi[idx] < 0.0:
                Fi[idx] = 0.0

        Fi[heavy_idx] = Fi[heavy_idx] + R24 * dW
        if Fi[heavy_idx] < 0.0:
            Fi[heavy_idx] = 0.0

    # =========================================================
    # SAÍDA
    # =========================================================
    Ftot_out = sum(Fi)

    if Ftot_out <= 0.0:
        prod.Calculate()
    else:
        prod.SetMolarFlow(Ftot_out)
        for i in range(len(Fi)):
            prod.SetOverallCompoundMolarFlow(i, Fi[i])
        prod.Calculate()