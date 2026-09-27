"""
vanco_pk.py - Vancomycin dosing engine (no UI code).

Every calculation lives here so it can be unit-tested and reused
(e.g. in a notebook, a CLI, or the Streamlit app in app.py).

Scope
-----
* Adults: AUC24-guided empiric dosing from population PK
  (Cockcroft-Gault CrCl -> Matzke clearance, weight-based Vd),
  per the 2020 ASHP/IDSA/PIDS/SIDP consensus guideline (AUC24/MIC 400-600).
* Obesity: body-weight selection for CrCl, reduced Vd, loading-dose caps.
* Renal replacement: intermittent hemodialysis (IHD) and CRRT
  weight-based regimens from the 2020 guideline.
* Pediatrics (>= 3 months): guideline mg/kg/day regimens + bedside Schwartz eGFR.

NOT a substitute for clinical judgement or an institutional protocol.
All defaults are editable constants below.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# ----------------------------------------------------------------------------
# Tunable constants (edit to match your institutional protocol)
# ----------------------------------------------------------------------------
AUC_TARGET_DEFAULT = 500.0             # mg*h/L
AUC_RANGE_ADULT = (400.0, 600.0)
AUC_RANGE_PEDS = (400.0, 600.0)        # guideline: 400 (up to 600-800 may be ok)

VD_L_PER_KG = 0.7                      # non-obese adults, total body weight
VD_L_PER_KG_OBESE = 0.5                # BMI >= 30, total body weight
CRCL_CAP_FOR_PK = 150.0                # mL/min, stops runaway clearance estimates

LOADING_CAP_MG = 3000.0
MAX_SINGLE_DOSE_MG = 3000.0
NEPHROTOX_DAILY_FLAG_MG = 4000.0       # daily doses > 4 g linked to more AKI
ROUND_TO_MG = 250.0
INTERVALS_H = (8, 12, 24, 48)

# Intermittent hemodialysis (2020 guideline). mg/kg actual body weight.
IHD_TABLE = {
    ("low", "post"):   {"load": (25, 25), "maint": (7.5, 7.5)},
    ("high", "post"):  {"load": (25, 25), "maint": (10, 10)},
    ("high", "intra"): {"load": (30, 35), "maint": (10, 15)},
}
IHD_PRE_LEVEL_TARGET = (15, 20)        # mg/L pre-dialysis

# CRRT (CVVH/CVVHD/CVVHDF at usual effluent rates 20-25 mL/kg/h)
CRRT_LOAD_MG_KG = (20, 25)
CRRT_MAINT_MG_KG_Q12 = (7.5, 10)

# Pediatrics (>= 3 months, normal renal function)
PEDS_MAX_DAILY_MG = 3600.0
PEDS_OBESE_LOAD_MG_KG = 20.0
PEDS_OBESE_MAINT_MG_KG_DAY = 60.0


# ----------------------------------------------------------------------------
# Body size & renal function
# ----------------------------------------------------------------------------
def ibw_kg(height_cm: float, sex: str) -> float:
    """Devine ideal body weight. For height < 152.4 cm returns the base weight."""
    inches_over_5ft = max(0.0, height_cm / 2.54 - 60.0)
    base = 50.0 if sex.upper().startswith("M") else 45.5
    return base + 2.3 * inches_over_5ft


def adjbw_kg(tbw: float, ibw: float, factor: float = 0.4) -> float:
    return ibw + factor * (tbw - ibw)


def bmi(tbw: float, height_cm: float) -> float:
    return tbw / (height_cm / 100.0) ** 2


def crcl_weight(tbw: float, ibw: float, method: str = "auto") -> tuple[float, str]:
    """Weight used in Cockcroft-Gault.

    auto: TBW if TBW < IBW; AdjBW if TBW > 120% IBW; otherwise IBW.
    """
    method = method.lower()
    if method == "tbw":
        return tbw, "TBW"
    if method == "ibw":
        return ibw, "IBW"
    if method == "adjbw":
        return adjbw_kg(tbw, ibw), "AdjBW"
    if tbw < ibw:
        return tbw, "TBW (TBW < IBW)"
    if tbw > 1.2 * ibw:
        return adjbw_kg(tbw, ibw), "AdjBW (TBW > 120% IBW)"
    return ibw, "IBW"


def cockcroft_gault(age: float, weight_kg: float, scr: float, sex: str) -> float:
    """CrCl in mL/min."""
    crcl = (140.0 - age) * weight_kg / (72.0 * scr)
    return crcl * (0.85 if sex.upper().startswith("F") else 1.0)


def bedside_schwartz(height_cm: float, scr: float) -> float:
    """Pediatric eGFR in mL/min/1.73 m^2 (SCr in mg/dL)."""
    return 0.413 * height_cm / scr


# ----------------------------------------------------------------------------
# Pharmacokinetics (one-compartment, intermittent infusion)
# ----------------------------------------------------------------------------
@dataclass
class PKParams:
    cl_l_h: float
    vd_l: float

    @property
    def ke(self) -> float:
        return self.cl_l_h / self.vd_l

    @property
    def t_half(self) -> float:
        return math.log(2) / self.ke


def adult_population_pk(crcl: float, tbw: float, obese: bool) -> PKParams:
    """Matzke: CL(mL/min) = 0.689*CrCl + 3.66 ; Vd weight-based on TBW."""
    crcl_eff = min(max(crcl, 0.0), CRCL_CAP_FOR_PK)
    cl = (0.689 * crcl_eff + 3.66) * 0.06   # -> L/h
    vd = (VD_L_PER_KG_OBESE if obese else VD_L_PER_KG) * tbw
    return PKParams(cl_l_h=cl, vd_l=vd)


def infusion_hours(dose_mg: float) -> float:
    """~1 g/h, rounded up to the next 0.5 h, minimum 1 h."""
    return max(1.0, math.ceil(dose_mg / 500.0) * 0.5)


def steady_state(dose_mg: float, tau_h: float, pk: PKParams,
                 tinf_h: float | None = None) -> dict:
    """Steady-state Cmax (end of infusion), Cmin (trough) and AUC24."""
    tinf = tinf_h if tinf_h is not None else infusion_hours(dose_mg)
    k = pk.ke
    rate = dose_mg / tinf
    cmax = (rate / pk.cl_l_h) * (1 - math.exp(-k * tinf)) / (1 - math.exp(-k * tau_h))
    cmin = cmax * math.exp(-k * (tau_h - tinf))
    auc24 = dose_mg * (24.0 / tau_h) / pk.cl_l_h
    return {"cmax": cmax, "cmin": cmin, "auc24": auc24, "tinf": tinf}


def round_dose(mg: float, step: float = ROUND_TO_MG) -> float:
    return max(step, round(mg / step) * step)


@dataclass
class Regimen:
    dose_mg: float
    tau_h: float
    auc24: float
    cmax: float
    cmin: float
    tinf_h: float
    in_range: bool
    notes: list[str] = field(default_factory=list)

    @property
    def daily_mg(self) -> float:
        return self.dose_mg * 24.0 / self.tau_h

    @property
    def label(self) -> str:
        return f"{self.dose_mg:,.0f} mg q{self.tau_h:g}h"


def evaluate_regimen(dose_mg: float, tau_h: float, pk: PKParams,
                     auc_range=AUC_RANGE_ADULT,
                     max_single=MAX_SINGLE_DOSE_MG) -> Regimen:
    ss = steady_state(dose_mg, tau_h, pk)
    notes = []
    if dose_mg > max_single:
        notes.append(f"exceeds max single dose ({max_single:,.0f} mg)")
    if dose_mg * 24 / tau_h > NEPHROTOX_DAILY_FLAG_MG:
        notes.append("daily dose > 4 g: higher AKI risk")
    lo, hi = auc_range
    return Regimen(dose_mg, tau_h, ss["auc24"], ss["cmax"], ss["cmin"],
                   ss["tinf"], lo <= ss["auc24"] <= hi and dose_mg <= max_single, notes)


def candidate_regimens(pk: PKParams, target_auc=AUC_TARGET_DEFAULT,
                       auc_range=AUC_RANGE_ADULT, intervals=INTERVALS_H,
                       step=ROUND_TO_MG, max_single=MAX_SINGLE_DOSE_MG) -> list[Regimen]:
    """One rounded regimen per dosing interval, aimed at the AUC target."""
    out = []
    for tau in intervals:
        ideal = target_auc * pk.cl_l_h * tau / 24.0
        out.append(evaluate_regimen(round_dose(ideal, step), tau, pk, auc_range, max_single))
    return out


def recommend(regs: list[Regimen], pk: PKParams, target_auc=AUC_TARGET_DEFAULT) -> Regimen | None:
    """Pick an in-range regimen whose interval best matches ~1.5x half-life,
    tie-breaking on closeness to the AUC target."""
    ok = [r for r in regs if r.in_range]
    if not ok:
        return None
    ideal_tau = min(48.0, max(8.0, 1.5 * pk.t_half))
    return min(ok, key=lambda r: (abs(math.log(r.tau_h / ideal_tau)),
                                  abs(r.auc24 - target_auc)))


def loading_dose(weight_kg: float, mg_per_kg: float = 25.0,
                 cap=LOADING_CAP_MG, step=ROUND_TO_MG) -> float:
    return min(cap, round_dose(weight_kg * mg_per_kg, step))


def concentration_curve(doses: list[tuple[float, float, float]], pk: PKParams,
                        t_end: float, n: int = 600) -> tuple[list[float], list[float]]:
    """Superposition of infusions. doses = [(start_h, amount_mg, tinf_h), ...]."""
    k, cl = pk.ke, pk.cl_l_h
    ts, cs = [], []
    for i in range(n + 1):
        t = t_end * i / n
        c = 0.0
        for t0, amt, tinf in doses:
            if t <= t0:
                continue
            r = amt / tinf
            if t <= t0 + tinf:
                c += r / cl * (1 - math.exp(-k * (t - t0)))
            else:
                c += r / cl * (1 - math.exp(-k * tinf)) * math.exp(-k * (t - t0 - tinf))
        ts.append(t)
        cs.append(c)
    return ts, cs


def build_schedule(maint_dose: float, tau: float, load_dose: float | None,
                   t_end: float) -> list[tuple[float, float, float]]:
    """Loading dose at t=0 (if any), first maintenance dose one interval later."""
    doses = []
    t = 0.0
    if load_dose:
        doses.append((0.0, load_dose, infusion_hours(load_dose)))
        t = tau
    while t < t_end:
        doses.append((t, maint_dose, infusion_hours(maint_dose)))
        t += tau
    return doses


# ----------------------------------------------------------------------------
# Special populations
# ----------------------------------------------------------------------------
def _rng_mg(weight, mgkg_range, step=ROUND_TO_MG, cap=None):
    lo, hi = (round_dose(weight * x, step) for x in mgkg_range)
    if cap:
        lo, hi = min(lo, cap), min(hi, cap)
    return lo, hi


def ihd_regimen(tbw: float, flux: str = "high", timing: str = "post") -> dict:
    key = (flux, timing)
    if key not in IHD_TABLE:
        raise ValueError("Intradialytic dosing is only defined for high-flux dialyzers.")
    row = IHD_TABLE[key]
    return {
        "load_mg_kg": row["load"], "maint_mg_kg": row["maint"],
        "load_mg": _rng_mg(tbw, row["load"], cap=LOADING_CAP_MG),
        "maint_mg": _rng_mg(tbw, row["maint"]),
        "pre_hd_target": IHD_PRE_LEVEL_TARGET,
    }


def crrt_regimen(tbw: float) -> dict:
    return {
        "load_mg_kg": CRRT_LOAD_MG_KG, "maint_mg_kg": CRRT_MAINT_MG_KG_Q12,
        "load_mg": _rng_mg(tbw, CRRT_LOAD_MG_KG, cap=LOADING_CAP_MG),
        "maint_mg": _rng_mg(tbw, CRRT_MAINT_MG_KG_Q12),
    }


def pediatric_regimen(age_years: float, tbw: float, obese: bool = False,
                      critically_ill: bool = False,
                      max_daily=PEDS_MAX_DAILY_MG) -> dict:
    """Initial empiric regimen for children >= 3 months with normal renal function."""
    if age_years < 0.25:
        raise ValueError("Infants < 3 months / neonates need PMA-based dosing - not covered.")
    if age_years >= 18:
        raise ValueError("Use adult dosing for age >= 18 years.")

    if obese:
        mgkg_day, intervals = (PEDS_OBESE_MAINT_MG_KG_DAY, PEDS_OBESE_MAINT_MG_KG_DAY), (6, 8)
    elif age_years < 12:
        mgkg_day, intervals = (60, 80), (6,)
    else:
        mgkg_day, intervals = (60, 70), (6, 8)

    options = []
    for tau in intervals:
        per_day = 24 / tau
        for mgkg in sorted(set(mgkg_day)):
            daily = min(tbw * mgkg, max_daily)
            dose = round(daily / per_day / 5) * 5          # round to 5 mg for kids
            if any(o["tau_h"] == tau and o["dose_mg"] == dose for o in options):
                continue                                    # identical after capping
            options.append({"mg_kg_day": mgkg, "tau_h": tau, "dose_mg": dose,
                            "daily_mg": dose * per_day,
                            "capped": tbw * mgkg > max_daily})

    load = None
    if obese:
        load = min(LOADING_CAP_MG, round(tbw * PEDS_OBESE_LOAD_MG_KG / 5) * 5)
    elif critically_ill:
        load = (min(LOADING_CAP_MG, round(tbw * 20 / 5) * 5),
                min(LOADING_CAP_MG, round(tbw * 35 / 5) * 5))
    return {"mg_kg_day": mgkg_day, "options": options, "loading": load,
            "auc_range": AUC_RANGE_PEDS}
