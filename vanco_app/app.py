"""
Vancomycin AUC-guided dosing calculator - Streamlit UI.

Run:  streamlit run app.py
All maths lives in vanco_pk.py.
"""
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import vanco_pk as vk

st.set_page_config(page_title="Vancomycin Dosing", page_icon="💉", layout="wide")

st.title("Vancomycin AUC-Guided Dosing Calculator")
st.caption("Based on the 2020 ASHP/IDSA/PIDS/SIDP consensus guideline (target AUC24/MIC 400-600, "
           "assuming MIC 1 mg/L). This tool supports clinical decisions and does not replace "
           "clinical judgement or your institution's protocol.")

# ----------------------------------------------------------------------------
# Sidebar - patient inputs
# ----------------------------------------------------------------------------
with st.sidebar:
    st.header("Patient")
    population = st.radio("Population", ["Adult (≥18 y)", "Pediatric (3 mo – <18 y)"], key="pop")
    is_peds = population.startswith("Pediatric")

    renal_opts = (["Stable renal function", "Unstable / AKI"] if is_peds else
                  ["Stable renal function", "Unstable / AKI",
                   "Intermittent hemodialysis", "CRRT"])
    renal = st.selectbox("Renal status", renal_opts, key="renal")

    if is_peds:
        age = st.number_input("Age (years)", 0.25, 17.99, 6.0, 0.25,
                              help="Use decimals for months, e.g. 0.5 = 6 months")
    else:
        age = st.number_input("Age (years)", 18, 110, 55)
    sex = st.radio("Sex", ["Male", "Female"], horizontal=True)
    weight = st.number_input("Actual body weight (kg)", 2.0, 350.0, 20.0 if is_peds else 80.0, 0.5, key=f"wt_{is_peds}")
    height = st.number_input("Height (cm)", 45.0, 230.0, 115.0 if is_peds else 175.0, 0.5)
    scr = st.number_input("Serum creatinine (mg/dL)", 0.1, 20.0, 0.4 if is_peds else 1.0, 0.1)

    if not is_peds:
        st.header("Options")
        scr_floor = st.selectbox("Round SCr up to", ["No rounding", "0.7", "0.8", "1.0"],
                                 help="Some protocols round low SCr up in elderly or low-muscle-mass patients.")
        wt_method = st.selectbox("CrCl weight", ["auto", "tbw", "ibw", "adjbw"],
                                 format_func=lambda m: {"auto": "Auto (TBW / IBW / AdjBW rule)",
                                                        "tbw": "Actual body weight", "ibw": "Ideal body weight",
                                                        "adjbw": "Adjusted body weight"}[m])
        target_auc = st.slider("Target AUC24 (mg·h/L)", 400, 600, 500, 25)
        load_mgkg = st.slider("Loading dose (mg/kg actual weight)", 20, 35, 25)
        max_single = st.number_input("Max single dose (mg)", 1000, 4000, int(vk.MAX_SINGLE_DOSE_MG), 250)
    else:
        peds_obese = st.checkbox("Obese (BMI ≥ 95th percentile for age)")
        peds_icu = st.checkbox("Critically ill (consider loading dose)")

bmi = vk.bmi(weight, height)

# ----------------------------------------------------------------------------
# Pediatric pathway
# ----------------------------------------------------------------------------
if is_peds:
    egfr = vk.bedside_schwartz(height, scr)
    c1, c2, c3 = st.columns(3)
    c1.metric("BMI", f"{bmi:.1f} kg/m²")
    c2.metric("eGFR (bedside Schwartz)", f"{egfr:.0f} mL/min/1.73m²")
    c3.metric("Dosing weight", f"{weight:.1f} kg (actual)")

    if renal == "Unstable / AKI" or egfr < 60:
        st.error("Reduced or unstable renal function: the standard pediatric mg/kg/day "
                 "regimens below assume normal renal function. Give a single dose, "
                 "then use measured levels to choose the next dose (consult pharmacy / ID).")

    reg = vk.pediatric_regimen(age, weight, obese=peds_obese, critically_ill=peds_icu)
    st.subheader("Initial empiric regimen")
    if reg["loading"]:
        if isinstance(reg["loading"], tuple):
            st.info(f"**Loading dose (optional, critically ill):** 20–35 mg/kg → "
                    f"{reg['loading'][0]:,.0f}–{reg['loading'][1]:,.0f} mg (max 3,000 mg)")
        else:
            st.info(f"**Loading dose:** 20 mg/kg actual body weight → {reg['loading']:,.0f} mg")

    df = pd.DataFrame(reg["options"])
    df["Regimen"] = df.apply(lambda r: f"{r.dose_mg:,.0f} mg q{r.tau_h:g}h", axis=1)
    df = df[["Regimen", "mg_kg_day", "daily_mg", "capped"]].rename(columns={
        "mg_kg_day": "mg/kg/day", "daily_mg": "Daily dose (mg)", "capped": f"Capped at {vk.PEDS_MAX_DAILY_MG:,.0f} mg/day"})
    st.dataframe(df, hide_index=True, width="stretch")

    st.markdown(f"""
**Monitoring**
- Target AUC24 **{reg['auc_range'][0]:.0f}–{reg['auc_range'][1]:.0f} mg·h/L**. Children clear
  vancomycin quickly and vary a lot, so confirm the AUC with **two levels** (or Bayesian software)
  within the first 24–48 h.
- Monitor SCr and urine output at least every 2–3 days (daily if critically ill or on other nephrotoxins).
""")
    if peds_obese:
        st.warning("Obese children: the loading dose uses actual weight, but maintenance "
                   "doses should be conservative. Check levels early.")
    st.caption("Neonates and infants under 3 months need dosing based on postmenstrual age and are not covered here.")
    st.stop()

# ----------------------------------------------------------------------------
# Adult pathway
# ----------------------------------------------------------------------------
ibw = vk.ibw_kg(height, sex)
adj = vk.adjbw_kg(weight, ibw)
obese = bmi >= 30

scr_used = scr if scr_floor == "No rounding" else max(scr, float(scr_floor))
cg_wt, cg_label = vk.crcl_weight(weight, ibw, wt_method)
crcl = vk.cockcroft_gault(age, cg_wt, scr_used, sex)

m = st.columns(5)
m[0].metric("BMI", f"{bmi:.1f}", "Obese" if obese else None, delta_color="off", delta_arrow="off")
m[1].metric("IBW", f"{ibw:.1f} kg")
m[2].metric("AdjBW", f"{adj:.1f} kg")
m[3].metric("CrCl (C-G)", f"{crcl:.0f} mL/min", cg_label, delta_color="off", delta_arrow="off")
m[4].metric("SCr used", f"{scr_used:.2f} mg/dL")

load_wt_mgkg = min(load_mgkg, 25) if obese else load_mgkg
ld = vk.loading_dose(weight, load_wt_mgkg)

# ---- Renal replacement therapy --------------------------------------------
if renal == "Intermittent hemodialysis":
    st.subheader("Intermittent hemodialysis")
    c1, c2 = st.columns(2)
    flux = c1.radio("Dialyzer", ["high", "low"], format_func=lambda f: f"{f.title()}-permeability (flux)",
                    horizontal=True)
    timing = c2.radio("Dose timing", ["post", "intra"] if flux == "high" else ["post"],
                      format_func=lambda t: {"post": "After HD", "intra": "During last hour of HD"}[t],
                      horizontal=True)
    r = vk.ihd_regimen(weight, flux, timing)

    def rng(lo, hi, fmt="{:g}"):
        return fmt.format(lo) if lo == hi else f"{fmt.format(lo)}–{fmt.format(hi)}"

    a, b = st.columns(2)
    a.success(f"**Loading dose:** {rng(*r['load_mg_kg'])} mg/kg → "
              f"**{rng(*r['load_mg'], '{:,.0f}')} mg**")
    b.success(f"**Maintenance:** {rng(*r['maint_mg_kg'])} mg/kg per HD session → "
              f"**{rng(*r['maint_mg'], '{:,.0f}')} mg**")
    st.markdown(f"""
- Check a **pre-dialysis level** before the 3rd or 4th session and adjust to keep it at
  **{r['pre_hd_target'][0]}–{r['pre_hd_target'][1]} mg/L** (this corresponds to an AUC of about 400–600).
- Hold or reduce the dose if the pre-HD level is above 25 mg/L. Give a larger dose if HD is delayed or the session is longer than usual.
- Doses are based on actual body weight. The loading dose is capped at 3,000 mg.
""")
    st.stop()

if renal == "CRRT":
    st.subheader("Continuous renal replacement therapy (CVVH / CVVHD / CVVHDF)")
    r = vk.crrt_regimen(weight)
    a, b = st.columns(2)
    a.success(f"**Loading dose:** 20–25 mg/kg → **{r['load_mg'][0]:,.0f}–{r['load_mg'][1]:,.0f} mg**")
    b.success(f"**Maintenance:** 7.5–10 mg/kg q12h → **{r['maint_mg'][0]:,.0f}–{r['maint_mg'][1]:,.0f} mg q12h**")
    st.markdown("""
- Assumes usual effluent rates (20–25 mL/kg/h). Clearance changes with the CRRT dose and with interruptions to the circuit.
- Check levels within 24 h, then use AUC or trough monitoring. Re-check whenever the CRRT settings change.
""")
    st.stop()

# ---- Population-PK AUC dosing ----------------------------------------------
pk = vk.adult_population_pk(crcl, weight, obese)
k = st.columns(4)
k[0].metric("Est. clearance", f"{pk.cl_l_h:.2f} L/h")
k[1].metric("Est. Vd", f"{pk.vd_l:.0f} L", f"{vk.VD_L_PER_KG_OBESE if obese else vk.VD_L_PER_KG} L/kg TBW",
            delta_color="off", delta_arrow="off")
k[2].metric("ke", f"{pk.ke:.3f} h⁻¹")
k[3].metric("Half-life", f"{pk.t_half:.1f} h")

if renal == "Unstable / AKI":
    st.error("**Unstable renal function:** the Cockcroft-Gault CrCl is unreliable when SCr is changing. "
             "Give the loading dose, then use levels to choose the next dose (a random level at "
             "about 24 h, or two levels). The regimens below are only an estimate.")

regs = vk.candidate_regimens(pk, target_auc, vk.AUC_RANGE_ADULT, max_single=max_single)
best = vk.recommend(regs, pk, target_auc)

left, right = st.columns([1, 1.3])
with left:
    st.subheader("Recommendation")
    st.info(f"**Loading dose:** {ld:,.0f} mg once "
            f"({load_wt_mgkg} mg/kg × {weight:.0f} kg actual, max 3,000 mg)"
            + (" - obesity: 20–25 mg/kg" if obese else ""))
    if best:
        st.success(f"**Maintenance:** {best.label} (infuse over {best.tinf_h:g} h), "
                   f"starting {best.tau_h:g} h after the loading dose\n\n"
                   f"Predicted AUC24 **{best.auc24:.0f}** · peak {best.cmax:.1f} · trough {best.cmin:.1f} mg/L")
        for n in best.notes:
            st.warning(n)
    else:
        st.warning("No standard interval gives an AUC within range at or below the max single dose. "
                   "Use the custom regimen below or dose based on levels.")

    st.subheader("Options by interval")
    tbl = pd.DataFrame([{
        "Regimen": r.label, "Daily (mg)": r.daily_mg, "AUC24": round(r.auc24),
        "Peak": round(r.cmax, 1), "Trough": round(r.cmin, 1),
        "In range": "✅" if r.in_range else "—", "Notes": "; ".join(r.notes)} for r in regs])
    st.dataframe(tbl, hide_index=True, width="stretch")

with right:
    st.subheader("Custom regimen / simulation")
    c1, c2, c3 = st.columns(3)
    dflt = best or regs[1]
    dose = c1.number_input("Dose (mg)", 250, 4000, int(dflt.dose_mg), 250)
    tau = c2.selectbox("Interval (h)", [6, 8, 12, 18, 24, 36, 48],
                       index=[6, 8, 12, 18, 24, 36, 48].index(int(dflt.tau_h)))
    use_ld = c3.checkbox("Include loading dose", True)
    cust = vk.evaluate_regimen(dose, tau, pk, max_single=max_single)
    st.write(f"Predicted steady state: **AUC24 {cust.auc24:.0f}** · peak {cust.cmax:.1f} · "
             f"trough {cust.cmin:.1f} mg/L · {cust.daily_mg:,.0f} mg/day")

    t_end = max(72.0, 4 * tau)
    sched = vk.build_schedule(dose, tau, ld if use_ld else None, t_end)
    ts, cs = vk.concentration_curve(sched, pk, t_end)
    fig = go.Figure()
    fig.add_hrect(y0=10, y1=20, fillcolor="green", opacity=0.08, line_width=0,
                  annotation_text="trough 10–20 (reference)", annotation_position="top left")
    fig.add_trace(go.Scatter(x=ts, y=cs, mode="lines", name="Predicted conc.", line=dict(width=2.5)))
    fig.update_layout(xaxis_title="Time (h)", yaxis_title="Vancomycin (mg/L)", height=380,
                      margin=dict(l=10, r=10, t=30, b=10), showlegend=False)
    st.plotly_chart(fig, width="stretch")

# ---- Flags & monitoring ------------------------------------------------------
flags = []
if obese:
    flags.append("BMI ≥ 30: Vd estimated at 0.5 L/kg TBW and CrCl uses AdjBW (auto). Population estimates are "
                 "less reliable in obesity, so get levels early.")
if crcl < 30:
    flags.append("CrCl < 30 mL/min: long half-life, so consider level-based dosing after the loading dose.")
if age >= 65 and scr < 0.8:
    flags.append("Elderly patient with low SCr: CrCl may be overestimated (consider SCr rounding).")
if crcl > 130:
    flags.append("Augmented renal clearance (CrCl > 130): may need more frequent dosing, so check levels early.")
if flags:
    st.subheader("Flags")
    for f in flags:
        st.warning(f)

with st.expander("Monitoring & method notes"):
    st.markdown(f"""
**Monitoring (2020 guideline)**
- Aim for an AUC24/MIC of 400–600 (MIC 1). Estimate the AUC from two levels (a peak 1–2 h after the end of the
  infusion plus a trough) or with Bayesian software, ideally within 24–48 h.
- Trough-only monitoring is no longer recommended for serious MRSA infections.
- Check SCr at least every 2–3 days, and daily if critically ill, on nephrotoxins, or in unstable renal function.

**Method**
- CrCl: Cockcroft-Gault (×0.85 female). Auto weight rule: TBW if TBW < IBW, AdjBW (IBW + 0.4 × excess) if
  TBW > 120% IBW, otherwise IBW. CrCl is capped at {vk.CRCL_CAP_FOR_PK:.0f} mL/min for the PK estimate.
- CL (Matzke) = 0.689 × CrCl + 3.66 mL/min. Vd = {vk.VD_L_PER_KG} L/kg TBW
  ({vk.VD_L_PER_KG_OBESE} L/kg if BMI ≥ 30).
- AUC24 = daily dose / CL. Peak and trough come from one-compartment intermittent-infusion equations.
- Infusion time is about 1 g/h (minimum 1 h). Doses are rounded to {vk.ROUND_TO_MG:.0f} mg.
""")
