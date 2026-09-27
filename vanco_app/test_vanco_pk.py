import math

import pytest

import vanco_pk as vk


def test_ibw_and_adjbw():
    assert vk.ibw_kg(175, "M") == pytest.approx(50 + 2.3 * (175 / 2.54 - 60), abs=1e-6)  # ~70.5
    assert vk.ibw_kg(150, "F") == 45.5          # below 5 ft -> base weight
    assert vk.adjbw_kg(120, 70) == pytest.approx(90)


def test_cockcroft_gault():
    # 60 y male, 70 kg, SCr 1.0 -> (80*70)/72 = 77.8
    assert vk.cockcroft_gault(60, 70, 1.0, "M") == pytest.approx(77.78, abs=0.01)
    assert vk.cockcroft_gault(60, 70, 1.0, "F") == pytest.approx(77.78 * 0.85, abs=0.01)


def test_crcl_weight_rule():
    assert vk.crcl_weight(60, 70)[1].startswith("TBW")
    assert vk.crcl_weight(75, 70)[1] == "IBW"
    w, lab = vk.crcl_weight(120, 70)
    assert lab.startswith("AdjBW") and w == pytest.approx(90)


def test_population_pk_and_auc():
    pk = vk.adult_population_pk(crcl=100, tbw=80, obese=False)
    assert pk.cl_l_h == pytest.approx((0.689 * 100 + 3.66) * 0.06)   # 4.3536 L/h
    assert pk.vd_l == pytest.approx(56)
    # AUC24 = daily dose / CL
    ss = vk.steady_state(1250, 12, pk)
    assert ss["auc24"] == pytest.approx(2500 / pk.cl_l_h)
    assert ss["tinf"] == 1.5


def test_steady_state_matches_superposition():
    pk = vk.adult_population_pk(80, 70, False)
    dose, tau = 1000, 12
    ss = vk.steady_state(dose, tau, pk)
    sched = vk.build_schedule(dose, tau, None, 24 * 12)
    ts, cs = vk.concentration_curve(sched, pk, 24 * 12, n=24 * 12 * 4)
    # the trough just before the last dose at t=276 h should approach Css,min
    idx = ts.index(276.0)
    assert cs[idx] == pytest.approx(ss["cmin"], rel=0.01)
    # peak at end of infusion (t = 264 + 1)
    assert cs[ts.index(265.0)] == pytest.approx(ss["cmax"], rel=0.01)


def test_recommendation_in_range():
    pk = vk.adult_population_pk(100, 80, False)
    regs = vk.candidate_regimens(pk)
    best = vk.recommend(regs, pk)
    assert best is not None and 400 <= best.auc24 <= 600
    assert best.tau_h == 12                     # t1/2 ~ 8.9 h -> q12h


def test_poor_renal_function_prefers_long_interval():
    pk = vk.adult_population_pk(20, 70, False)
    best = vk.recommend(vk.candidate_regimens(pk), pk)
    assert best.tau_h >= 24


def test_loading_dose_cap():
    assert vk.loading_dose(80, 25) == 2000
    assert vk.loading_dose(150, 25) == 3000


def test_ihd_and_crrt():
    r = vk.ihd_regimen(80, "high", "post")
    assert r["load_mg"] == (2000, 2000) and r["maint_mg"] == (750, 750)
    with pytest.raises(ValueError):
        vk.ihd_regimen(80, "low", "intra")
    c = vk.crrt_regimen(80)
    assert c["load_mg"] == (1500, 2000) and c["maint_mg"] == (500, 750)  # 600->500, 800->750


def test_pediatrics():
    r = vk.pediatric_regimen(6, 20)
    doses = {(o["mg_kg_day"], o["tau_h"]): o["dose_mg"] for o in r["options"]}
    assert doses[(60, 6)] == 300 and doses[(80, 6)] == 400
    big = vk.pediatric_regimen(15, 70)
    assert all(o["daily_mg"] <= vk.PEDS_MAX_DAILY_MG for o in big["options"])
    ob = vk.pediatric_regimen(10, 60, obese=True)
    assert ob["loading"] == 1200
    with pytest.raises(ValueError):
        vk.pediatric_regimen(0.1, 4)
    assert vk.bedside_schwartz(120, 0.5) == pytest.approx(99.12)
