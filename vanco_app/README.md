# Vancomycin AUC-Guided Dosing Calculator (Streamlit)

## Run
```bash
pip install -r requirements.txt
streamlit run app.py        # opens http://localhost:8501
pytest -q                   # 10 unit tests for the PK engine
```

## Files
| File | Purpose |
|---|---|
| `vanco_pk.py` | Pure-Python dosing engine (no UI). All tunable constants are at the top. |
| `app.py` | Streamlit interface |
| `test_vanco_pk.py` | Unit tests |

## What it covers
- **Adults (stable renal function):** Cockcroft-Gault CrCl → Matzke CL (0.689·CrCl + 3.66 mL/min), Vd 0.7 L/kg TBW.
  Rounded regimens for q8, q12, q24 and q48 h targeting AUC24 400–600, a recommended option, and a custom-regimen
  concentration–time simulation that includes the loading dose.
- **Obesity (BMI ≥ 30):** AdjBW for CrCl (auto rule), Vd 0.5 L/kg TBW, loading dose limited to 20–25 mg/kg with a 3 g cap.
- **Unstable renal function / AKI:** loading dose, then level-based dosing (with a warning).
- **Intermittent HD:** low-/high-permeability dialyzers, post- or intradialytic dosing, pre-HD target 15–20 mg/L.
- **CRRT:** 20–25 mg/kg load, then 7.5–10 mg/kg q12h.
- **Pediatrics (3 mo to <18 y):** bedside Schwartz eGFR. Ages 3 mo to <12 y: 60–80 mg/kg/day q6h. Ages ≥12 y: 60–70 mg/kg/day q6–8h.
  Obese children: 20 mg/kg load, then 60 mg/kg/day. Maximum 3,600 mg/day.

## Ideas for later
- Levels-based AUC (Sawchuk-Zaske from a peak and a trough)
- Bayesian estimation (MAP) using a published population model
- A PDF or printable dosing note

> Clinical decision support only. Verify against your institutional protocol before clinical use.
