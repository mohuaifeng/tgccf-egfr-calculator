# TG-CCF Pediatric eGFR Calculator

Version: v1.5 layout-refined static prototype

This is a static, browser-based prototype calculator for pediatric eGFR estimation with optional transportability-guided center intercept calibration.

## Intended users

This prototype is intended for clinicians or researchers for research display and external validation support. Patients or caregivers should not upload patient-level data by themselves.

## Patient eligibility and input ranges

This public prototype is restricted to:

- Age: 0–18 years
- Height: 63–189 cm
- Serum creatinine: 9–273.1 μmol/L, equivalent to 0.10–3.09 mg/dL

Values outside these ranges are unavailable in this prototype. Ages below 0.33 years are allowed by the public age range but fall below the lower spline boundary used by the model and are flagged.

The model was developed using data from five Chinese hospitals. The overall median measured GFR was 110.2 mL/min/1.73 m², indicating that the development population had largely preserved kidney function. Results should be interpreted cautiously in populations with severe kidney impairment or clinical profiles outside the model development range.

## Files

- `index.html` — webpage structure
- `style.css` — visual styling
- `app.js` — eGFR equation, natural spline, unit conversion, and calibration logic
- `example_calibration.csv` — template for local mGFR calibration samples
- `example_batch_base_egfr.csv` — template for batch base eGFR calculation
- `LICENSE` — MIT license placeholder for open distribution

## Local use

Open `index.html` directly in Chrome or Edge. No server, R session, Python environment, database, or internet connection is required for local use.

For smoother development, open the folder in VS Code and use the Live Server extension.

## Serum creatinine unit conversion

The deployed equation uses serum creatinine in μmol/L. If the user enters serum creatinine in mg/dL, the calculator converts it internally using:

```text
Scr_μmol/L = Scr_mg/dL × 88.4
```

The same conversion is available for uploaded calibration and batch CSV files.

## Base model

The base equation implemented in `app.js` is:

```text
ln(mGFR_hat_base) = 4.283243
                  + 0.184304 * ln(height / Scr)
                  + 0.046885 * NS1(age)
                  + 0.336329 * NS2(age)
                  + 0.046227 * NS3(age)
```

Natural spline configuration:

```r
splines::ns(age, knots = c(5, 10), Boundary.knots = c(0.33, 20))
```

Units:

- age: years
- height: cm
- Scr: μmol/L after internal conversion
- mGFR/eGFR: mL/min/1.73 m²

## Local calibration

Upload a CSV with these columns:

```csv
age,height,Scr,mGFR
```

The public prototype applies calibration when at least 10 valid local samples are uploaded.

The calculator estimates:

```text
w_hat_naive = mean(log(mGFR) - log(mGFR_hat_base))
```

For 10 ≤ n < 20, empirical-Bayes shrinkage is applied:

```text
lambda = tau2 / (tau2 + sigma2 / n)
w_hat_used = lambda * w_hat_naive
```

Fixed EB parameters:

```text
tau2 = 0.01007562
sigma2 = 0.0345
```

For n ≥ 20, ordinary intercept-only calibration is used:

```text
w_hat_used = w_hat_naive
```

The calibrated estimate is:

```text
eGFR_calibrated = eGFR_base * exp(w_hat_used)
```

## Batch base eGFR calculation

Upload a CSV with these columns:

```csv
id,age,height,Scr
```

The `id` column is optional. The batch module calculates uncalibrated base eGFR only and allows the calculated results to be downloaded as a CSV file. Center-calibrated eGFR requires local mGFR calibration samples in the local calibration module.

## Recommendation rule in this prototype

This first public prototype uses center-level intercept shift magnitude:

- `|w_hat_naive| < 0.03`: calibration not recommended
- `0.03 ≤ |w_hat_naive| < 0.08`: review before routine calibration
- `|w_hat_naive| ≥ 0.08`: center calibration may be considered

Performance gains such as ΔRMSE, ΔP30, and ΔP10 require an independent local validation set. Clinicians may apply calibration according to the Recommendation panel. Researchers with an independent local validation set should use the supplementary article code to evaluate ΔRMSE, ΔP30, and ΔP10 before routine implementation.

## Important validation note

Before public release, verify that JavaScript `window.TGCCF.rCompatibleNsAge(age)` matches the R output of:

```r
splines::ns(
  c(0, 0.33, 1, 5, 10, 18),
  knots = c(5, 10),
  Boundary.knots = c(0.33, 20)
)
```

## Disclaimer

This prototype is intended for research display and external validation support only. It should not replace clinical judgment or be used as a standalone medical decision-making tool.

Patients or caregivers should not upload patient-level data by themselves. The tool is intended for use by clinicians or researchers.

All calculations are performed locally in the user's browser. This static webpage does not upload, store, or transmit patient-level data.
