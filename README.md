# Monkeypox Two-Population SEIQR-SEI Model

This workspace contains a coupled human-rodent compartmental model for the 2022-2025 mpox (monkeypox) outbreaks. The C++ implementation uses **continuous correction with dynamic parameter adjustment** to fit observed data while maintaining smooth epidemiological trajectories. It features a full two-population SEIQR-SEI system with endemic rodent reservoir dynamics and generates 21-day trend-based rolling forecasts.

## Repository Layout

- `monkeypox_seir.cpp` — core simulation with RK4 integration, continuous correction, dynamic parameter adjustment, and 21-day trend-based forecasts
- `observed_data.csv` — historical mpox case counts (daily), expects `Date,Infectious` columns in `MM-DD-YYYY` format
- `monkeypox_fitted_prediction.csv` — generated output with 13 columns including S_h, E_h, I_h, Q_h, R_h, S_r, E_r, I_r, Total_Infected_h, Observed, MA7, MA21, and Forecast_21d
- `monkeypox_prediction.csv` — 21-day future forecast beyond last observed date
- `index.html` — Interactive Chart.js dashboard with year filtering, multiple views, and forecast visualization

## Model Overview

### Compartmental Structure

**Human compartments (SEIQR model):**
- S_h: Susceptible
- E_h: Exposed (incubating)
- I_h: Infectious (undetected/community transmission)
- Q_h: Quarantined/Detected (isolated cases)
- R_h: Recovered (immune)

**Rodent compartments (SEI model):**
- S_r: Susceptible rodents (endemic reservoir)
- E_r: Exposed rodents
- I_r: Infectious rodents (maintains ~8 individuals at equilibrium)

**Transmission dynamics:**
- Zoonotic spillover: β₁·I_r (rodent-to-human)
- Human-to-human: β₂·I_h
- Rodent intra-species: β₃·I_r
- Dual exposure pathways: E → I (undetected) and E → Q (detected/quarantined)

### Model Approach

**Continuous Correction with Period-Specific Parameters:**
- **RK4 Integration:** Fourth-order Runge-Kutta integration of full SEIQR-SEI differential equations
- **Dynamic β Adjustment:** Transmission rates (β₁, β₂) automatically calibrated per outbreak period based on observed growth/decline
- **Continuous Correction:** 25% daily correction rate targeting 4× observed values to account for underreporting
- **Cubic Spline Smoothing:** Dual-pass smoothing (window-5 + window-3) applied to fitted model only
- **Endemic Rodent Dynamics:** Rodent population maintains natural equilibrium through homeostasis (β₃=0.08, endemic I_r≈8)

### Mathematical References

Model structure and baseline parameter magnitudes follow the extended SEIR framework in Md. Rasel et al., *Symmetry* 2022 (Symmetry 14, 2545) [[link](https://www.mdpi.com/2073-8994/14/12/2545#symmetry-14-02545-t001)]. Clinical duration choices (incubation, infectious, isolation) align with Victoria Department of Health mpox guidance [[link](https://www.health.vic.gov.au/infectious-diseases/mpox-monkeypox)]. See inline comments in `monkeypox_seir.cpp` for parameter-level citations.

## Data Preparation

1. Place observed daily counts in `observed_data.csv` sorted by date (CDC mpox surveillance feed: https://www.cdc.gov/monkeypox/data-research/cases/index.html).
2. Ensure dates use `MM-DD-YYYY` strings; the parser auto-computes day offsets from the first entry.
3. (Optional) Provide custom prediction horizons in `monkeypox_prediction.csv` if you plan to extend the tool beyond fitting.

## Building and Running

All commands assume Windows PowerShell in the workspace root:

```powershell
# Navigate to project directory
cd c:\Users\Admin\Documents\GitHub\approximationthingies
# Compile the C++ model
g++ monkeypox_seir.cpp -std=c++20 -O2 -o monkeypox_seir
# Run the simulation
.\monkeypox_seir
```

Execution steps:

1. Compile with a modern C++ compiler (tested with `g++`/MinGW, C++20).
2. Run the generated `monkeypox_seir.exe` to:
   - Load observed data and compute 7-day and 21-day moving averages.
   - Detect outbreak periods where 21-day MA exceeds threshold.
   - Generate fitted CSV outputs with moving averages.
3. Open `index.html` in any browser to inspect results (use controls to filter by year, change view, and toggle scale).

## GitHub Pages Deployment

To deploy the visualization on GitHub Pages:

1. **Commit the required files** to your repository:
   ```powershell
   git add index.html styles.css monkeypox_fitted_prediction.csv monkeypox_prediction.csv observed_data.csv
   git commit -m "Add SEIQR-SEIR visualization with 21-day forecasts"
   git push origin main
   ```

2. **Enable GitHub Pages**:
   - Go to your repository **Settings** → **Pages**
   - Under "Source", select **Deploy from a branch**
   - Select branch: **main** and folder: **/ (root)**
   - Click **Save**

3. **Access your visualization**:
   - Your site will be available at: `https://minhlol.github.io/approximationthingies/`
   - GitHub Pages deployment typically takes 1-2 minutes
## Forecasting System

### 21-Day Rolling Forecasts

The model generates **trend-based 21-day forecasts** for every historical point using the smoothed fitted model:

- **Forecast Method:** Weighted average of 70% future fitted value (21 days ahead) + 30% average trajectory over 21-day window
- **Smooth by Design:** Uses already-smoothed fitted model values, producing naturally smooth forecast curves
- **Date-Aware Display:** Forecast values displayed at their predicted date (21 days ahead), not generation date
- **No Additional Correction:** Forecasts inherit the correction already applied during model fitting

Forecasts allow retrospective validation by comparing predicted values with actual observations 21 days later.
- **Conditional Dynamics**: When I_r ≈ 0, switches to simplified human-only model (no zoonotic spillover)

Forecasts are displayed 21 days ahead of their generation point, allowing comparison with actual observations.

### Key Features

- **Date-aware shifting**: Forecast values displayed at their predicted date, not generation date
- **Year filtering**: Properly handles cross-year forecasts (e.g., late 2022 forecasts into early 2023)
- **Epidemic phase detection**: Automatically adapts to growth vs. decline dynamics

## Output Artifacts

- `monkeypox_fitted_prediction.csv`: 14 columns including all compartments, moving averages (7-day, 21-day), and 21-day rolling forecasts
- `monkeypox_prediction.csv`: 21-day future forecast extending beyond last observed date
- `index.html`: Interactive visualization with:
  - Three view modes: Comparison (observed vs model vs forecast), Human compartments, Rodent compartments
  - Year filtering (2022-2025)
## Model Parameters

Key epidemiological parameters (daily rates):
- **Human transmission**: β₂ = 0.00035 (baseline), dynamically adjusted per outbreak (range: 0.36-5.50)
- **Zoonotic spillover**: β₁ = 0.0002 (baseline), dynamically adjusted (coupled to β₂)
- **Incubation rates**: α₁ = 0.2 (E→I, ~5 days), α₂ = 2.0 (E→Q)
- **Recovery/isolation**: ν = 0.83 (I→R, ~14 days), τ = 0.52 (Q→R, ~15 days)
- **Removal rate**: 4.607/day (combined α₁+α₂+ν+δ_h+μ_h)
- **Rodent dynamics**: β₃ = 0.08, α₃ = 2.0, δ_r = 0.4, endemic equilibrium: E_r≈3, I_r≈8

Populations:
- N_h = 1,000,000 (human)
- N_r = 100,000 (rodent reservoir)
- Initial: E_h=10, I_h=5, Q_h=1; E_r=2000, I_r=2000, α₃ = 0.20, target I_r = 1.5% of population

Populations:
- N_h = 1,000,000 (human)
- N_r = 100,000 (rodent reservoir)

## Extending the Model

- **Modify forecast horizon**: Change `forecast_days` from 21 to other values (7, 14, 30, 60 days)
- **Adjust correction rate**: Modify `correction_rate` (currently 0.25) or target multiplier (currently 4×) to change fit tightness
- **Rodent dynamics**: Adjust β₃ (0.08) or δ_r (0.4) to explore different endemic equilibrium levels
- **Parameter sensitivity**: Test different baseline β₁/β₂ values for intervention scenarios
- **Outbreak detection**: Adjust MA21 threshold (currently 5 cases) for different epidemic definitions
- **Smoothing strength**: Change cubic spline window sizes (currently 5 and 3) for more/less smoothness

## References

1. Md. Rasel, M. S. H. Hawlader, et al., "A Mathematical Model for Monkeypox Transmission," *Symmetry* 14(12):2545, 2022. https://www.mdpi.com/2073-8994/14/12/2545#symmetry-14-02545-t001
2. Victoria Department of Health, "Mpox (monkeypox)," 2023 guidance. https://www.health.vic.gov.au/infectious-diseases/mpox-monkeypox
3. Centers for Disease Control and Prevention, "Mpox (monkeypox) Data and Statistics." https://www.cdc.gov/monkeypox/data-research/cases/index.html

These sources ground the compartment choices, time scales, and intervention assumptions documented in the code.
