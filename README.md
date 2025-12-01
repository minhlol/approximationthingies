# Monkeypox Two-Population SEIQR-SEIR Model

This workspace contains a coupled human-rodent compartmental model for the 2022-2025 mpox (monkeypox) outbreaks. The C++ implementation uses a **three-stage hybrid approach**: (1) moving average smoothing, (2) human compartment fitting, and (3) rodent SEIR dynamics refinement. It generates both CSV outputs and an interactive browser-based visualization with 21-day rolling forecasts.

## Repository Layout

- `monkeypox_seir.cpp` — core simulation with three-stage modeling, outbreak detection, and 21-day rolling forecasts
- `observed_data.csv` — historical mpox case counts (daily), expects `Date,Infectious` columns in `MM-DD-YYYY` format
- `monkeypox_fitted_prediction.csv` — generated output with 14 columns including S_h, E_h, I_h, Q_h, R_h, S_r, E_r, I_r, moving averages, and 21-day forecasts
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

**Rodent compartments (SEIR model):**
- S_r: Susceptible rodents
- E_r: Exposed rodents
- I_r: Infectious rodents (reservoir)

**Transmission dynamics:**
- Zoonotic spillover: β₁·I_r (rodent-to-human)
- Human-to-human: β₂·I_h
- Rodent intra-species: β₃·I_r
- Dual exposure pathways: E → I (undetected) and E → Q (detected/quarantined)

### Three-Stage Approach

1. **Stage 1 - Moving Average Smoothing:** Exponential smoothing (α=0.25) of 21-day moving average to create stable epidemic trajectory
2. **Stage 2 - Human Compartment Construction:** Build S_h, E_h, I_h, Q_h, R_h from smoothed prediction using compartmental flow equations
3. **Stage 3 - Rodent SEIR Refinement:** Integrate full rodent dynamics with population homeostasis (immigration) and endemic maintenance (replenishment when I_r drops below 50% target)

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

4. **Note**: Make sure all CSV files are committed to the repository. GitHub Pages serves static files only, so the visualization loads data via fetch() from the same directory.

## Forecasting System

### 21-Day Rolling Forecasts

The model generates **trend-aware 21-day forecasts** for every historical point:

- **Growth/Stable Phases** (trend > -2 cases/day): Uses 90% model trajectory, capturing epidemic momentum
- **Rapid Decline** (trend < -2 cases/day): Uses full SEIQR-SEIR dynamics with homeostasis
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
  - Linear/logarithmic scale toggle
  - Summary statistics cards
- Console logs: Outbreak detection, RMSE per period, rodent dynamics, forecast generation progress

## Model Parameters

Key epidemiological parameters (daily rates):
- **Human transmission**: β₂ = 0.00012 (fitted to outbreak dynamics)
- **Zoonotic spillover**: β₁ = 0.00008
- **Incubation rates**: α₁ = 0.143 (E→I), α₂ = 0.067 (E→Q)
- **Recovery/isolation**: ν = 0.067 (I→R), τ = 0.143 (Q→R)
- **Rodent dynamics**: β₃ = 0.0003, α₃ = 0.20, target I_r = 1.5% of population

Populations:
- N_h = 1,000,000 (human)
- N_r = 100,000 (rodent reservoir)

## Extending the Model

- **Modify forecast horizon**: Change `forecast_days` from 21 to other values (7, 14, 30, 90 days)
- **Adjust rodent maintenance**: Modify `target_I_r_forecast` and replenishment rates to explore endemic persistence
- **Parameter sensitivity**: Test different β values to model intervention scenarios
- **Outbreak detection**: Adjust MA21 threshold (currently 5 cases) for different epidemic definitions
- **Smoothing factors**: Change exponential smoothing α (currently 0.25) for more/less responsiveness

## References

1. Md. Rasel, M. S. H. Hawlader, et al., "A Mathematical Model for Monkeypox Transmission," *Symmetry* 14(12):2545, 2022. https://www.mdpi.com/2073-8994/14/12/2545#symmetry-14-02545-t001
2. Victoria Department of Health, "Mpox (monkeypox)," 2023 guidance. https://www.health.vic.gov.au/infectious-diseases/mpox-monkeypox
3. Centers for Disease Control and Prevention, "Mpox (monkeypox) Data and Statistics." https://www.cdc.gov/monkeypox/data-research/cases/index.html

These sources ground the compartment choices, time scales, and intervention assumptions documented in the code.
