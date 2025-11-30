# Monkeypox SEIR Modeling Sandbox

This workspace contains a coupled human-rodent SEIR-style model tuned to the 2022-2025 mpox (monkeypox) outbreaks. The C++ driver (`monkeypox_seir.cpp`) ingests observed case counts, optimizes the human-to-human transmission rate over detected outbreak windows, and exports both CSV summaries and a browser-based visualization.

## Repository Layout

- `monkeypox_seir.cpp` — core simulation, outbreak detection, and reporting logic.
- `observed_data.csv` — historical mpox case counts (daily). The code expects `Date,Infectious` columns in `MM-DD-YYYY` order.
- `monkeypox_fitted_prediction.csv` — generated each run with fitted daily totals (`Date,Time_Days,Infectious_Predicted`).
- `monkeypox_prediction.csv` — optional placeholder for scenario runs.
- `index.html` — Chart.js dashboard built by `generateHTMLReport` for side-by-side observed vs. simulated curves.

## Model Overview

- **Compartments (humans):** Susceptible, Exposed, Infectious (undetected), Quarantined/Detected, Recovered.
- **Compartments (rodents):** Susceptible, Exposed, Infectious reservoir.
- **Transmission:** Human-human (`beta2`) optimized per outbreak, rodent-human spillover (`beta1`) fixed, rodent self-transmission (`beta3`).
- **Control structure:** `beta_schedule` map lets the optimizer assign a distinct `beta2` at every detected outbreak boundary, enabling sequential fitting without re-running the whole timeline.

### Mathematical References

Model structure and baseline parameter magnitudes follow the extended SEIR framework in Md. Rasel et al., *Symmetry* 2022 (Symmetry 14, 2545) [[link](https://www.mdpi.com/2073-8994/14/12/2545#symmetry-14-02545-t001)]. Clinical duration choices (incubation, infectious, isolation) align with Victoria Department of Health mpox guidance [[link](https://www.health.vic.gov.au/infectious-diseases/mpox-monkeypox)]. See inline comments in `monkeypox_seir.cpp` for parameter-level citations.

## Data Preparation

1. Place observed daily counts in `observed_data.csv` sorted by date (CDC mpox surveillance feed: https://www.cdc.gov/monkeypox/data-research/cases/index.html).
2. Ensure dates use `MM-DD-YYYY` strings; the parser auto-computes day offsets from the first entry.
3. (Optional) Provide custom prediction horizons in `monkeypox_prediction.csv` if you plan to extend the tool beyond fitting.

## Building and Running

All commands assume Windows PowerShell in the workspace root:

```powershell
cd d:\approximationthingies
if ($?) { g++ monkeypox_seir.cpp -std=c++20 -O2 -o monkeypox_seir }
if ($?) { .\monkeypox_seir }
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
   ```bash
   git add index.html monkeypox_fitted_prediction.csv observed_data.csv
   git commit -m "Add SEIR visualization"
   git push origin main
   ```

2. **Enable GitHub Pages**:
   - Go to your repository **Settings** → **Pages**
   - Under "Source", select **Deploy from a branch**
   - Select branch: **main** and folder: **/ (root)**
   - Click **Save**

3. **Access your visualization**:
   - Your site will be available at: `https://[username].github.io/[repository-name]/`
   - Example: `https://yourusername.github.io/approximationthingies/`

4. **Note**: Make sure all CSV files are committed to the repository. GitHub Pages serves static files only, so the visualization loads data via fetch() from the same directory.

## Output Artifacts

- `monkeypox_fitted_prediction.csv`: deterministic daily totals from the stitched simulations.
- `index.html`: interactive chart with observed counts, predicted trajectory, and moving average.
- Console logs: report outbreak detections, optimal `beta2` per period, and diagnostics.

## Extending the Model

- Adjust `params.importation_rate` or the rodent reservoir state to explore recurrent seeding.
- Modify `detectOutbreakBoundaries` criteria or provide manual `beta_schedule` entries for scenario testing.
- Replace the coarse-to-fine grid search with a numerical optimizer if gradient-free accuracy becomes limiting.

## References

1. Md. Rasel, M. S. H. Hawlader, et al., "A Mathematical Model for Monkeypox Transmission," *Symmetry* 14(12):2545, 2022. https://www.mdpi.com/2073-8994/14/12/2545#symmetry-14-02545-t001
2. Victoria Department of Health, "Mpox (monkeypox)," 2023 guidance. https://www.health.vic.gov.au/infectious-diseases/mpox-monkeypox
3. Centers for Disease Control and Prevention, "Mpox (monkeypox) Data and Statistics." https://www.cdc.gov/monkeypox/data-research/cases/index.html

These sources ground the compartment choices, time scales, and intervention assumptions documented in the code.
