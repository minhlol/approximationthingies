# GitHub Pages Deployment Checklist

## Files Required for Visualization

Make sure these files are committed and pushed:

- ✅ `index.html` - Interactive visualization webpage with Chart.js
- ✅ `styles.css` - Styling for the visualization
- ✅ `monkeypox_seir.cpp` - C++ source code for SEIQR-SEIR model (optional for Pages, but good for reference)
- ✅ `monkeypox_fitted_prediction.csv` - Main data file with 14 columns (~1289 lines, includes 21-day forecasts)
- ✅ `monkeypox_prediction.csv` - 21-day future forecast (91 lines)
- ✅ `observed_data.csv` - Original observed case data (1288 points)
- ✅ `README.md` - Complete documentation
- ✅ `ARCHITECTURE.md` - Model architecture and flow diagrams
- ✅ `DEPLOYMENT.md` - This deployment guide

## Deployment Steps

1. **Check that all required files are committed:**
   ```powershell
   git status
   git add index.html styles.css monkeypox_fitted_prediction.csv monkeypox_prediction.csv observed_data.csv README.md DEPLOYMENT.md ARCHITECTURE.md
   git commit -m "Update SEIQR-SEIR model with 21-day rolling forecasts"
   git push origin main
   ```

2. **Enable GitHub Pages:**
   - Go to repository Settings → Pages
   - Source: Deploy from branch
   - Branch: main
   - Folder: / (root)
   - Click Save

3. **Wait 1-2 minutes** for deployment to complete

3. **Access your site:**
   - URL: `https://minhlol.github.io/approximationthingies/`
   - Check "Actions" tab to see deployment status
   - Deployment typically completes in 1-2 minutes

## Troubleshooting

### If data doesn't load on GitHub Pages:

1. **Check browser console (F12)** for specific error messages

2. **Verify files are in repository:**
   - Go to your repo on GitHub.com
   - Make sure you can see both `monkeypox_fitted_prediction.csv` and `monkeypox_prediction.csv`
   - Click on `monkeypox_fitted_prediction.csv` to verify it has ~1289 lines
   - Click on `monkeypox_prediction.csv` to verify it has ~91 lines

3. **Check file paths:**
   - CSV files MUST be in the same directory as index.html (root of repository)
   - Path is case-sensitive on GitHub Pages (use exact filename)
   - Both historical and forecast CSVs are loaded by index.html

4. **Verify CSV format:**
   - Open `monkeypox_fitted_prediction.csv` in a text editor
   - First line should be: `Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h,Observed,MA7,MA21,Forecast_90d`
   - Data lines should have 14 comma-separated values (including 21-day forecast column)
   - Open `monkeypox_prediction.csv` and verify it has: `Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h`

5. **Check .gitignore:**
   - Make sure `.gitignore` does NOT exclude `*.csv` files
   - Our .gitignore only excludes `.exe`, `.o`, `.obj`, `.bak`, and `.vscode/`

6. **Clear browser cache:**
   - Hard refresh: Ctrl+F5 (Windows) or Cmd+Shift+R (Mac)
   - Or open in incognito/private window

## Testing Locally

Before deploying, test locally:

```powershell
# Start local server
.\server.ps1

# Open browser to http://localhost:8080/
# You should see the visualization with data
```

## Common Issues

**404 Error on CSV files:**
- Files not committed to git
- Files in wrong directory
- Case mismatch in filename

**CORS Error:**
- Shouldn't happen on GitHub Pages (same origin)
- If testing locally, must use a server (not file://)

**Blank Chart:**
- Check browser console for JavaScript errors
- Verify CSV has correct format
- Check that Chart.js CDN is loading

**Data Shows but is All Zeros:**
- CSV files might be empty or corrupted
- Re-run C++ program to regenerate both CSV files
- Check that fitted prediction has 14 columns and forecast has 10 columns per row

**Forecast Line Not Showing:**
- Verify `Forecast_90d` column exists in `monkeypox_fitted_prediction.csv` (14th column)
- Check that `monkeypox_prediction.csv` is loaded (check browser console)
- Forecast line is purple/dashed, displayed 21 days ahead of generation point
- When filtering by year, forecasts from previous year may predict into selected year

## Visualization Features

Once deployed, your visualization includes:

- **Three View Modes:**
  - Comparison: Observed vs Model vs 21-Day Forecast
  - Human Compartments: S_h, E_h, I_h, Q_h, R_h dynamics
  - Rodent Compartments: S_r, E_r, I_r reservoir dynamics

- **Year Filtering:** 2022, 2023, 2024, 2025, or All Years

- **Scale Toggle:** Linear or Logarithmic Y-axis

- **Summary Statistics:** Total cases, peak values, final compartment states

- **Interactive Legend:** Click to show/hide individual data series
