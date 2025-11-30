# GitHub Pages Deployment Checklist

## Files Required for Visualization

Make sure these files are committed and pushed:

- ✅ `index.html` - The visualization webpage
- ✅ `monkeypox_fitted_prediction.csv` - Main data file (106 KB)
- ✅ `observed_data.csv` - Original observed data (18 KB)
- ✅ `README.md` - Documentation

## Deployment Steps

1. **Check that CSV files are committed:**
   ```bash
   git status
   git add index.html monkeypox_fitted_prediction.csv observed_data.csv
   git commit -m "Add visualization and data files"
   git push origin main
   ```

2. **Enable GitHub Pages:**
   - Go to repository Settings → Pages
   - Source: Deploy from branch
   - Branch: main
   - Folder: / (root)
   - Click Save

3. **Wait 1-2 minutes** for deployment to complete

4. **Access your site:**
   - URL will be: `https://[username].github.io/[repository-name]/`
   - Check "Actions" tab to see deployment status

## Troubleshooting

### If data doesn't load on GitHub Pages:

1. **Check browser console (F12)** for specific error messages

2. **Verify files are in repository:**
   - Go to your repo on GitHub.com
   - Make sure you can see `monkeypox_fitted_prediction.csv` in the file list
   - Click on it to verify it has content (should show ~1288 lines)

3. **Check file paths:**
   - CSV files MUST be in the same directory as index.html
   - Path is case-sensitive on GitHub Pages (use exact filename)

4. **Verify CSV format:**
   - Open `monkeypox_fitted_prediction.csv` in a text editor
   - First line should be: `Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h,Observed,MA7,MA21`
   - Data lines should have 13 comma-separated values

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
- CSV file might be empty or corrupted
- Re-run C++ program to regenerate CSV
- Check that CSV has 13 columns per row
