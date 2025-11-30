#include <bits/stdc++.h>
using namespace std;

// SEIR Model Parameters (Extended Human-Rodent Model)
// Reference models: Md. Rasel et al., Symmetry 2022, doi:10.3390/sym14122545 (Table 1)
// Additional mpox clinical/literature context: Victoria Department of Health, mpox guidance
struct Parameters
{
    // Human Parameters (calibrated against Symmetry 2022 model assumptions)
    double theta_h; // Human recruitment rate
    double beta1;   // Human-rodent contact rate
    double beta2;   // Human-human contact rate
    double alpha1;  // Proportion of infected humans from exposed humans (E_h -> I_h)
    double alpha2;  // Proportion of suspected cases detected (E_h -> Q_h)
    double phi;     // Proportion not detected after diagnosis (Q_h -> S_h)
    double tau;     // Progression from isolated to recovered (Q_h -> R_h)
    double nu;      // Humans recovery rate (I_h -> R_h)
    double mu_h;    // Natural human death rate
    double delta_h; // Disease induced human death rate

    // Rodent Parameters (consistent with the reservoir compartment in Symmetry 2022)
    double theta_r; // Rodent recruitment rate
    double beta3;   // Rodent-rodent contact rate
    double alpha3;  // Proportion of infected rodents from exposed rodents (E_r -> I_r)
    double mu_r;    // Natural rodent death rate
    double delta_r; // Disease induced rodent death rate

    // Multiple Outbreak Parameters (scenario-based control input)
    double importation_rate;           // Daily external infections introduced (travelers etc.)
    map<double, double> beta_schedule; // Map of (Day -> Beta2 Value) to change contact rate over time
};

// State variables
struct State
{
    // Human Compartments
    double S_h; // Susceptible
    double E_h; // Exposed
    double I_h; // Infectious (Undetected)
    double Q_h; // Quarantined/Detected
    double R_h; // Recovered

    // Rodent Compartments
    double S_r; // Susceptible
    double E_r; // Exposed
    double I_r; // Infectious

    double t; // Time (days from start)
};

struct Observation
{
    string dateStr;
    double time; // Days from start
    double infectious;
    double movingAvg7Day;
};

struct DailyPrediction
{
    string dateStr;
    double time;
    double infectious;
};

// Date Helper Functions
time_t parseDate(const string &dateStr)
{
    tm tm = {};
    stringstream ss(dateStr);
    char delimiter;
    int year, month, day;
    // Expected format: MM-DD-YYYY
    ss >> month >> delimiter >> day >> delimiter >> year;
    tm.tm_year = year - 1900;
    tm.tm_mon = month - 1;
    tm.tm_mday = day;
    tm.tm_hour = 12; // Noon to avoid DST issues
    return mktime(&tm);
}

string formatDate(time_t date)
{
    tm *tm = localtime(&date);
    stringstream ss;
    ss << setw(2) << setfill('0') << (tm->tm_mon + 1) << "-"
       << setw(2) << setfill('0') << tm->tm_mday << "-"
       << (tm->tm_year + 1900);
    return ss.str();
}

double daysBetween(time_t start, time_t end)
{
    return difftime(end, start) / (60 * 60 * 24);
}

time_t addDays(time_t start, double days)
{
    return start + (time_t)(days * 24 * 60 * 60);
}

// Function to read observed data from CSV
vector<Observation> readObservations(const string &filename)
{
    vector<Observation> data;
    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "Error: Could not open file " << filename << endl;
        return data;
    }

    string line;

    // Skip header
    getline(file, line);

    time_t startDate = 0;
    bool first = true;

    while (getline(file, line))
    {
        // Handle both comma and tab delimiters by replacing them with spaces
        replace(line.begin(), line.end(), ',', ' ');
        replace(line.begin(), line.end(), '\t', ' ');

        stringstream ss(line);
        string dateVal;
        double infectiousVal;

        if (ss >> dateVal >> infectiousVal)
        {
            Observation obs;
            obs.dateStr = dateVal;
            obs.infectious = infectiousVal;

            time_t currentDate = parseDate(dateVal);

            if (first)
            {
                startDate = currentDate;
                obs.time = 0.0;
                first = false;
            }
            else
            {
                obs.time = daysBetween(startDate, currentDate);
            }
            data.push_back(obs);
        }
    }
    return data;
}

// Function to run simulation for a specific period with given initial state
vector<State> simulatePeriod(Parameters params, State initialState, double t_end, double dt)
{
    vector<State> history;
    State current = initialState;
    history.push_back(current);

    while (current.t < t_end - 1e-9)
    {
        double remaining = t_end - current.t;
        if (remaining <= 1e-9)
            break;
        double step = min(dt, remaining);

        // Calculate current populations
        double N_h = current.S_h + current.E_h + current.I_h + current.Q_h + current.R_h;
        double N_r = current.S_r + current.E_r + current.I_r;

        if (N_h < 1.0)
            N_h = 1.0;
        if (N_r < 1.0)
            N_r = 1.0;

        // Determine current Human-Human Transmission Rate (beta2) based on schedule
        double current_beta2 = params.beta2;
        if (!params.beta_schedule.empty())
        {
            auto it = params.beta_schedule.upper_bound(current.t);
            if (it != params.beta_schedule.begin())
            {
                current_beta2 = prev(it)->second;
            }
        }

        // Human Derivatives
        double rodent_spillover = params.beta1 * current.I_r / N_r;
        double human_spread = current_beta2 * current.I_h / N_h;
        double force_infection_h = rodent_spillover + human_spread;
        double dS_h = params.theta_h - force_infection_h * current.S_h - params.mu_h * current.S_h + params.phi * current.Q_h;
        double dE_h = force_infection_h * current.S_h - (params.alpha1 + params.alpha2 + params.mu_h) * current.E_h + params.importation_rate;
        double dI_h = params.alpha1 * current.E_h - (params.mu_h + params.delta_h + params.nu) * current.I_h;
        double dQ_h = params.alpha2 * current.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * current.Q_h;
        double dR_h = params.nu * current.I_h + params.tau * current.Q_h - params.mu_h * current.R_h;

        // Rodent Derivatives
        double force_infection_r = (params.beta3 * current.I_r) / N_r;
        double dS_r = params.theta_r - force_infection_r * current.S_r - params.mu_r * current.S_r;
        double dE_r = force_infection_r * current.S_r - (params.mu_r + params.alpha3) * current.E_r;
        double dI_r = params.alpha3 * current.E_r - (params.mu_r + params.delta_r) * current.I_r;

        // Update State
        current.S_h += dS_h * step;
        current.E_h += dE_h * step;
        current.I_h += dI_h * step;
        current.Q_h += dQ_h * step;
        current.R_h += dR_h * step;

        current.S_r += dS_r * step;
        current.E_r += dE_r * step;
        current.I_r += dI_r * step;

        current.t += step;
        history.push_back(current);
    }
    return history;
}

// Calculate Sum of Squared Errors (SSE)
double calculateSSE(const vector<Observation> &observed, const vector<State> &predicted)
{
    double sse = 0.0;

    for (const auto &obs : observed)
    {
        // Find closest predicted point in time
        auto it = lower_bound(predicted.begin(), predicted.end(), obs.time,
                              [](const State &s, double t)
                              {
                                  return s.t < t;
                              });

        if (it != predicted.end())
        {
            // Compare observed infectious against Total Active Human Infections (I_h + Q_h)
            // Assuming observed data captures both or represents the active burden
            double total_infectious_h = it->I_h + it->Q_h;
            double error = obs.infectious - total_infectious_h;
            sse += error * error;
        }
    }
    return sse;
}

// Function to generate HTML Visualization
void generateHTMLReport(const vector<Observation> &observed, const vector<DailyPrediction> &predicted)
{
    ofstream htmlFile("index.html");

    htmlFile << "<!DOCTYPE html>\n<html>\n<head>\n<title>Monkeypox SEIR Prediction</title>\n";
    htmlFile << "<script src='https://cdn.jsdelivr.net/npm/chart.js'></script>\n";
    htmlFile << "<style>body { font-family: sans-serif; } .controls { text-align: center; margin-bottom: 20px; }</style>\n";
    htmlFile << "</head>\n<body>\n";
    htmlFile << "<div style='width: 80%; margin: auto;'>\n";
    htmlFile << "<h2 style='text-align: center;'>Monkeypox Outbreak Prediction (SEIR Model)</h2>\n";

    // Add Year Selector
    htmlFile << "<div class='controls'>\n";
    htmlFile << "<label for='yearSelect'>Select Year: </label>\n";
    htmlFile << "<select id='yearSelect' onchange='updateChart()'>\n";
    htmlFile << "<option value='all'>All Years</option>\n";
    htmlFile << "<option value='2022'>2022</option>\n";
    htmlFile << "<option value='2023'>2023</option>\n";
    htmlFile << "<option value='2024'>2024</option>\n";
    htmlFile << "<option value='2025'>2025</option>\n";
    htmlFile << "</select>\n";
    htmlFile << "</div>\n";

    htmlFile << "<canvas id='seirChart'></canvas>\n</div>\n";
    htmlFile << "<script>\n";

    // Generate Labels (Dates from prediction)
    htmlFile << "const allLabels = [";
    for (const auto &pred : predicted)
    {
        htmlFile << "'" << pred.dateStr << "',";
    }
    htmlFile << "];\n";

    // Generate Predicted Data
    htmlFile << "const allPredictedData = [";
    for (const auto &pred : predicted)
    {
        htmlFile << pred.infectious << ",";
    }
    htmlFile << "];\n";

    // Generate 7-Day Moving Average Data (Observed)
    htmlFile << "const allMovingAvg = [";
    size_t maIdx = 0;
    for (const auto &pred : predicted)
    {
        // Advance index if observed time is behind current prediction time
        while (maIdx < observed.size() && observed[maIdx].time < pred.time - 0.5)
        {
            maIdx++;
        }

        bool found = false;
        if (maIdx < observed.size())
        {
            double timeDiff = fabs(observed[maIdx].time - pred.time);
            if (timeDiff < 0.5)
            {
                htmlFile << observed[maIdx].movingAvg7Day << ",";
                found = true;
            }
        }

        if (!found)
        {
            htmlFile << "null,";
        }
    }
    htmlFile << "];\n";

    // Generate Observed Data (Aligned)
    htmlFile << "const allObservedData = [";
    size_t obsIdx = 0;
    for (const auto &pred : predicted)
    {
        // Advance index if observed time is behind current prediction time
        while (obsIdx < observed.size() && observed[obsIdx].time < pred.time - 0.5)
        {
            obsIdx++;
        }

        bool found = false;
        // Check if current observed data point matches the prediction time (within 0.5 days)
        if (obsIdx < observed.size())
        {
            double timeDiff = fabs(observed[obsIdx].time - pred.time);
            if (timeDiff < 0.5)
            {
                htmlFile << observed[obsIdx].infectious << ",";
                found = true;
            }
        }

        if (!found)
        {
            htmlFile << "null,";
        }
    }
    htmlFile << "];\n";

    // JavaScript for Chart and Filtering
    htmlFile << R"(
    const ctx = document.getElementById('seirChart').getContext('2d');
    let myChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: allLabels,
            datasets: [{
                label: 'Observed Cases',
                data: allObservedData,
                borderColor: 'red',
                backgroundColor: 'red',
                pointRadius: 3,
                pointHoverRadius: 5,
                showLine: false,
                spanGaps: false
            }, {
                label: 'Predicted Cases (SEIR)',
                data: allPredictedData,
                borderColor: 'blue',
                backgroundColor: 'rgba(0, 0, 255, 0.1)',
                fill: true,
                pointRadius: 0,
                borderWidth: 2
            }, {
                label: '7-Day Moving Average',
                data: allMovingAvg,
                borderColor: 'orange',
                backgroundColor: 'rgba(255, 165, 0, 0.1)',
                borderWidth: 2,
                pointRadius: 0,
                fill: false
            }]
        },
        options: {
            responsive: true,
            interaction: {
                mode: 'index',
                intersect: false,
            },
            scales: {
                x: { title: { display: true, text: 'Date' } },
                y: { title: { display: true, text: 'Infectious Cases' } }
            }
        }
    });

    function updateChart() {
        const year = document.getElementById('yearSelect').value;
        
        let filteredLabels = [];
        let filteredObserved = [];
        let filteredPredicted = [];
        let filteredMovingAvg = [];

        for (let i = 0; i < allLabels.length; i++) {
            if (year === 'all' || allLabels[i].includes(year)) {
                filteredLabels.push(allLabels[i]);
                filteredObserved.push(allObservedData[i]);
                filteredPredicted.push(allPredictedData[i]);
                filteredMovingAvg.push(allMovingAvg[i]);
            }
        }

        myChart.data.labels = filteredLabels;
        myChart.data.datasets[0].data = filteredObserved;
        myChart.data.datasets[1].data = filteredPredicted;
        myChart.data.datasets[2].data = filteredMovingAvg;
        myChart.update();
    }
    )";

    htmlFile << "</script>\n</body>\n</html>";

    htmlFile.close();
    cout << "Visualization saved to 'index.html'. Open this file in your browser." << endl;
}
// Function to detect outbreak start points based on 7-day Moving Average trends
vector<double> detectOutbreakBoundaries(const vector<Observation> &data)
{
    vector<double> boundaries;
    boundaries.push_back(0.0); // Start of simulation

    if (data.size() < 30)
        return boundaries;

    double last_boundary = 0.0;
    double min_gap = 60.0; // Minimum days between detected outbreaks

    // Scan for valleys in the 7-day Moving Average
    // We start a bit later to avoid initial boundary effects
    for (size_t i = 20; i < data.size() - 21; ++i)
    {
        if (data[i].time - last_boundary < min_gap)
            continue;

        double current_ma = data[i].movingAvg7Day;

        // 1. Check if this is a local minimum (or flat bottom) in a +/- 7 day window
        // This ensures we are at the bottom of a curve
        bool is_valley = true;
        for (int k = -7; k <= 7; ++k)
        {
            int idx = static_cast<int>(i) + k;
            if (idx < 0 || idx >= static_cast<int>(data.size()))
                continue;

            if (data[idx].movingAvg7Day < current_ma)
            {
                is_valley = false;
                break;
            }
        }

        if (is_valley)
        {
            // 2. Check for sustained increase in the NEXT 21 days
            // We want to see the MA rise significantly to confirm an outbreak
            // Bounds already checked by loop condition (i < data.size() - 21)
            double future_ma = data[i + 21].movingAvg7Day;

            // Criteria:
            // - Relative increase: > 20% (1.2x)
            // - Absolute increase: > 2 cases (to ignore noise at near-zero levels)
            if (future_ma > current_ma * 1.2 && (future_ma - current_ma) > 2.0)
            {

                double boundary_time = data[i].time;
                boundaries.push_back(boundary_time);
                last_boundary = boundary_time;

                cout << "Detected outbreak start (7-day MA Valley) at Day " << boundary_time
                     << " (" << data[i].dateStr << ") - MA: " << current_ma << " -> " << future_ma << endl;

                // Skip forward to avoid detecting the same valley multiple times
                i += 20;
            }
        }
    }

    return boundaries;
}

int main()
{
    // 1. Load observed data
    string dataFile = "observed_data.csv";
    vector<Observation> observedData = readObservations(dataFile);

    if (observedData.empty())
    {
        cerr << "Error: No data loaded from " << dataFile << endl;
        return 1;
    }

    cout << "Loaded " << observedData.size() << " data points." << endl;
    cout << "Start Date: " << observedData[0].dateStr << endl;

    // Calculate 7-Day Moving Average for Observed Data
    for (size_t i = 0; i < observedData.size(); ++i)
    {
        double sum = 0.0;
        int count = 0;
        for (int j = 0; j < 7; ++j)
        {
            if ((int)i - j >= 0)
            {
                sum += observedData[i - j].infectious;
                count++;
            }
        }
        observedData[i].movingAvg7Day = (count > 0) ? sum / count : 0.0;
    }

    time_t startDate = parseDate(observedData[0].dateStr);

    // 2. Setup fixed parameters
    Parameters params;

    // Human Parameters (incubation/infectious periods aligned with Symmetry 2022 Table 1 and Vic Health mpox notes)
    params.mu_h = 1.0 / (70.0 * 365.0);      // Natural death rate (70 year life expectancy)
    params.theta_h = params.mu_h * 100000.0; // Recruitment balances death for N=100k
    params.delta_h = 0.0;                    // Disease-induced death (CFR ~0% for 2022 clade II outbreak)

    params.alpha1 = 1.0 / 8.0;  // E_h -> I_h (Incubation period: 8 days, median for mpox)
    params.alpha2 = 1.0 / 10.0; // E_h -> Q_h (Detection/isolation delay: ~10 days)

    params.nu = 1.0 / 21.0;  // I_h -> R_h (Infectious period: ~21 days for mpox)
    params.tau = 1.0 / 14.0; // Q_h -> R_h (Recovery from isolation: ~14 days)
    params.phi = 0.001;      // Q_h -> S_h (Very small isolation escape/reinfection rate)

    params.beta1 = 0.0001; // Human-Rodent spillover contact rate (very low, minimal role)

    // Rodent Parameters (lifespan and transmission adapted from Symmetry 2022 reservoir settings)
    params.mu_r = 1.0 / (1.5 * 365.0);      // Rodent life ~1.5 years (small rodents)
    params.theta_r = params.mu_r * 10000.0; // Maintains stable rodent population
    params.delta_r = 0.03;                  // Disease-induced rodent death rate
    params.alpha3 = 1.0 / 7.0;              // Rodent incubation ~7 days
    params.beta3 = 0.01;                    // Rodent-rodent transmission (low endemic maintenance)

    // Multiple Outbreak Settings
    params.importation_rate = 0.0; // Set to > 0 to allow constant external seeding
    params.beta2 = 0.05;           // Initialize default beta2 value

    // 3. Sequential Optimization for each Year
    cout << "Optimizing Beta2 separately for each detected outbreak period..." << endl;

    // Define time boundaries based on detected outbreaks
    vector<double> boundaries = detectOutbreakBoundaries(observedData);
    boundaries.push_back(observedData.back().time + 1.0); // End of data

    map<double, double> optimized_schedule;

    // Initial State
    State currentState;
    currentState.t = 0;
    double N_h_initial = 10000.0;                        // Effective susceptible population in contact network
    currentState.I_h = observedData[0].infectious * 8.0; // Underreporting factor ~8x
    // Initialize with exposed population reflecting early transmission chains
    currentState.E_h = currentState.I_h * 12.0; // Large exposed cohort from early undetected spread
    currentState.Q_h = 0.0;
    currentState.R_h = 0.0;
    currentState.S_h = N_h_initial - currentState.I_h - currentState.E_h;

    // Rodent Initial (reservoir population)
    double N_r_initial = 10000.0;
    currentState.I_r = 10.0; // Endemic rodent infections (minimal stable reservoir)
    currentState.E_r = 5.0;  // Exposed rodents
    currentState.S_r = N_r_initial - currentState.I_r - currentState.E_r;

    vector<State> fullHistory;

    for (size_t i = 0; i < boundaries.size() - 1; ++i)
    {
        double t_start = boundaries[i];
        double t_end = boundaries[i + 1];

        // Re-initialize I_h based on observed data at t_start to prevent error propagation
        // SKIP reset for the very first period (i=0) as initial state is already properly set
        if (i > 0)
        {
            // Find observed data point closest to t_start
            auto obsIt = lower_bound(observedData.begin(), observedData.end(), t_start,
                                     [](const Observation &obs, double t)
                                     { return obs.time < t; });

            if (obsIt != observedData.end() && fabs(obsIt->time - t_start) < 5.0)
            {
                // Reset infectious count to observed data (Data Assimilation)
                // Use Moving Average to avoid outliers/noise at the boundary
                double observed_I = obsIt->movingAvg7Day;
                if (observed_I <= 0.1)
                    observed_I = obsIt->infectious; // Fallback if MA is not ready or 0

                // Ensure a minimum seed to allow Beta2 to work if cases are low
                if (observed_I < 1.0)
                    observed_I = 1.0;

                // FULL STATE RESET for independent yearly approximation
                currentState.I_h = observed_I;
                currentState.E_h = observed_I * 2.0; // Seed exposed population (2x infectious)
                currentState.Q_h = 0.0;
                currentState.R_h = 0.0;
                currentState.S_h = N_h_initial - (currentState.I_h + currentState.E_h);

                // Reset Rodents with higher background infection
                currentState.I_r = 100.0;
                currentState.E_r = 50.0;
                currentState.S_r = N_r_initial - currentState.I_r - currentState.E_r;
            }
        }

        cout << "Optimizing Period " << i << " (" << t_start << " to " << t_end << ")" << endl;
        cout << "  Initial I_h (from Data/MA): " << currentState.I_h << endl;

        double best_b = 0.0;
        double min_sse_period = numeric_limits<double>::max();
        State bestEndState = currentState;

        // Iterative Grid Search for high precision (Coarse-to-Fine)
        double current_step = 0.1;
        double search_start = 0.0;
        double search_end = 5.0; // Increased upper bound to allow higher transmission rates

        // Passes: 0.1 down to 1e-7
        for (int pass = 0; pass < 8; ++pass)
        {
            for (double b = search_start; b <= search_end + (current_step / 100.0); b += current_step)
            {
                if (b < 0.0)
                    continue;

                params.beta_schedule = optimized_schedule;
                params.beta_schedule[t_start] = b;

                // Simulate ONLY this period
                auto predictions = simulatePeriod(params, currentState, t_end, 0.1);

                double sse = 0.0;
                int count = 0;
                for (const auto &obs : observedData)
                {
                    if (obs.time >= t_start && obs.time < t_end)
                    {
                        auto it = lower_bound(predictions.begin(), predictions.end(), obs.time,
                                              [](const State &s, double t)
                                              { return s.t < t; });
                        if (it != predictions.end())
                        {
                            double total_infectious_h = it->I_h + it->Q_h;
                            double error = obs.infectious - total_infectious_h;
                            sse += error * error;
                            count++;
                        }
                    }
                }

                if (count > 0 && sse < min_sse_period)
                {
                    min_sse_period = sse;
                    best_b = b;
                    bestEndState = predictions.back();
                }
            }
            // Refine search range around the best found so far
            search_start = max(0.0, best_b - current_step);
            search_end = best_b + current_step;
            current_step /= 10.0;
        }
        if (min_sse_period == numeric_limits<double>::max())
        {
            best_b = (i > 0) ? optimized_schedule[boundaries[i - 1]] : 0.05;
            cout << "Period " << t_start << " to " << t_end << ": No data. Keeping Beta: " << best_b << endl;
            // Run simulation with default beta to get end state
            params.beta_schedule[t_start] = best_b;
            auto res = simulatePeriod(params, currentState, t_end, 0.1);
            bestEndState = res.back();
        }
        else
        {
            cout << "Period " << t_start << " to " << t_end << ": Best Beta: " << fixed << setprecision(7) << best_b << " (SSE: " << min_sse_period << ")" << endl;
        }

        optimized_schedule[t_start] = best_b;

        // Run the best simulation for this period again to append to history
        params.beta_schedule = optimized_schedule;
        auto periodResults = simulatePeriod(params, currentState, t_end, 0.1);

        // Append to full history (avoiding duplicates at boundaries)
        if (fullHistory.empty())
        {
            fullHistory = periodResults;
        }
        else
        {
            fullHistory.insert(fullHistory.end(), periodResults.begin() + 1, periodResults.end());
        }

        // Update currentState for next iteration
        currentState = bestEndState;
    }

    cout << "Optimization Complete." << endl;
    params.beta_schedule = optimized_schedule;

    // 5. Process Daily Results
    vector<DailyPrediction> dailyResults;
    auto rawResults = fullHistory; // Use the stitched history

    double final_t_max = rawResults.back().t;

    // Resample to daily intervals
    for (int day = 0; day <= (int)final_t_max; ++day)
    {
        // Find state at t = day
        auto it = lower_bound(rawResults.begin(), rawResults.end(), (double)day,
                              [](const State &s, double t)
                              { return s.t < t; });

        if (it != rawResults.end())
        {
            DailyPrediction dp;
            dp.time = it->t;
            dp.infectious = it->I_h + it->Q_h; // Total Active Human Cases
            dp.dateStr = formatDate(addDays(startDate, dp.time));
            dailyResults.push_back(dp);
        }
    }

    // Save final prediction
    ofstream outFile("monkeypox_fitted_prediction.csv");
    outFile << "Date,Time_Days,Infectious_Predicted\n";
    for (const auto &p : dailyResults)
    {
        outFile << p.dateStr << "," << p.time << "," << p.infectious << "\n";
    }
    outFile.close();
    cout << "Fitted model saved to 'monkeypox_fitted_prediction.csv'" << endl;

    // 6. Generate Visualization
    generateHTMLReport(observedData, dailyResults);

    return 0;
}
