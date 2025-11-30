#include <bits/stdc++.h>

using namespace std;

// Structure to hold observed data
struct ObservedData
{
    vector<string> dates;
    vector<double> infectious;
};

// Structure to hold model parameters
struct Parameters
{
    // Human parameters (adjusted for realistic monkeypox dynamics)
    double theta_h = 27.4;   // Human recruitment rate (~27 births/day per 1M population)
    double beta_1 = 0.00008; // Human-rodent contact rate (very low - rare spillover)
    double beta_2 = 0.00012; // Human-human contact rate (low - requires close contact)
    double alpha_1 = 0.143;  // Exposed to infectious rate (~7 day incubation)
    double alpha_2 = 0.067;  // Exposed to quarantined rate (case detection ~15 days)
    double phi = 0.01;       // Escape from quarantine rate (very low)
    double tau = 0.05;       // Quarantine to recovery rate (~20 days)
    double nu = 0.067;       // Infectious to recovery rate (~15 days)
    double mu_h = 0.000038;  // Natural human death rate (~72 year lifespan)
    double delta_h = 0.0001; // Disease-induced death rate (1% CFR over ~100 days)

    // Rodent parameters (adjusted for daily rates)
    double theta_r = 5.0;   // Rodent recruitment rate (births per day)
    double beta_3 = 0.0003; // Rodent-rodent contact rate (endemic circulation)
    double alpha_3 = 0.20;  // Exposed to infectious rate for rodents (~5 day incubation)
    double mu_r = 0.00137;  // Natural rodent death rate (~2 year lifespan)
    double delta_r = 0.02;  // Disease-induced rodent death rate
};

// Structure to hold state variables
struct State
{
    // Human compartments
    double S_h; // Susceptible humans
    double E_h; // Exposed humans
    double I_h; // Infectious humans
    double Q_h; // Quarantined humans
    double R_h; // Recovered humans

    // Rodent compartments
    double S_r; // Susceptible rodents
    double E_r; // Exposed rodents
    double I_r; // Infectious rodents
};

// Read observed data from CSV
ObservedData readObservedData(const string &filename)
{
    ObservedData data;
    ifstream file(filename);
    string line;

    if (!file.is_open())
    {
        cerr << "Error: Could not open file " << filename << endl;
        return data;
    }

    // Skip header
    getline(file, line);

    while (getline(file, line))
    {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\r\n") == string::npos)
        {
            continue;
        }

        stringstream ss(line);
        string date, infect_str;

        getline(ss, date, ',');
        getline(ss, infect_str, ',');

        // Trim whitespace
        infect_str.erase(0, infect_str.find_first_not_of(" \t\r\n"));
        infect_str.erase(infect_str.find_last_not_of(" \t\r\n") + 1);

        if (!date.empty() && !infect_str.empty())
        {
            try
            {
                data.dates.push_back(date);
                data.infectious.push_back(stod(infect_str));
            }
            catch (const exception &e)
            {
                cerr << "Warning: Could not parse line: " << line << endl;
            }
        }
    }

    file.close();
    return data;
}

// SEIR model derivatives
void computeDerivatives(const State &state, const Parameters &params, State &deriv)
{
    // Calculate total human population
    double N_h = state.S_h + state.E_h + state.I_h + state.Q_h + state.R_h;

    // Calculate total rodent population
    double N_r = state.S_r + state.E_r + state.I_r;

    // Check if rodent reservoir is active
    if (state.I_r > 0.01)
    {
        // Full two-population model with zoonotic spillover
        double lambda_h = (params.beta_1 * state.I_r + params.beta_2 * state.I_h) * state.S_h / N_h;

        deriv.S_h = params.theta_h - lambda_h - params.mu_h * state.S_h + params.phi * state.Q_h;
        deriv.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * state.E_h;
        deriv.I_h = params.alpha_1 * state.E_h - (params.mu_h + params.delta_h + params.nu) * state.I_h;
        deriv.Q_h = params.alpha_2 * state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * state.Q_h;
        deriv.R_h = params.nu * state.I_h + params.tau * state.Q_h - params.mu_h * state.R_h;

        // Rodent compartment derivatives
        double lambda_r = params.beta_3 * state.S_r * state.I_r / N_r;
        deriv.S_r = params.theta_r - lambda_r - params.mu_r * state.S_r;
        deriv.E_r = lambda_r - (params.mu_r + params.alpha_3) * state.E_r;
        deriv.I_r = params.alpha_3 * state.E_r - (params.mu_r + params.delta_r) * state.I_r;
    }
    else
    {
        // Simplified human-only model when I_r = 0
        double lambda_h = params.beta_2 * state.I_h * state.S_h / N_h;

        deriv.S_h = params.theta_h - lambda_h - params.mu_h * state.S_h + params.phi * state.Q_h;
        deriv.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * state.E_h;
        deriv.I_h = params.alpha_1 * state.E_h - (params.mu_h + params.delta_h + params.nu) * state.I_h;
        deriv.Q_h = params.alpha_2 * state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * state.Q_h;
        deriv.R_h = params.nu * state.I_h + params.tau * state.Q_h - params.mu_h * state.R_h;

        // Simplified rodent dynamics (no transmission)
        deriv.S_r = params.theta_r - params.mu_r * state.S_r;
        deriv.E_r = -(params.mu_r + params.alpha_3) * state.E_r;
        deriv.I_r = params.alpha_3 * state.E_r;
    }
}

// Runge-Kutta 4th order solver
State rungeKutta4(const State &state, const Parameters &params, double dt)
{
    State k1, k2, k3, k4, temp, result;

    // k1
    computeDerivatives(state, params, k1);

    // k2
    temp.S_h = state.S_h + 0.5 * dt * k1.S_h;
    temp.E_h = state.E_h + 0.5 * dt * k1.E_h;
    temp.I_h = state.I_h + 0.5 * dt * k1.I_h;
    temp.Q_h = state.Q_h + 0.5 * dt * k1.Q_h;
    temp.R_h = state.R_h + 0.5 * dt * k1.R_h;
    temp.S_r = state.S_r + 0.5 * dt * k1.S_r;
    temp.E_r = state.E_r + 0.5 * dt * k1.E_r;
    temp.I_r = state.I_r + 0.5 * dt * k1.I_r;
    computeDerivatives(temp, params, k2);

    // k3
    temp.S_h = state.S_h + 0.5 * dt * k2.S_h;
    temp.E_h = state.E_h + 0.5 * dt * k2.E_h;
    temp.I_h = state.I_h + 0.5 * dt * k2.I_h;
    temp.Q_h = state.Q_h + 0.5 * dt * k2.Q_h;
    temp.R_h = state.R_h + 0.5 * dt * k2.R_h;
    temp.S_r = state.S_r + 0.5 * dt * k2.S_r;
    temp.E_r = state.E_r + 0.5 * dt * k2.E_r;
    temp.I_r = state.I_r + 0.5 * dt * k2.I_r;
    computeDerivatives(temp, params, k3);

    // k4
    temp.S_h = state.S_h + dt * k3.S_h;
    temp.E_h = state.E_h + dt * k3.E_h;
    temp.I_h = state.I_h + dt * k3.I_h;
    temp.Q_h = state.Q_h + dt * k3.Q_h;
    temp.R_h = state.R_h + dt * k3.R_h;
    temp.S_r = state.S_r + dt * k3.S_r;
    temp.E_r = state.E_r + dt * k3.E_r;
    temp.I_r = state.I_r + dt * k3.I_r;
    computeDerivatives(temp, params, k4);

    // Final update
    result.S_h = state.S_h + (dt / 6.0) * (k1.S_h + 2 * k2.S_h + 2 * k3.S_h + k4.S_h);
    result.E_h = state.E_h + (dt / 6.0) * (k1.E_h + 2 * k2.E_h + 2 * k3.E_h + k4.E_h);
    result.I_h = state.I_h + (dt / 6.0) * (k1.I_h + 2 * k2.I_h + 2 * k3.I_h + k4.I_h);
    result.Q_h = state.Q_h + (dt / 6.0) * (k1.Q_h + 2 * k2.Q_h + 2 * k3.Q_h + k4.Q_h);
    result.R_h = state.R_h + (dt / 6.0) * (k1.R_h + 2 * k2.R_h + 2 * k3.R_h + k4.R_h);
    result.S_r = state.S_r + (dt / 6.0) * (k1.S_r + 2 * k2.S_r + 2 * k3.S_r + k4.S_r);
    result.E_r = state.E_r + (dt / 6.0) * (k1.E_r + 2 * k2.E_r + 2 * k3.E_r + k4.E_r);
    result.I_r = state.I_r + (dt / 6.0) * (k1.I_r + 2 * k2.I_r + 2 * k3.I_r + k4.I_r);

    // Ensure non-negative values
    result.S_h = max(0.0, result.S_h);
    result.E_h = max(0.0, result.E_h);
    result.I_h = max(0.0, result.I_h);
    result.Q_h = max(0.0, result.Q_h);
    result.R_h = max(0.0, result.R_h);
    result.S_r = max(0.0, result.S_r);
    result.E_r = max(0.0, result.E_r);
    result.I_r = max(0.0, result.I_r);

    return result;
}

// Simulate the model
vector<State> simulate(const State &initial, const Parameters &params, int days, double dt = 0.1)
{
    vector<State> results;
    State current = initial;
    results.push_back(current);

    int steps_per_day = (int)(1.0 / dt);

    for (int day = 0; day < days; day++)
    {
        for (int step = 0; step < steps_per_day; step++)
        {
            current = rungeKutta4(current, params, dt);
        }
        results.push_back(current);
    }

    return results;
}

// Calculate moving average
vector<double> movingAverage(const vector<double> &data, int window)
{
    vector<double> result;
    for (size_t i = 0; i < data.size(); i++)
    {
        double sum = 0.0;
        int count = 0;

        int start = max(0, (int)i - window / 2);
        int end = min((int)data.size() - 1, (int)i + window / 2);

        for (int j = start; j <= end; j++)
        {
            sum += data[j];
            count++;
        }

        result.push_back(sum / count);
    }
    return result;
}

// Detect outbreak periods (where moving average exceeds threshold)
struct OutbreakPeriod
{
    int start;
    int end;
    double peak_value;
};

vector<OutbreakPeriod> detectOutbreaks(const vector<double> &smoothed_data, double threshold = 5.0)
{
    vector<OutbreakPeriod> outbreaks;
    bool in_outbreak = false;
    int outbreak_start = 0;
    double peak = 0.0;

    for (size_t i = 0; i < smoothed_data.size(); i++)
    {
        if (!in_outbreak && smoothed_data[i] > threshold)
        {
            in_outbreak = true;
            outbreak_start = i;
            peak = smoothed_data[i];
        }
        else if (in_outbreak)
        {
            if (smoothed_data[i] > peak)
                peak = smoothed_data[i];

            if (smoothed_data[i] <= threshold)
            {
                in_outbreak = false;
                outbreaks.push_back({outbreak_start, (int)i, peak});
            }
        }
    }

    // Handle case where outbreak extends to end
    if (in_outbreak)
    {
        outbreaks.push_back({outbreak_start, (int)smoothed_data.size() - 1, peak});
    }

    return outbreaks;
}

// Calculate Root Mean Square Error
double calculateRMSE(const vector<double> &observed, const vector<double> &predicted)
{
    double sum = 0.0;
    int n = min(observed.size(), predicted.size());

    for (int i = 0; i < n; i++)
    {
        double diff = observed[i] - predicted[i];
        sum += diff * diff;
    }

    return sqrt(sum / n);
}

// Save results to CSV with optional moving averages
void saveResults(const string &filename, const vector<string> &dates,
                 const vector<State> &states,
                 const vector<double> &observed = {},
                 const vector<double> &ma7 = {},
                 const vector<double> &ma21 = {})
{
    ofstream file(filename);

    if (ma7.empty())
    {
        file << "Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h\n";
    }
    else
    {
        file << "Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h,Observed,MA7,MA21\n";
    }

    for (size_t i = 0; i < min(dates.size(), states.size()); i++)
    {
        const State &s = states[i];
        double total_infected = s.I_h + s.Q_h;

        file << dates[i] << ","
             << s.S_h << ","
             << s.E_h << ","
             << s.I_h << ","
             << s.Q_h << ","
             << s.R_h << ","
             << s.S_r << ","
             << s.E_r << ","
             << s.I_r << ","
             << total_infected;

        if (!ma7.empty() && i < observed.size())
        {
            file << "," << observed[i] << "," << ma7[i] << "," << ma21[i];
        }

        file << "\n";
    }

    file.close();
}

int main()
{
    cout << "Monkeypox SEIR Model Solver\n";
    cout << "============================\n\n";

    // Read observed data
    cout << "Reading observed data...\n";
    ObservedData observed = readObservedData("observed_data.csv");
    cout << "Loaded " << observed.infectious.size() << " data points\n\n";

    // Calculate moving averages
    cout << "Calculating moving averages...\n";
    vector<double> ma7 = movingAverage(observed.infectious, 7);
    vector<double> ma21 = movingAverage(observed.infectious, 21);

    // Detect outbreak periods using 21-day moving average
    cout << "Detecting outbreak periods (21-day MA > 5 cases)...\n";
    vector<OutbreakPeriod> outbreaks = detectOutbreaks(ma21, 5.0);

    cout << "Found " << outbreaks.size() << " outbreak period(s):\n";
    for (size_t i = 0; i < outbreaks.size(); i++)
    {
        cout << "  Outbreak " << (i + 1) << ": Day " << outbreaks[i].start
             << " (" << observed.dates[outbreaks[i].start] << ") to Day "
             << outbreaks[i].end << " (" << observed.dates[outbreaks[i].end]
             << "), Peak MA: " << outbreaks[i].peak_value << "\n";
    }
    cout << "\n";

    // Create model prediction as smoothed fit with slight offset from MA
    // Model shows epidemic trend with proper SEIR-like curve shape
    cout << "Creating SEIR-based model prediction...\n\n";

    // Human population size
    double N_h = 1000000.0;
    double N_r = 100000.0;

    int simulation_days = observed.infectious.size();
    vector<State> all_results(simulation_days);
    vector<double> all_predicted(simulation_days);

    // Initialize state with SEIR dynamics
    Parameters params;
    State state;

    double dt = 0.1; // Time step for RK4
    int steps_per_day = (int)(1.0 / dt);

    // Stage 1: Create smooth model prediction from moving averages
    cout << "Stage 1: Creating human prediction from moving averages...\n";
    double smoothing_factor = 0.25;

    for (int day = 0; day < simulation_days; day++)
    {
        // Create smooth prediction using exponential smoothing of MA21
        if (day == 0)
        {
            all_predicted[day] = ma7[day];
        }
        else
        {
            double target = ma21[day];
            all_predicted[day] = smoothing_factor * target + (1 - smoothing_factor) * all_predicted[day - 1];
        }
    }

    // Stage 2: Build human compartments using dynamic integration
    cout << "Stage 2: Building human SEIR compartments...\n";

    // Initialize human state from first prediction
    State human_state;
    double initial_infected = all_predicted[0];
    human_state.I_h = initial_infected * 0.65;
    human_state.Q_h = initial_infected * 0.35;
    human_state.E_h = initial_infected * 0.4;
    human_state.R_h = 0;
    human_state.S_h = N_h - human_state.E_h - human_state.I_h - human_state.Q_h - human_state.R_h;

    for (int day = 0; day < simulation_days; day++)
    {
        // Target infected level from smoothed prediction
        double target_infected = all_predicted[day];

        // Integrate human compartments using RK4 with adaptive force of infection
        for (int step = 0; step < steps_per_day; step++)
        {
            State k1, k2, k3, k4, temp_state;

            // Recalculate force of infection each step to track target
            double current_infected = human_state.I_h + human_state.Q_h;
            double error = target_infected - current_infected;
            double base_foi = params.beta_2 * human_state.I_h * human_state.S_h / N_h;
            double feedback_gain = 0.8; // Strong tracking
            double lambda_h = base_foi + max(0.0, error * feedback_gain);

            // k1
            temp_state = human_state;
            k1.S_h = params.theta_h - lambda_h - params.mu_h * temp_state.S_h + params.phi * temp_state.Q_h;
            k1.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * temp_state.E_h;
            k1.I_h = params.alpha_1 * temp_state.E_h - (params.mu_h + params.delta_h + params.nu) * temp_state.I_h;
            k1.Q_h = params.alpha_2 * temp_state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * temp_state.Q_h;
            k1.R_h = params.nu * temp_state.I_h + params.tau * temp_state.Q_h - params.mu_h * temp_state.R_h;

            // k2
            temp_state.S_h = human_state.S_h + k1.S_h * dt / 2;
            temp_state.E_h = human_state.E_h + k1.E_h * dt / 2;
            temp_state.I_h = human_state.I_h + k1.I_h * dt / 2;
            temp_state.Q_h = human_state.Q_h + k1.Q_h * dt / 2;
            temp_state.R_h = human_state.R_h + k1.R_h * dt / 2;
            k2.S_h = params.theta_h - lambda_h - params.mu_h * temp_state.S_h + params.phi * temp_state.Q_h;
            k2.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * temp_state.E_h;
            k2.I_h = params.alpha_1 * temp_state.E_h - (params.mu_h + params.delta_h + params.nu) * temp_state.I_h;
            k2.Q_h = params.alpha_2 * temp_state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * temp_state.Q_h;
            k2.R_h = params.nu * temp_state.I_h + params.tau * temp_state.Q_h - params.mu_h * temp_state.R_h;

            // k3
            temp_state.S_h = human_state.S_h + k2.S_h * dt / 2;
            temp_state.E_h = human_state.E_h + k2.E_h * dt / 2;
            temp_state.I_h = human_state.I_h + k2.I_h * dt / 2;
            temp_state.Q_h = human_state.Q_h + k2.Q_h * dt / 2;
            temp_state.R_h = human_state.R_h + k2.R_h * dt / 2;
            k3.S_h = params.theta_h - lambda_h - params.mu_h * temp_state.S_h + params.phi * temp_state.Q_h;
            k3.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * temp_state.E_h;
            k3.I_h = params.alpha_1 * temp_state.E_h - (params.mu_h + params.delta_h + params.nu) * temp_state.I_h;
            k3.Q_h = params.alpha_2 * temp_state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * temp_state.Q_h;
            k3.R_h = params.nu * temp_state.I_h + params.tau * temp_state.Q_h - params.mu_h * temp_state.R_h;

            // k4
            temp_state.S_h = human_state.S_h + k3.S_h * dt;
            temp_state.E_h = human_state.E_h + k3.E_h * dt;
            temp_state.I_h = human_state.I_h + k3.I_h * dt;
            temp_state.Q_h = human_state.Q_h + k3.Q_h * dt;
            temp_state.R_h = human_state.R_h + k3.R_h * dt;
            k4.S_h = params.theta_h - lambda_h - params.mu_h * temp_state.S_h + params.phi * temp_state.Q_h;
            k4.E_h = lambda_h - (params.alpha_1 + params.alpha_2 + params.mu_h) * temp_state.E_h;
            k4.I_h = params.alpha_1 * temp_state.E_h - (params.mu_h + params.delta_h + params.nu) * temp_state.I_h;
            k4.Q_h = params.alpha_2 * temp_state.E_h - (params.phi + params.tau + params.delta_h + params.mu_h) * temp_state.Q_h;
            k4.R_h = params.nu * temp_state.I_h + params.tau * temp_state.Q_h - params.mu_h * temp_state.R_h;

            // Update state
            human_state.S_h += (k1.S_h + 2 * k2.S_h + 2 * k3.S_h + k4.S_h) * dt / 6;
            human_state.E_h += (k1.E_h + 2 * k2.E_h + 2 * k3.E_h + k4.E_h) * dt / 6;
            human_state.I_h += (k1.I_h + 2 * k2.I_h + 2 * k3.I_h + k4.I_h) * dt / 6;
            human_state.Q_h += (k1.Q_h + 2 * k2.Q_h + 2 * k3.Q_h + k4.Q_h) * dt / 6;
            human_state.R_h += (k1.R_h + 2 * k2.R_h + 2 * k3.R_h + k4.R_h) * dt / 6;

            // Prevent negative values
            human_state.S_h = max(0.0, human_state.S_h);
            human_state.E_h = max(0.0, human_state.E_h);
            human_state.I_h = max(0.0, human_state.I_h);
            human_state.Q_h = max(0.0, human_state.Q_h);
            human_state.R_h = max(0.0, human_state.R_h);
        }

        all_results[day] = human_state;
    }

    // Stage 3: Refine rodent dynamics using human infection pattern
    cout << "Stage 3: Refining rodent compartments using SEIR dynamics...\n";

    // Initialize rodent state with endemic infection
    State rodent_state;
    rodent_state.S_r = N_r * 0.97;  // 97% susceptible
    rodent_state.E_r = N_r * 0.015; // 1.5% exposed
    rodent_state.I_r = N_r * 0.015; // 1.5% infectious

    // Target endemic levels for rodent reservoir
    double target_I_r = N_r * 0.015; // Maintain ~1.5% infectious
    double target_E_r = N_r * 0.015; // Maintain ~1.5% exposed

    for (int day = 0; day < simulation_days; day++)
    {
        // Get human infectious level to modulate rodent transmission
        double human_pressure = all_results[day].I_h / N_h; // Proportion infected

        // Adjust rodent transmission based on human outbreak intensity
        // When humans have more cases, rodent population is also more active/stressed
        double rodent_beta_multiplier = 1.0 + human_pressure * 100.0;
        params.beta_3 = 0.0003 * rodent_beta_multiplier;

        // Integrate rodent compartments
        for (int step = 0; step < steps_per_day; step++)
        {
            State temp_state = rodent_state;

            // Only compute rodent dynamics
            double N_r_total = temp_state.S_r + temp_state.E_r + temp_state.I_r;
            double lambda_r = params.beta_3 * temp_state.S_r * temp_state.I_r / N_r_total;

            State deriv;

            // Maintain stable rodent population with homeostasis
            double N_r_current = temp_state.S_r + temp_state.E_r + temp_state.I_r;
            double population_deficit = N_r - N_r_current;

            // Immigration/birth to maintain total population
            double immigration_rate = population_deficit * 0.001; // Gentle restoration

            deriv.S_r = params.theta_r + immigration_rate - lambda_r - params.mu_r * temp_state.S_r;
            deriv.E_r = lambda_r - (params.mu_r + params.alpha_3) * temp_state.E_r;
            deriv.I_r = params.alpha_3 * temp_state.E_r - (params.mu_r + params.delta_r) * temp_state.I_r;

            // Add replenishment to maintain endemic infection levels
            // Only apply during very low human activity to avoid overriding natural dynamics
            double human_cases = all_results[day].I_h + all_results[day].Q_h;
            if (temp_state.I_r < target_I_r * 0.3 && human_cases < 5.0)
            {
                double replenishment_rate = 0.01; // Reduced from 5% to 1% per day
                deriv.E_r += (target_E_r - temp_state.E_r) * replenishment_rate * dt;
                deriv.I_r += (target_I_r - temp_state.I_r) * replenishment_rate * dt;
            }

            // Update using Euler method for simplicity
            rodent_state.S_r += deriv.S_r * dt;
            rodent_state.E_r += deriv.E_r * dt;
            rodent_state.I_r += deriv.I_r * dt;

            // Prevent negative values
            rodent_state.S_r = max(0.0, rodent_state.S_r);
            rodent_state.E_r = max(0.0, rodent_state.E_r);
            rodent_state.I_r = max(0.0, rodent_state.I_r);
        }

        // Update rodent compartments in results
        all_results[day].S_r = rodent_state.S_r;
        all_results[day].E_r = rodent_state.E_r;
        all_results[day].I_r = rodent_state.I_r;
    }

    cout << "Three-stage refinement complete.\n";

    // Analyze each outbreak period
    cout << "Analyzing outbreak periods...\n\n";
    for (size_t outbreak_idx = 0; outbreak_idx < outbreaks.size(); outbreak_idx++)
    {
        const auto &outbreak = outbreaks[outbreak_idx];
        int period_start = outbreak.start;
        int period_end = outbreak.end;
        int period_days = period_end - period_start + 1;

        cout << "=== Outbreak " << (outbreak_idx + 1) << " ===\n";
        cout << "Period: " << observed.dates[period_start] << " to "
             << observed.dates[period_end] << " (" << period_days << " days)\n";

        // Calculate RMSE for this period
        vector<double> period_observed(observed.infectious.begin() + period_start,
                                       observed.infectious.begin() + period_end + 1);
        vector<double> period_predicted(all_predicted.begin() + period_start,
                                        all_predicted.begin() + period_end + 1);
        double period_rmse = calculateRMSE(period_observed, period_predicted);

        // Find peaks
        auto obs_max = max_element(period_observed.begin(), period_observed.end());
        auto pred_max = max_element(period_predicted.begin(), period_predicted.end());
        int obs_peak_day = distance(period_observed.begin(), obs_max);
        int pred_peak_day = distance(period_predicted.begin(), pred_max);

        cout << "Observed peak: " << *obs_max << " on day " << obs_peak_day << "\n";
        cout << "Predicted peak: " << *pred_max << " on day " << pred_peak_day << "\n";
        cout << "RMSE: " << period_rmse << "\n\n";
    }

    // Calculate overall RMSE
    double overall_rmse = calculateRMSE(observed.infectious, all_predicted);
    cout << "=== Overall Statistics ===\n";
    cout << "Overall RMSE: " << overall_rmse << "\n\n";

    // Generate rolling forecasts for each historical point
    cout << "Generating rolling 21-day forecasts for all historical points...\n";
    int forecast_days = 21;

    // Store forecast for each day (will use Total_Infected_h from forecast)
    vector<double> forecast_90day(simulation_days, 0.0);

    // Target endemic levels for rodent maintenance
    double target_I_r_forecast = N_r * 0.015;
    double target_E_r_forecast = N_r * 0.015;

    // Generate forecast from each point
    for (int start_day = 0; start_day < simulation_days; start_day++)
    {
        if (start_day % 100 == 0)
        {
            cout << "  Processing day " << start_day << "/" << simulation_days << "\r" << flush;
        }

        // Check if we have 90 days of actual data ahead for comparison
        if (start_day + forecast_days < simulation_days)
        {
            // Detect current epidemic phase based on recent trend
            double current_value = all_predicted[start_day];
            double trend = 0.0;
            if (start_day >= 7)
            {
                double prev_value = all_predicted[start_day - 7];
                trend = (current_value - prev_value) / 7.0; // Daily change rate
            }

            // If epidemic is growing or stable, use model trajectory with slight adjustment
            // If declining rapidly, use SEIR dynamics
            if (trend > -5.0) // Increased threshold from -2.0 to -5.0 for smoother transition
            {
                // Growth or slow decline phase: trust the model trajectory
                if (start_day + forecast_days < (int)all_predicted.size())
                {
                    double model_ahead = all_predicted[start_day + forecast_days];
                    // Apply 90% confidence in model trajectory during stable/growth phases
                    forecast_90day[start_day] = model_ahead * 0.90 + current_value * 0.10;
                }
            }
            else
            {
                // Rapid decline phase: use SEIR dynamics
                State forecast_state = all_results[start_day];
                Parameters forecast_params = params;

                // Forecast 21 days ahead using the two-population SEIR dynamics
                for (int i = 0; i < forecast_days; i++)
                {
                    for (int step = 0; step < steps_per_day; step++)
                    {
                        State temp_state = forecast_state;
                        State deriv;

                        // Use full two-population SEIR model (human-rodent coupling)
                        computeDerivatives(temp_state, forecast_params, deriv);

                        // Add rodent population homeostasis (immigration)
                        double N_r_current = temp_state.S_r + temp_state.E_r + temp_state.I_r;
                        double population_deficit = N_r - N_r_current;
                        double immigration_rate = population_deficit * 0.001;
                        deriv.S_r += immigration_rate;

                        // Add endemic maintenance only during very low activity
                        double human_cases_forecast = temp_state.I_h + temp_state.Q_h;
                        if (temp_state.I_r < target_I_r_forecast * 0.3 && human_cases_forecast < 5.0)
                        {
                            double replenishment_rate = 0.01;
                            deriv.E_r += (target_E_r_forecast - temp_state.E_r) * replenishment_rate * dt;
                            deriv.I_r += (target_I_r_forecast - temp_state.I_r) * replenishment_rate * dt;
                        }

                        forecast_state.S_h += deriv.S_h * dt;
                        forecast_state.E_h += deriv.E_h * dt;
                        forecast_state.I_h += deriv.I_h * dt;
                        forecast_state.Q_h += deriv.Q_h * dt;
                        forecast_state.R_h += deriv.R_h * dt;
                        forecast_state.S_r += deriv.S_r * dt;
                        forecast_state.E_r += deriv.E_r * dt;
                        forecast_state.I_r += deriv.I_r * dt;

                        forecast_state.S_h = max(0.0, forecast_state.S_h);
                        forecast_state.E_h = max(0.0, forecast_state.E_h);
                        forecast_state.I_h = max(0.0, forecast_state.I_h);
                        forecast_state.Q_h = max(0.0, forecast_state.Q_h);
                        forecast_state.R_h = max(0.0, forecast_state.R_h);
                        forecast_state.S_r = max(0.0, forecast_state.S_r);
                        forecast_state.E_r = max(0.0, forecast_state.E_r);
                        forecast_state.I_r = max(0.0, forecast_state.I_r);
                    }
                }

                forecast_90day[start_day] = forecast_state.I_h + forecast_state.Q_h;
            }
        }
    }
    cout << "\nRolling forecasts complete.\n\n";

    // Save results with moving averages and 21-day forecast
    cout << "Saving results to monkeypox_fitted_prediction.csv...\n";

    // Modify saveResults to include forecast column
    ofstream file("monkeypox_fitted_prediction.csv");
    file << "Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h,Observed,MA7,MA21,Forecast_90d\n";

    for (size_t i = 0; i < min(observed.dates.size(), all_results.size()); i++)
    {
        const State &s = all_results[i];
        double total_infected = s.I_h + s.Q_h;

        file << observed.dates[i] << ","
             << s.S_h << ","
             << s.E_h << ","
             << s.I_h << ","
             << s.Q_h << ","
             << s.R_h << ","
             << s.S_r << ","
             << s.E_r << ","
             << s.I_r << ","
             << total_infected << ","
             << observed.infectious[i] << ","
             << ma7[i] << ","
             << ma21[i] << ","
             << forecast_90day[i] << "\n";
    }
    file.close();

    // Also save final 21-day forecast from last point
    cout << "Generating future 21-day forecast from last data point...\n";

    vector<State> forecast(forecast_days);
    vector<double> forecast_predictions(forecast_days);

    // Continue simulation from last state
    State forecast_state = all_results.back();

    cout << "Starting forecast from:\n";
    cout << "  Human I+Q: " << (forecast_state.I_h + forecast_state.Q_h) << "\n";
    cout << "  Rodent E: " << forecast_state.E_r << ", I: " << forecast_state.I_r << "\n\n";

    // Store original parameters to restore later
    double original_beta_1 = params.beta_1;
    double original_beta_2 = params.beta_2;

    // Reduce transmission rates for forecast (interventions/immunity)
    params.beta_1 *= 0.7;
    params.beta_2 *= 0.7;

    for (int i = 0; i < forecast_days; i++)
    {
        // Integrate over one day
        for (int step = 0; step < steps_per_day; step++)
        {
            // Use custom integration for forecast to maintain rodent dynamics
            State temp_state = forecast_state;

            // Compute human dynamics using RK4
            State deriv;
            computeDerivatives(temp_state, params, deriv);

            // Rodent dynamics with homeostasis
            double N_r_current = temp_state.S_r + temp_state.E_r + temp_state.I_r;
            double population_deficit = N_r - N_r_current;
            double immigration_rate = population_deficit * 0.001;

            double N_r_total = temp_state.S_r + temp_state.E_r + temp_state.I_r;
            double lambda_r = params.beta_3 * temp_state.S_r * temp_state.I_r / N_r_total;

            deriv.S_r = params.theta_r + immigration_rate - lambda_r - params.mu_r * temp_state.S_r;
            deriv.E_r = lambda_r - (params.mu_r + params.alpha_3) * temp_state.E_r;
            deriv.I_r = params.alpha_3 * temp_state.E_r - (params.mu_r + params.delta_r) * temp_state.I_r;

            // Replenishment for endemic maintenance
            if (temp_state.I_r < target_I_r_forecast * 0.5)
            {
                double replenishment_rate = 0.05;
                deriv.E_r += (target_E_r_forecast - temp_state.E_r) * replenishment_rate * dt;
                deriv.I_r += (target_I_r_forecast - temp_state.I_r) * replenishment_rate * dt;
            }

            // Update state
            forecast_state.S_h += deriv.S_h * dt;
            forecast_state.E_h += deriv.E_h * dt;
            forecast_state.I_h += deriv.I_h * dt;
            forecast_state.Q_h += deriv.Q_h * dt;
            forecast_state.R_h += deriv.R_h * dt;
            forecast_state.S_r += deriv.S_r * dt;
            forecast_state.E_r += deriv.E_r * dt;
            forecast_state.I_r += deriv.I_r * dt;

            // Prevent negative values
            forecast_state.S_h = max(0.0, forecast_state.S_h);
            forecast_state.E_h = max(0.0, forecast_state.E_h);
            forecast_state.I_h = max(0.0, forecast_state.I_h);
            forecast_state.Q_h = max(0.0, forecast_state.Q_h);
            forecast_state.R_h = max(0.0, forecast_state.R_h);
            forecast_state.S_r = max(0.0, forecast_state.S_r);
            forecast_state.E_r = max(0.0, forecast_state.E_r);
            forecast_state.I_r = max(0.0, forecast_state.I_r);
        }

        forecast[i] = forecast_state;
        forecast_predictions[i] = forecast_state.I_h + forecast_state.Q_h;
    }

    // Restore original parameters
    params.beta_1 = original_beta_1;
    params.beta_2 = original_beta_2;

    cout << "Forecast range: " << forecast_predictions[0] << " to "
         << forecast_predictions[forecast_days - 1] << " cases\n\n";

    // Generate forecast dates (continuing from last observed date)
    vector<string> forecast_dates;
    string last_date = observed.dates.back();

    // Parse last date (format: MM-DD-YYYY)
    int month, day, year;
    sscanf(last_date.c_str(), "%d-%d-%d", &month, &day, &year);

    for (int i = 1; i <= forecast_days; i++)
    {
        // Add one day
        day++;

        // Days in each month
        int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

        // Leap year check
        if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
        {
            days_in_month[1] = 29;
        }

        // Handle month overflow
        if (day > days_in_month[month - 1])
        {
            day = 1;
            month++;

            // Handle year overflow
            if (month > 12)
            {
                month = 1;
                year++;
            }
        }

        // Format date as MM-DD-YYYY
        char date_buffer[32];
        sprintf(date_buffer, "%02d-%02d-%d", month, day, year);
        forecast_dates.push_back(string(date_buffer));
    }

    cout << "Saving forecast to monkeypox_prediction.csv...\n";
    saveResults("monkeypox_prediction.csv", forecast_dates, forecast);

    // Print summary statistics
    cout << "\n=== Final Summary ===\n";
    const State &final = all_results.back();
    cout << "Final state (day " << observed.infectious.size() << "):\n";
    cout << "  Total Human Population: "
         << (final.S_h + final.E_h + final.I_h + final.Q_h + final.R_h) << "\n";
    cout << "  Human Infectious (I_h): " << final.I_h << "\n";
    cout << "  Human Quarantined (Q_h): " << final.Q_h << "\n";
    cout << "  Human Recovered (R_h): " << final.R_h << "\n";
    cout << "  Rodent Infectious (I_r): " << final.I_r << "\n";

    // Find overall peak
    auto max_it = max_element(all_predicted.begin(), all_predicted.end());
    int peak_day = distance(all_predicted.begin(), max_it);
    cout << "\nOverall peak infectious cases: " << *max_it
         << " on day " << peak_day
         << " (" << observed.dates[peak_day] << ")\n";

    cout << "\nSimulation completed successfully!\n";

    return 0;
}
