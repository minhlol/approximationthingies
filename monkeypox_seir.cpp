#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

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
    // Human parameters (from reference table with exact precision)
    double theta_h = 0.029;  // Human recruitment rate
    double beta_1 = 0.00025; // Human-rodent contact rate (2.5 × 10^-4)
    double beta_2 = 0.0006;  // Human-human contact rate (6 × 10^-4)
    double alpha_1 = 0.20;   // Proportion of infected humans from exposed humans
    double alpha_2 = 2.00;   // Proportion of suspected cases detected
    double phi = 2.00;       // Proportion not detected after diagnosis
    double tau = 0.52;       // Progression from isolated class to recovered class
    double nu = 0.83;        // Humans recovery rate
    double mu_h = 1.50;      // Natural (human) death rate
    double delta_h = 0.20;   // Disease (human) induced death rate

    // Rodent parameters (from reference table with exact precision)
    double theta_r = 0.20; // Rodents recruitment rate
    double beta_3 = 0.027; // Rodent-rodent contact rate
    double alpha_3 = 2.00; // Proportion of infected rodents from exposed rodents
    double mu_r = 0.020;   // Natural (rodents) death rate (2 × 10^-2)
    double delta_r = 0.50; // Disease (rodents) induced death rate
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

    // Human compartment derivatives
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

// Save results to CSV
void saveResults(const string &filename, const vector<string> &dates,
                 const vector<State> &states)
{
    ofstream file(filename);

    file << "Date,S_h,E_h,I_h,Q_h,R_h,S_r,E_r,I_r,Total_Infected_h\n";

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
             << total_infected << "\n";
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

    // Initialize parameters
    Parameters params;

    // Initialize state with reasonable values based on early 2022 monkeypox outbreak
    State initial;
    // Human population (normalized to 1, representing proportions)
    // Starting with very small initial infections consistent with outbreak start
    initial.S_h = 0.999999;  // Nearly all susceptible
    initial.E_h = 0.0000005; // Very few exposed
    initial.I_h = 0.0000004; // Initial infectious cases
    initial.Q_h = 0.0000001; // Initial quarantined
    initial.R_h = 0.0;       // No recovered initially

    // Rodent population (normalized to 1, representing proportions)
    initial.S_r = 0.9999;  // Most rodents susceptible
    initial.E_r = 0.00005; // Small exposed rodent population
    initial.I_r = 0.00005; // Small infected rodent population (reservoir)

    cout << "Initial Conditions:\n";
    cout << "  Human Susceptible: " << initial.S_h << "\n";
    cout << "  Human Exposed: " << initial.E_h << "\n";
    cout << "  Human Infectious: " << initial.I_h << "\n";
    cout << "  Rodent Susceptible: " << initial.S_r << "\n";
    cout << "  Rodent Exposed: " << initial.E_r << "\n";
    cout << "  Rodent Infectious: " << initial.I_r << "\n\n";

    // Simulate
    cout << "Running simulation...\n";
    int simulation_days = observed.infectious.size();
    vector<State> results = simulate(initial, params, simulation_days, 0.1);

    // Extract predicted infectious cases
    vector<double> predicted_infectious;
    for (const auto &state : results)
    {
        predicted_infectious.push_back(state.I_h + state.Q_h);
    }

    // Calculate RMSE
    double rmse = calculateRMSE(observed.infectious, predicted_infectious);
    cout << "RMSE: " << rmse << "\n\n";

    // Save results
    cout << "Saving results to monkeypox_fitted_prediction.csv...\n";
    saveResults("monkeypox_fitted_prediction.csv", observed.dates, results);

    // Extend prediction for 30 more days
    cout << "Generating 30-day forecast...\n";
    int forecast_days = 30;
    State last_state = results.back();
    vector<State> forecast = simulate(last_state, params, forecast_days, 0.1);

    // Generate forecast dates
    vector<string> forecast_dates;
    for (int i = 1; i <= forecast_days; i++)
    {
        forecast_dates.push_back("Day+" + to_string(i));
    }

    cout << "Saving forecast to monkeypox_prediction.csv...\n";
    saveResults("monkeypox_prediction.csv", forecast_dates, forecast);

    // Print summary statistics
    cout << "\n=== Simulation Summary ===\n";
    const State &final = results.back();
    cout << "Final state (day " << simulation_days << "):\n";
    cout << "  Total Human Population: "
         << (final.S_h + final.E_h + final.I_h + final.Q_h + final.R_h) << "\n";
    cout << "  Human Infectious (I_h): " << final.I_h << "\n";
    cout << "  Human Quarantined (Q_h): " << final.Q_h << "\n";
    cout << "  Human Recovered (R_h): " << final.R_h << "\n";
    cout << "  Rodent Infectious (I_r): " << final.I_r << "\n";

    // Find peak
    auto max_it = max_element(predicted_infectious.begin(), predicted_infectious.end());
    int peak_day = distance(predicted_infectious.begin(), max_it);
    cout << "\nPeak infectious cases: " << *max_it
         << " on day " << peak_day
         << " (" << observed.dates[peak_day] << ")\n";

    cout << "\nSimulation completed successfully!\n";

    return 0;
}
