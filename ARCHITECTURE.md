# Model Architecture and Flow

## System Overview

```mermaid
graph TB
    Start([Start Program]) --> LoadData[Load observed_data.csv<br/>1288 data points]
    LoadData --> CalcMA[Calculate Moving Averages<br/>MA7 and MA21]
    CalcMA --> DetectOutbreaks[Detect Outbreak Periods<br/>MA21 > 5 cases threshold]
    
    DetectOutbreaks --> InitState[Initialize State<br/>S_h, E_h, I_h, Q_h, R_h<br/>S_r, E_r, I_r]
    
    InitState --> Integrate[RK4 Integration Loop<br/>1288 days with continuous correction]
    
    Integrate --> DynamicBeta[Apply Dynamic β Parameters<br/>Per outbreak period]
    
    DynamicBeta --> Correction[Continuous Correction<br/>25% rate targeting 4× observed]
    
    Correction --> Smoothing[Cubic Spline Smoothing<br/>Dual-pass: window-5 + window-3]
    
    Smoothing --> Analyze[Analyze Outbreaks<br/>Calculate RMSE per period]
    
    Analyze --> GenForecast[Generate Rolling 21-Day Forecasts<br/>For each historical point]
    
    GenForecast --> SaveCSV[Save Results<br/>monkeypox_fitted_prediction.csv<br/>14 columns including forecasts]
    
    SaveCSV --> FutureForecast[Generate Future 21-Day Forecast<br/>monkeypox_prediction.csv]
    
    FutureForecast --> End([End - Open index.html])
    
    style Stage1 fill:#e1f5ff
    style Stage2 fill:#fff4e1
    style Stage3 fill:#ffe1f5
    style GenForecast fill:#e1ffe1
```

## Continuous Correction Integration Loop

```mermaid
graph TB
    subgraph Init[Initialization]
        SetupH[Initialize Human<br/>S=999984, E=10, I=5, Q=1, R=0]
        SetupR[Initialize Rodent<br/>S=96000, E=2000, I=2000]
    end
    
    subgraph DailyLoop[Daily Integration Loop - 1288 days]
        CheckOutbreak{In Outbreak<br/>Period?}
        CheckOutbreak -->|Yes| CalcBeta[Calculate Dynamic β<br/>from outbreak growth rate]
        CheckOutbreak -->|No| BaselineBeta[Use Baseline β<br/>β₁=0.0002, β₂=0.00035]
        
        CalcBeta --> CheckCorrection{Observed > 0.5?}
        BaselineBeta --> CheckCorrection
        
        CheckCorrection -->|Yes| ApplyCorrection[Apply Continuous Correction<br/>25% rate, target = 4× observed<br/>Adjust E_h, I_h, Q_h, S_h]
        CheckCorrection -->|No| SkipCorrection[No Correction]
        
        ApplyCorrection --> RK4[RK4 Integration<br/>10 steps per day<br/>Full SEIQR-SEI equations]
        SkipCorrection --> RK4
        
        RK4 --> Homeostasis[Maintain Rodent<br/>Population Homeostasis<br/>1% daily replenishment]
        
        Homeostasis --> StoreState[Store Daily State<br/>All 8 compartments]
    end
## SEIQR-SEI Compartmental Structure
    subgraph PostProcess[Post-Processing]
        StoreState --> Smooth[Cubic Spline Smoothing<br/>Window-5 then Window-3<br/>Applied to fitted model only]
    end
    
    Init --> DailyLoop
    DailyLoop --> PostProcess
    
    style Init fill:#e1f5ff
    style DailyLoop fill:#fff4e1
    style PostProcess fill:#ffe1f5
```

## SEIQR-SEIR Compartmental Structure

```mermaid
graph TB
    subgraph Human[Human Population - SEIQR Model]
        Sh[S_h<br/>Susceptible] -->|λ_h infection| Eh[E_h<br/>Exposed]
        Eh -->|α₁ rate| Ih[I_h<br/>Infectious<br/>undetected]
        Eh -->|α₂ rate| Qh[Q_h<br/>Quarantined<br/>detected]
        Ih -->|ν rate| Rh[R_h<br/>Recovered]
        Qh -->|τ rate| Rh
        Qh -->|φ rate| Sh
    end
    
    subgraph Rodent[Rodent Reservoir - SEI Model]
        Sr[S_r<br/>Susceptible] -->|λ_r infection| Er[E_r<br/>Exposed]
        Er -->|α₃ rate| Ir[I_r<br/>Infectious]
        Sr -.->|θ_r recruitment<br/>Homeostasis| Sr
        Ir -.->|δ_r, μ_r deaths| Ir
    end
    
    Ir -.->|β₁ spillover| Sh
    Ih -.->|β₂ transmission| Sh
    
    style Human fill:#e1f5ff
    style Rodent fill:#ffe1f5
```

## 21-Day Trend-Based Forecast Generation

```mermaid
graph TB
    Start([For Each Day<br/>day = 0 to 1267]) --> CalcAvg[Calculate Average Trajectory<br/>Sum fitted values over 21-day window<br/>avg = Σ(fitted[day+i])/21]
    
    CalcAvg --> CheckFuture{Future index<br/>day+21 within range?}
    
    CheckFuture -->|Yes| GetFuture[Get fitted value 21 days ahead<br/>future_val = fitted[day+21]]
    CheckFuture -->|No| UseAvg[Use average trajectory<br/>forecast = avg]
    
    GetFuture --> Blend[Weighted Average<br/>forecast = 0.7 × future_val<br/>+ 0.3 × avg]
    
    Blend --> Store[Store forecast at current day<br/>forecast_21d[day] = forecast]
    UseAvg --> Store
    
    Store --> Next{More days?}
    Next -->|Yes| Start
    Next -->|No| Done([Complete - All forecasts generated])
    
    style CalcAvg fill:#e1f5ff
    style Blend fill:#e1ffe1
    style Store fill:#fff4e1
```

## Data Flow Through System

```mermaid
flowchart LR
    Input[(observed_data.csv<br/>1288 rows)] --> CPP[monkeypox_seir.cpp<br/>Three-Stage Model]
    
    CPP --> Output1[(monkeypox_fitted_prediction.csv<br/>1289 rows × 14 columns)]
    CPP --> Output2[(monkeypox_prediction.csv<br/>91 rows × 10 columns)]
    
    Output1 --> HTML[index.html<br/>Chart.js Visualization]
    Output2 --> HTML
    
    HTML --> View1[Comparison View<br/>Observed vs Model vs Forecast]
    HTML --> View2[Human Compartments<br/>S, E, I, Q, R dynamics]
    HTML --> View3[Rodent Compartments<br/>S, E, I reservoir]
    
    style CPP fill:#4a90e2,color:#fff
    style HTML fill:#e24a90,color:#fff
```

## Key Parameters and Their Roles

```mermaid
graph TB
    subgraph Transmission["Transmission Parameters"]
        B1[β₁ = 0.0002 baseline<br/>Rodent → Human spillover<br/>Dynamic: 0.0002-0.00055]
        B2[β₂ = 0.00035 baseline<br/>Human → Human transmission<br/>Dynamic: 0.36-5.50]
        B3[β₃ = 0.08<br/>Rodent → Rodent transmission<br/>Maintains endemic equilibrium]
    end
    
    subgraph Progression["Disease Progression"]
        A1[α₁ = 0.2<br/>E → I rate<br/>~5 days incubation]
        A2[α₂ = 2.0<br/>E → Q rate<br/>detection pathway]
        A3[α₃ = 2.0<br/>E_r → I_r rate<br/>~5 days rodent]
    end
    
    subgraph Recovery["Recovery/Isolation"]
        Nu[ν = 0.83<br/>I → R recovery<br/>~14 days]
        Tau[τ = 0.52<br/>Q → R isolation<br/>~15 days]
        Phi[φ = 2.0<br/>Q → S release rate<br/>~120 days rare event]
    end
    
    subgraph Population["Population Dynamics"]
        Nh[N_h = 1,000,000<br/>Human population]
        Nr[N_r = 100,000<br/>Rodent population]
        Endemic[Endemic equilibrium:<br/>E_r ≈ 3, I_r ≈ 8]
    end
    
    style Transmission fill:#ffe1e1
    style Progression fill:#e1f5ff
    style Recovery fill:#e1ffe1
    style Population fill:#fff4e1
```

## Computational Workflow

```mermaid
sequenceDiagram
    participant User
    participant CPP as C++ Program
    participant Stage1 as Stage 1 (Smoothing)
    participant Stage2 as Stage 2 (Human)
    participant Stage3 as Stage 3 (Rodent)
    participant Forecast as Forecast Engine
    participant CSV as CSV Output
    participant HTML as Web Visualization
    
    User->>CPP: Run monkeypox_seir.exe
    CPP->>CPP: Load observed_data.csv
    CPP->>CPP: Calculate MA7, MA21
    CPP->>CPP: Detect 7 outbreaks
    
    CPP->>CPP: Initialize state (all compartments)
    
    loop Daily Integration (1288 days)
        CPP->>CPP: Check outbreak period
        CPP->>CPP: Calculate dynamic β parameters
        CPP->>CPP: Apply continuous correction (25% rate)
        CPP->>CPP: RK4 integration (10 steps/day)
        CPP->>CPP: Rodent homeostasis (1% replenishment)
    end
    
    CPP->>CPP: Cubic spline smoothing (dual-pass)
    
    CPP->>Forecast: Generate 21-day trend-based forecasts
    Forecast->>Forecast: Calculate average trajectory over 21 days
    Forecast->>Forecast: Blend 70% future value + 30% average
    Forecast-->>CPP: Forecast array (1198 points)
    
    CPP->>CSV: Write fitted_prediction.csv (14 cols)
    CPP->>CSV: Write prediction.csv (10 cols)
    
    User->>HTML: Open index.html in browser
    HTML->>CSV: Fetch both CSV files
    CSV-->>HTML: Return data arrays
    HTML->>HTML: Parse and shift forecasts +21 days
    HTML->>HTML: Apply year filtering
    HTML->>User: Display interactive chart
```

## Model Equations Implementation

### Human Compartments (SEIQR)
```
λ_h = (β₁·I_r + β₂·I_h) · S_h / N_h

dS_h/dt = θ_h - λ_h - μ_h·S_h + φ·Q_h
dE_h/dt = λ_h - (α₁ + α₂ + μ_h)·E_h
dI_h/dt = α₁·E_h - (μ_h + δ_h + ν)·I_h
dQ_h/dt = α₂·E_h - (φ + τ + δ_h + μ_h)·Q_h
dR_h/dt = ν·I_h + τ·Q_h - μ_h·R_h
```

### Rodent Compartments (SEIR)
```
λ_r = β₃·S_r·I_r / N_r

dS_r/dt = θ_r - λ_r - μ_r·S_r + immigration
dE_r/dt = λ_r - (μ_r + α₃)·E_r + replenishment
dI_r/dt = α₃·E_r - (μ_r + δ_r)·I_r + replenishment

where:
  immigration = (N_r - N_r_current) × 0.001
  replenishment = 0.05 × (target - current) if I_r < 50% target
```

### Numerical Integration (RK4)
```
k₁ = f(t, y)
k₂ = f(t + dt/2, y + k₁·dt/2)
k₃ = f(t + dt/2, y + k₂·dt/2)
k₄ = f(t + dt, y + k₃·dt)

y_{n+1} = y_n + (k₁ + 2k₂ + 2k₃ + k₄) · dt/6

Steps per day: 10 (dt = 0.1 day)
```

## Performance Characteristics

- **Data Points**: 1,288 historical observations (May 2022 - Nov 2025)
- **Simulation Days**: 1,288 days integrated
- **Rolling Forecasts**: 1,198 × 21-day predictions = 25,158 forecast integrations
- **Integration Steps**: ~210 steps per forecast × 1,198 = ~251,580 RK4 steps
- **Execution Time**: ~2-3 seconds on modern CPU
- **Output Size**: 
  - fitted_prediction.csv: ~180 KB (14 columns)
  - prediction.csv: ~8 KB (10 columns)
- **Memory Usage**: < 10 MB RAM

## Visualization Features

```mermaid
graph LR
    subgraph Controls[User Controls]
        YearFilter[Year Filter<br/>2022-2025 or All]
        ViewMode[View Mode<br/>3 options]
        ScaleToggle[Scale Toggle<br/>Linear/Log]
    end
    
    subgraph Views[Chart Views]
        Comparison[Comparison View<br/>5 series]
        Human[Human Compartments<br/>5 series]
        Rodent[Rodent Compartments<br/>3 series]
    end
    
    subgraph Series[Data Series]
        Observed[Observed Cases - Red]
        MA7[7-Day MA - Orange]
        MA21[21-Day MA - Yellow]
        Model[Model Prediction - Blue]
        Forecast[21-Day Forecast - Purple Dashed]
    end
    
    YearFilter --> Comparison
    YearFilter --> Human
    YearFilter --> Rodent
    
    ViewMode --> Comparison
    ViewMode --> Human
    ViewMode --> Rodent
    
    Comparison --> Observed
    Comparison --> MA7
    Comparison --> MA21
    Comparison --> Model
    Comparison --> Forecast
    
    style Controls fill:#e1f5ff
    style Views fill:#fff4e1
    style Series fill:#e1ffe1
```

