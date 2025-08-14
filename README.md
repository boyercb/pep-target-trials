# pep-target-trials

Replication files for:

Boyer, C. & Lipsitch, M. (2025). "Emulating target trials of postexposure vaccines using observational data". American Journal of Epidemiology, 194(7), 2037-2046. https://doi.org/10.1093/aje/kwae350

## Abstract

Postexposure vaccination has the potential to prevent or modify the course of clinical disease among those exposed to a pathogen. However, due to logistical constraints, postexposure vaccine trials have been difficult to implement in practice. In place of trials, investigators have used observational data to estimate the effectiveness or optimal timing window for postexposure vaccines, but the relationship between these analyses and those that would be conducted in a trial is often unclear. Here, we define several possible target trials for postexposure vaccination and show how, under certain conditions, they can be emulated using observational data. We emphasize the importance of the incubation period and the timing of vaccination in trial design and emulation. As an example, we specify a protocol for postexposure vaccination against mpox and provide a step-by-step description of how to emulate it using data from a health care database or contact tracing program. We further illustrate some of the benefits of the target trial approach through simulation.

## Repository Structure

```
├── __master_run.R          # Main execution script for all simulations
├── 0_data/                 # Simulation results and datasets
│   ├── sim1.rds           # Simulation scenario 1 results (created)
│   ├── sim2.rds           # Simulation scenario 2 results (created)
│   ├── sim3.rds           # Simulation scenario 3 results (created)
│   ├── sim4.rds           # Simulation scenario 4 results (created)
│   └── sim5.rds           # Simulation scenario 5 results (created)
├── 1_code/                 # R source code
│   ├── 0_packages.R       # Required R packages
│   ├── 1_functions.R      # Core simulation and analysis functions
│   ├── 2_sims.R          # Main simulation scenarios
│   ├── 3_plots.R         # Figure generation code
│   └── 4_tables.R         # Table generation code
├── 2_tables/              # Generated LaTeX tables
│   ├── sim_results_hr.tex # Hazard ratio simulation results
│   └── sim_results_rr.tex # Relative risk simulation results
├── 3_figures/             # Generated figures and plots
│   ├── sim.pdf           # Main simulation results figure
│   ├── sim_hr.pdf        # Hazard ratio results
│   ├── sim_rr.pdf        # Relative risk results
│   ├── sim_overlap.pdf   # Overlap diagnostics
│   ├── sim_hetx.pdf      # Heterogeneous effects
│   └── dist.pdf          # Distribution plots
└── 4_manuscripts/         # Manuscript files
    ├── main.tex          # Main manuscript LaTeX source
    ├── supplement.tex    # Supplementary materials
    ├── pep.bib          # Bibliography file
    ├── setup.tex        # LaTeX setup and packages
    └── submissions/      # Journal submission versions
```

## Getting Started

### Prerequisites

This project requires R (version 4.0 or higher) with the following packages:
- `data.table` - Fast data manipulation
- `tidyverse` - Data science toolkit
- `survival` - Survival analysis
- `splines` - Spline functions
- `progressr` - Progress reporting
- `scales` - Scale functions for visualization
- `patchwork` - Composing plots
- `knitr` & `kableExtra` - Table formatting
- `boot` - Bootstrap methods
- `parallel` - Parallel computing
- `tictoc` - Timing functions

### Running the Analysis

1. **Setup**: Install required packages by running:
   ```r
   source("1_code/0_packages.R")
   ```

2. **Full Replication**: Execute all simulations and generate results:
   ```r
   source("__master_run.R")
   ```

3. **Individual Components**: Run specific parts:
   ```r
   # Load functions
   source("1_code/1_functions.R")
   
   # Run simulations
   source("1_code/2_sims.R")
   
   # Generate plots
   source("1_code/3_plots.R")

   # Generate tables
   source("1_code/4_tables.R")
   ```

### Simulation Parameters

The main simulations can be customized by modifying parameters in `__master_run.R`:
- `N_OBS`: Sample size per simulation (default: 1000)
- `N_SIMS`: Number of simulation replicates (default: 1000)
- `R`: Number of bootstrap replicates (default: 500)
- `rerun_sim1` through `rerun_sim4`: Control which scenarios to run

## Key Simulations

### Scenario 1: Null Effect (VE = 0%)
Demonstrates bias detection when no vaccine effect exists, testing the ability of different methods to correctly identify null effects.

### Scenario 2: Moderate Effect (VE = 40%)
Evaluates estimation accuracy and bias for a realistic vaccine effectiveness level, comparing naive approaches with target trial emulation.

### Scenario 3: Time-Varying Effectiveness
Models scenarios where vaccine effectiveness depends on the timing of administration post-exposure, using a logistic decay function.

### Scenario 4: Heterogeneous Effects
Examines performance under varying incubation period distributions and different vaccine effectiveness levels (0%, 40%, 80%).

## Statistical Methods

The project implements and compares several analytical approaches:

1. **Naive Methods**:
   - "Leave when vaccinated" approach
   - "Move when vaccinated" approach

2. **Time-to-Event Methods**:
   - Time-varying Cox proportional hazards models

3. **Target Trial Emulation**:
   - Sequential trial emulation avoiding immortal time bias
   - Proper alignment of time zero, eligibility, and treatment assignment

## Output Files

### Tables
- `2_tables/sim_results_hr.tex`: Simulation results for hazard ratio estimates
- `2_tables/sim_results_rr.tex`: Simulation results for relative risk estimates

### Figures
- Simulation performance comparisons across different scenarios
- Bias, coverage, and efficiency metrics for each method
- Distribution plots showing incubation periods and vaccination timing
- Target trial emulation diagnostic plots

## Citation

If you use this code or methodology, please cite:

```
Boyer, C. & Lipsitch, M. (2025). "Emulating target trials of postexposure vaccines using 
observational data". American Journal of Epidemiology, 194(7), 2037-2046. 
https://doi.org/10.1093/aje/kwae350
```

## Contact

For questions about the code or analysis, please contact:
- Christopher Boyer: cboyer@hsph.harvard.edu

## License

This project is licensed under the terms specified in the LICENSE file.
