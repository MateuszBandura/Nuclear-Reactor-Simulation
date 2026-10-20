# Nuclear Reactor Simulation

Interactive web application for simulating nuclear reactor point-kinetics dynamics with thermal feedback.

## Overview

This project provides a Dash-based simulator that models:
- neutron density evolution,
- delayed neutron precursor groups,
- core temperature response,
- external reactivity insertion and temperature feedback.

The UI allows selecting reactor and thermal parameters, running repeated simulations, and comparing current results with the previous run.

## Features

- Three predefined reactor types (`PWR (A)`, `Badawczy (B)`, `Szybki (C)`)
- Configurable nominal power, thermal time constant, temperature coefficient, inlet temperature, and temperature rise
- Reactivity control via slider (`Δk/k`)
- Automatic state carry-over between runs
- Manual reset to original initial conditions
- Interactive plots:
  - Reactor power vs time
  - Core temperature vs time
  - Relative delayed neutron precursor concentrations

## Technology Stack

- Python
- Dash + Plotly
- NumPy
- SciPy (`solve_ivp`)
- Dash Bootstrap Components

## Project Structure

```text
.
├── main.py      # Application code and reactor model
└── README.md    # Project documentation
```

## Installation

1. Create and activate a virtual environment (recommended).
2. Install dependencies:

```bash
pip install dash dash-bootstrap-components plotly numpy scipy
```

## Running the Application

From the project root:

```bash
python main.py
```

Then open the local Dash URL shown in the terminal (typically `http://127.0.0.1:8050`).

## How to Use

1. Choose reactor and thermal parameters from dropdowns.
2. Set external reactivity using the slider.
3. Click **Uruchom Symulację** to run simulation.
4. Review updated plots and resulting initial values for the next run.
5. Click **Reset** to return to original conditions.

## Model Notes

- The simulation uses point-kinetics equations with six delayed neutron precursor groups.
- Temperature feedback is modeled as:
  - `ρ_total = ρ_ext + α_T (T - T_ref)`
- ODE integration is performed with SciPy’s `Radau` solver.