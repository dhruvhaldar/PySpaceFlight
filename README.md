# PySpaceFlight: SD2900 Fundamentals of Spaceflight Simulator

![Python](https://img.shields.io/badge/python-3.x-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-active-brightgreen.svg)

## Project Overview

**PySpaceFlight** is a comprehensive Python-based simulation tool developed for **[ KTH Royal Institute of Technology : SD2900 Fundamentals of Spaceflight](https://www.kth.se/student/kurser/kurs/SD2900?l=en)** course . It demonstrates key concepts from the syllabus, including:

*   **Launcher Dynamics**: Simulation of a multi-stage rocket launch from Earth, accounting for thrust, aerodynamic drag, gravity losses, and mass variation.
*   **Orbital Mechanics**: Calculation of Keplerian orbital elements (semi-major axis, eccentricity, inclination, etc.) from state vectors.
*   **Orbital Maneuvers**: Computation of Delta-V requirements for Hohmann transfers (e.g., Low Earth Orbit to Geostationary Orbit).

The simulator numerically integrates the equations of motion for a rocket ascent and visualizes the resulting trajectory and orbit.

## Features

*   **Multi-Stage Rocket Simulation**: Configurable rocket stages with specific impulse (ISP), thrust, mass, and aerodynamic properties.
*   **Physics Engine**: Includes a standard atmosphere model, gravity variation with altitude, and drag calculations.
*   **Guidance Algorithm**: Implements a gravity turn maneuver to achieve orbit.
*   **Orbital Analysis**: Automatically calculates orbital elements upon engine cutoff.
*   **Maneuver Planning**: Calculates the necessary burns to transfer from the parking orbit to a target Geostationary Orbit (GEO).
*   **Visualization**: Generates detailed plots of the flight profile and orbital geometry.

## Syllabus Coverage

This project directly addresses the following Intended Learning Outcomes (ILOs) from the SD2900 syllabus:

*   *Demonstrate basic methodology and understanding of spaceflight, including launcher dynamics, orbital mechanics, manoeuvres and relative motion in orbit.*
*   *Demonstrate the ability to model, simulate, predict and evaluate spacecraft behaviour from launching to rendezvous with other spacecraft.*
*   *Demonstrate the ability to analyse and critically evaluate various technical solutions for a geocentric space mission.*

## Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/dhruvhaldar/PySpaceFlight
    cd PySpaceFlight
    ```

2.  Create and activate a virtual environment:

    *   **Windows:**
        ```bash
        python -m venv venv
        .\venv\Scripts\activate.ps1
        ```

    *   **macOS/Linux:**
        ```bash
        python3 -m venv venv
        source venv/bin/activate
        ```

3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Run the main simulation script:

*   **Windows:**
    ```powershell
    python .\src\main.py
    ```

*   **macOS/Linux:**
    ```bash
    python3 ./src/main.py
    ```

This will execute the simulation and generate the following plots in the `images/` directory:

### Visualizations

1.  **Ascent Profile**: Altitude, Velocity, Dynamic Pressure, and Acceleration vs. Time.
    ![Ascent Profile](images/ascent_profile.png)

2.  **Trajectory**: 2D Ascent Trajectory (Altitude vs. Downrange Distance).
    ![Trajectory](images/trajectory.png)

3.  **Orbit Visualization**: Visualization of Earth, Parking Orbit, Transfer Orbit, and GEO.
    ![Orbit Visualization](images/orbit_viz.png)

### Example Output

```text
Starting simulation...
Stage separation! Ignition of stage 2 at t=153.00s
Mission complete (fuel exhausted or target reached)

Final Orbital Elements:
a: 10681311.7873
e: 0.3862
i: 180.0000
Omega: 0.0000
omega: 110.3479
nu: 60.8394
period: 10986.3643

Calculating Hohmann Transfer from 1277.91 km to GEO (35786 km)...
Delta-V 1 (LEO burn): 2173.54 m/s
Delta-V 2 (GEO circularization): 1370.75 m/s
Total Delta-V: 3544.29 m/s
Transfer Time: 5.43 hours
```

## Structure

*   `src/`: Source code.
    *   `main.py`: Main simulation script.
    *   `rocket.py`: Rocket physics and dynamics classes.
    *   `orbit.py`: Orbital mechanics calculations.
*   `tests/`: Unit tests.
    *   `test_physics.py`: Verification of physics calculations.
*   `images/`: Generated plots.

## License

MIT License
