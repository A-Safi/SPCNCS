# Dynamic Predictive Control for Networked Control Systems (Delay + Packet Loss)

This repository contains MATLAB code accompanying the paper:

**“The design of dynamic predictive control for networked control systems subject to latency and packet loss”**  
Ali Safi, Mostafa Nasiri, Reza Farasat

The project implements a **dynamic predictive control** scheme for **networked control systems (NCSs)** operating over a **UDP-like communication channel** subject to:
- **Random measurement latency** (time-varying delay), and
- **Packet dropout** (data loss without acknowledgment).

The controller is based on an **LQTI (Linear Quadratic Tracker with Integrator)** structure and extends classical networked predictive control by predicting **both**:
1. the plant state evolution, and  
2. the controller internal dynamics (integrator state),

so that stable regulation and tracking can be maintained under delays and packet losses.

---

## Repository contents (expected)

Typical files used by the simulation include:

- `*.m` MATLAB scripts and functions implementing:
  - plant/controller discretization
  - predictive estimation logic under delay/dropout
  - simulation and plotting
- `StateSpace.mat` plant state-space matrices (must include at least `A` and `B`)
- `Inverted_Pendulum.m` nonlinear pendulum dynamics (used by the integrator)
- `RungeKutta.m` Runge–Kutta integration routine

If any of these are named differently in your repo, update the filenames in the sections below.

---

## Method summary

### Plant model
The plant is modeled in discrete time as:
\[
x_p(k+1) = A_p x_p(k) + B_p u_p(k), \quad y_p(k)=C_p x_p(k)
\]

### Dynamic LQTI controller
The controller is represented in discrete-time state-space form:
\[
x_c(k+1)=A_c x_c(k)+B_c u_c(k), \quad y_c(k)=C_c x_c(k)+D_c u_c(k)
\]

The LQTI structure embeds integral action to reduce steady-state tracking error.

### Network effects
A UDP-like network is emulated such that:
- each measurement packet experiences a random delay, and
- packets may be dropped with a probability dependent on a single quality parameter `M`.

The parameter **M** is interpreted as a network-quality proxy:
- Minimum delay: `M`
- Maximum delay: `2M`
- Drop probability: `10M%`

### Prediction logic
When measurements are delayed or lost, the predictor:
- reuses the most recent available measurement if no new packet arrives, and
- increases the prediction depth according to the number of missing samples.

The simulation includes two parallel branches:
1. **Baseline**: no packet loss / no predictive compensation logic (reference)
2. **Predictive NCS**: delayed + dropped measurements handled via prediction

---

## Requirements

- MATLAB (Control System Toolbox recommended)
- `StateSpace.mat` available in the working directory

Optional (for stability LMIs discussed in the paper):
- YALMIP toolbox
- SDP solver (e.g., SDPT3)

---

## Quick start

1. Place all MATLAB scripts/functions in a single folder.
2. Ensure `StateSpace.mat` exists in the same folder and contains:
   - `A` (continuous-time state matrix)
   - `B` (continuous-time input matrix)
3. Run the main simulation script from MATLAB.

If you are using the refactored script that runs all cases automatically, it will simulate:

- `M = 1`
- `M = 2`
- `M = 3`
- `M = 4`

and produce:
- a dropout/delay plot per `M`, and
- a response comparison plot per `M`.

---

## Key parameters

Typical parameters used in the simulations:
- Sampling time: `T = 0.01 s`
- Simulation time: e.g. `SimTime = 3 s`
- Initial state example: `x0 = [0; 0; 0.04; 0]`

Network severity is swept by `M`:
- `M=1`: low delay / low loss
- `M=4`: high delay / high loss (baseline often becomes unstable)

---

## Reproducing paper-style results

To replicate the figures that compare “with prediction” vs “without prediction”:
1. Run the simulation for `M=1..4`.
2. For each `M`, inspect:
   - state trajectories (arm angle `θ`, pendulum angle `α`)
   - control torque command

The expected qualitative trend reported in the paper:
- predictive controller maintains stability and acceptable performance as `M` increases,
- non-predictive controller becomes oscillatory and may destabilize as delay/loss increases.

---

## Notes on stability analysis (LMIs)

The paper models the overall closed-loop networked system as a **discrete-time switched system**, with switching driven by delay/dropout events. Stability can be certified by a **Common Quadratic Lyapunov Function (CQLF)** leading to LMIs of the form:

\[
\Lambda_\delta^T P \Lambda_\delta - P < 0,\quad P>0
\]

Because these matrices can become large for bigger `M`, the repository focuses primarily on simulation; the LMI derivations and solver setup may be provided separately (or can be added via YALMIP scripts if desired).

---

## Citation

If you use this code in academic work, please cite the associated paper:

Ali Safi, Mostafa Nasiri, Reza Farasat,  
**“The design of dynamic predictive control for networked control systems subject to latency and packet loss.”**

---

## License

Add your preferred license here (e.g., MIT, BSD-3-Clause, GPL-3.0). If no license is included, default copyright restrictions apply.

---

## Contact

- Ali Safi: ali_safi@alumni.iust.ac.ir  
- Mostafa Nasiri: m.nasiri@iut.ac.ir  
- Reza Farasat: r.ferasat@gut.ac.ir
