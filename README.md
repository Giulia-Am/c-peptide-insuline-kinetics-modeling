# Pharmacokinetic & Enzymatic Kinetics Modeling Suite (MATLAB)

This repository contains a collection of MATLAB scripts developed for the analysis, simulation, and parameter estimation of pharmacokinetic (PK) and enzymatic kinetics models. The scripts cover liear and nonlinear compartmental simulation, input-output parameter estimation, model identification, enzymatic reaction kinetics, population analysis,multiple regression and deconvolution-based estimation of insulin secretion.  
## Folder Structure

All files — main scripts, custom function files (.m), and datasets/data-loading scripts — are kept together in a single flat folder, with no subfolder separation between scripts, functions, and data. This is intentional: since MATLAB only needs every file to be on the current path, keeping everything in one directory is the simplest way to guarantee that all scripts run without additional configuration.

Simply download/clone the entire folder as-is and set it as your MATLAB current working directory (or add it to the path) — no reorganization is needed.

---

## 1. Main Scripts

| # | Script (topic) | Description |
|---|---|---|
| 1 | **PK Models – LTI Systems** | One- and multi-compartment PK models (IV bolus, extravascular absorption, two- and three-compartment models), Bode analysis, impulse response, non-compartmental parameters (AUC, AUMC, MRT, CL), optimal sampling design, repeated dosing simulation. |
| 2 | **Nonlinear Compartmental Models** | Comparison of linear vs. nonlinear (Michaelis-Menten) absorption kinetics, sensitivity analysis, and repeated dosing under nonlinear absorption. |
| 3 | **C-Peptide I/O Estimation (v1)** | Weighted Least Squares (WLS) fitting of bi-exponential C-peptide kinetics across 7 subjects, plus population analysis (NAD, NPD, STS). |
| 4 | **C-Peptide I/O Estimation (v2)** | Compartmental (mechanistic) LS/WLS estimation of k01, k12, k21, V1; model comparison, Monte Carlo simulation, and Bootstrap uncertainty analysis. |
| 5 | **Enzymatic Kinetics Simulation** | ODE-based simulation of enzyme-substrate kinetics under different enzyme/substrate ratios, carbonic anhydrase kinetics, and competitive/uncompetitive/non-competitive inhibition mechanisms. |
| 6 | **Multiple Linear Regression (Van Cauter Model)** | Population-based regression of PK parameters (Volume, half-lives, Fraction) against anthropometric covariates, stepwise variable selection, confidence intervals, and individual parameter prediction. |
| 7 | **Deconvolution ** | Reconstruction of insulin secretion from plasma C-peptide via raw deconvolution and regularized deconvolution (Twomey, GCV, Maximum Likelihood criteria), plus Bayesian estimation (MCMC). |


> Scripts 1 and 7 are fully self-contained: their helper functions (`optimal_sampling_obj`, `fun_MM`) are defined as local functions at the bottom of the file and require no additional setup.

---

## 2. Required Custom Functions

The following user-defined `.m` function files **must be present in the working directory** for the corresponding scripts to run. These are not built-in MATLAB functions.

| Function file | Required by | Purpose |
|---|---|---|
| `residui.m` | C-Peptide I/O Estimation (v1) | Computes weighted residuals for the bi-exponential model, used by `lsqnonlin`. |
| `residui_new.m` | C-Peptide I/O Estimation (v2) | Computes weighted residuals for the compartmental model (parameters k01, k12, k21, V1). |
| `discrepanza.m` | Deconvolution (Exercise 7) | Objective function implementing Twomey's discrepancy criterion for regularization parameter (γ) selection. |
| `discrepanza_post.m` | Deconvolution (Exercise 7) | Computes the regularized solution, RSS, and residuals for a given γ. |
| `crossvalidation.m` | Deconvolution (Exercise 7) | Objective function implementing the Generalized Cross-Validation (GCV) criterion. |
| `mlikelihood.m` | Deconvolution (Exercise 7) | Objective function implementing the Maximum Likelihood criterion for γ selection. |
| `fun_enzimi.m` | Enzymatic Kinetics Simulation | ODE system for basic enzyme-substrate kinetics (with product recombination). |
| `fun_enzimi_co2.m` | Enzymatic Kinetics Simulation | ODE system for carbonic anhydrase kinetics (no recombination, k₋₂ = 0). |
| `fun_enzimi_co2_competitiva.m` | Enzymatic Kinetics Simulation | ODE system for competitive inhibition. |
| `fun_enzimi_co2_anticompetitiva.m` | Enzymatic Kinetics Simulation | ODE system for uncompetitive inhibition. |
| `fun_enzimi_co2_noncompetitiva.m` | Enzymatic Kinetics Simulation | ODE system for non-competitive inhibition. |
| `stepwise.m` | Multiple Linear Regression (Van Cauter Model) | Custom forward stepwise regression routine (returns R², selected coefficients, selected regressor matrix, and selected indices). **Not** the built-in Statistics Toolbox function — signature and outputs are project-specific. |

---

## 3. Required Datasets

| File | Required by | Description |
|---|---|---|
| `DatiCPsog1.dat` … `DatiCPsog7.dat` (7 files) | C-Peptide I/O Estimation (v1 and v2) | Time-concentration measurements of synthetic C-peptide for 7 individual subjects. Expected format: two columns (time, concentration), 4 header lines, tab-delimited. |
| `dati_es7.m` | Deconvolution (Exercise 7) | Data-loading script; must define the variables `DATA` (time/concentration matrix) and `yb` (basal concentration) used for the deconvolution analysis. |
| `dati_vc.m` | Multiple Linear Regression (Van Cauter Model) | Data-loading script; must define the variable `dati`, a 207×11 matrix containing 7 covariate columns (health status, gender, age, height, weight, BMI, BSA) followed by 4 kinetic parameter columns (Volume, short half-life, long half-life, Fraction). |

---

## 4. Requirements

- **MATLAB** (developed and tested on a recent release)
- **Control System Toolbox** — `ss`, `tf`, `impulse`, `lsim`, `residue`, `bode`
- **Optimization Toolbox** — `lsqnonlin`, `fmincon`, `fminsearch`
- **Statistics and Machine Learning Toolbox** — `tinv`, `trnd`, `prctile`, `kstest`, `makedist`, `mvnrnd`, `gamrnd`, `cov`

---

## 5. Usage Notes

- Some scripts (C-Peptide estimation, Van Cauter regression) prompt for user input via the MATLAB console (`input(...)`) to select a subject number or enter personal anthropometric parameters — run these interactively, not as background jobs.
- All figures are generated automatically; no manual plotting steps are required.
- Scripts assume a clean workspace and will run `clear`/`close all` at the start — save any unrelated work before executing.
- Random-seed-dependent sections (Monte Carlo, Bootstrap, Bayesian/MCMC estimation) will produce slightly different numerical results on each run unless a fixed seed (`rng(...)`) is set beforehand.

---

## 6. Quick Start

1. Keep all files — scripts, functions, and data — in the same folder (no reorganization needed).
2. Open MATLAB and set this folder as the current working directory (or add it to the path).
3. Run the desired script (e.g., type its filename without the `.m` extension in the Command Window, or open it in the Editor and press **Run**).
4. Follow any console prompts (subject number, anthropometric data) when requested.

