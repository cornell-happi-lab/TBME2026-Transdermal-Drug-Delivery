
MATLAB analysis of bovine serum albumin (BSA) release from PPy-containing hydrogels in Franz diffusion cells. Per-sample BSA concentrations are corrected for sampling and reported as cumulative permeation **Q (µg/cm²)** versus time for each hydrogel geometry, under iontophoretic (**CIP**) and passive (**PASS**) delivery. The script also reports steady-state flux, lag time and apparent permeability.

## Contents

1. [Files](#files)
2. [Quick start](#quick-start)
3. [Requirements](#requirements)
4. [Input data](#input-data)
5. [Pipeline](#pipeline)
6. [Outputs](#outputs)
7. [Settings reference](#settings-reference)
8. [Data notes](#data-notes)
9. [Known issues and caveats](#known-issues-and-caveats)

---

## Files

| File | What it is |
|---|---|
| [`BSA_Drug_Delivery_Results_Cam_Ready.m`](BSA_Drug_Delivery_Results_Cam_Ready.m) | The analysis script. It prints all results to the console and makes no figures. |
| [`NNLS_Cam_Ready_results.csv`](NNLS_Cam_Ready_results.csv) | Per-sample BSA concentrations from the UV-Vis measurements (180 samples). This is the script's only input. |
| `LICENSE` | License |

---

## Quick start

```matlab
cd('<path to this repository>')
BSA_Drug_Delivery_Results_Cam_Ready
```

From VS Code, open the file and press **F5** (or run **MATLAB: Run File**). The script finds the CSV in its own folder, so it runs from any working directory.

It prints three result blocks: Q at 180 min, then the steady-state flux, lag time and Papp, first as mean ± SD across replicates and then as fit standard errors. Nothing is written to disk.

---

## Requirements

| Requirement | Why |
|---|---|
| MATLAB (developed on R2024b) | Uses tables and `readtable` (R2019a or later). No toolboxes needed. |
| `NNLS_Cam_Ready_results.csv` | Must be in the same folder as the script. |

---

## Input data

`NNLS_Cam_Ready_results.csv` has one row per receptor sample:

```
file,geom,t_min,isPass,Vr_mL,C_BSA_ugmL
T0_B3C.csv,C,0,0,9.850898,-2.1199809779715002
T180_B4H2_Pass_DIL.csv,H2,180,1,9.850898,418.74182540075049
```

| Column | Meaning |
|---|---|
| `file` | The UV-Vis spectrum file the row came from (see naming below) |
| `geom` | Hydrogel geometry (see table below) |
| `t_min` | Sampling time: 0, 15, 30, 60, 120 or 180 min |
| `isPass` | 1 = passive (PASS); 0 = iontophoretic (CIP) |
| `Vr_mL` | Receptor volume (see [Receptor volume](#receptor-volume)) |
| `C_BSA_ugmL` | BSA concentration in the receptor, µg/mL, already corrected for dilution |

Values are written at full precision so the results reproduce exactly. That precision isn't meaningful for the measurement.

### Geometries and delivery modes

| Code | Printed as | Geometry |
|---|---|---|
| `R` | `Plain` | Plain |
| `C` | `C` | Circle |
| `I` | `I` | Interdigital |
| `H2` | `H2` | Hilbert 2nd order |
| `H3` | `H3` | Hilbert 3rd order |

| Code | Delivery |
|---|---|
| CIP (`isPass` = 0) | Iontophoretic (current-driven) |
| PASS (`isPass` = 1) | Passive diffusion |

Each geometry × delivery mode has **3 replicate Franz cells**, sampled at 6 time points: 10 groups × 3 cells × 6 = 180 samples.

### Sample file names

```
T<minutes>_B<cell#><GEOM>[_Pass][_DIL[n]].csv

T0_B1R.csv               t = 0 min,   cell B1, Plain,  CIP,  undiluted
T120_B3C_DIL.csv         t = 120 min, cell B3, Circle, CIP,  diluted
T180_B4H2_Pass_DIL.csv   t = 180 min, cell B4, H2,     PASS, diluted
```

| Part | Meaning |
|---|---|
| `T<n>` | Sampling time in minutes |
| `B<n><GEOM>` | Franz cell number + geometry. A cell is identified by geometry + delivery mode + `B<n>`. |
| `_Pass` | Passive delivery. Files without it are iontophoretic. |
| `_DIL`, `_DIL2`, … | The sample was diluted before measurement. Higher numbers are re-dilutions of the same time point. Concentrations in the CSV are already corrected for dilution. |

### Receptor volume

Every cell was fitted with an adapter. The two adapter shapes give different receptor volumes, and `Vr_mL` follows the adapter actually used:
- 10.255033 mL: `Rec` adapter;
- 9.850898 mL: `Cir` adapter, used for all Circle (C) samples and for the B4 PASS cells of H2, H3 and I.

### Sampling

At each time point, 1 mL is drawn from the receptor and replaced with fresh PBS. The analysis adds the BSA removed by earlier draws back in when it builds Q(t).

---

## Pipeline

```
NNLS_Cam_Ready_results.csv
                │
group by geometry × condition × time (mean ± SD)
                │
sampling-corrected cumulative conc ─► Q (µg/cm²)
                │
console: Q(180) + steady-state Jss / τ / Papp (SD and SE tables)
```

### Stage 1: Cumulative release with sampling correction

For each geometry × condition, the replicates at each time point are averaged (mean ± SD). Each sampling removes `Vs = 1 mL` from the receptor, so the cumulative amount at sample *j* is:

```
M_j      = Vr · C̄_j + Vs · Σ_{i<j} C̄_i
Var(M_j) = Vr² · s_j² + Vs² · Σ_{i<j} s_i²
C_cum,j  = M_j / Vr
Q_j      = Vr · C_cum,j / A        (A = 1.7549275 cm²)
```

`Vr` is the median receptor volume of the group. Error bars are the propagated SD.

### Stage 2: Derived metrics

- **Q(180):** the cumulative Q at t = 180 min, or the nearest time point (noted in the output).
- **Steady state.** The script finds the straightest stretch of Q(t), fits a line to it, and calculates flux, lag time and Papp from that line. It does this in two ways:

  **a. Per-replicate mean ± SD** (printed first; reported in the manuscript's summary table)
  1. Each of the 3 cells gets its own Q(t), built from its own samples with its own receptor volume.
  2. Each cell gets its own steady-state window (`pick_ss_window_longest`): the **longest** run of consecutive time points whose linear fit has R² ≥ **0.957**, with ties going to the higher R². 2-point windows are allowed, so a cell with no qualifying 3-point window still gets a fit (R² = 1). This fallback was used for one cell each in C PASS and C CIP, and two cells in H3 CIP.
  3. Jss, τ, Papp and R² are computed per cell, then reported as mean ± SD across the 3 cells (n − 1). The printout lists each cell's window.

  **b. Fit on the averaged curve, ± standard error** (printed second; reported in the Supplementary Information)
  1. Each group's averaged Q(t) is used.
  2. Every window of 3 or more consecutive points is fitted with `polyfit`. Of the windows with R² > 0.957, the one with the **highest R²** wins, and a tie goes to the longer window (`pick_ss_window`).
  3. Slope and intercept SEs come from the `polyfit` covariance. The SE of τ is propagated from both by the delta method, including their covariance. Jss and Papp SEs scale from the slope SE.

  Because the window rules and the data differ (3 separate cells vs one averaged curve), the two tables give different values.

  Both tables calculate:

  | Quantity | Formula | Units |
  |---|---|---|
  | Steady-state flux | `Jss = slope × 60` | µg cm⁻² h⁻¹ |
  | Lag time | `τ = max(−intercept / slope, 0)`: where the line crosses Q = 0, clipped at 0 | min |
  | Apparent permeability | `Papp = Jss / C0 / 3600`, with `C0` = 20,000 µg/mL (donor loading) | cm/s |

  **Interpreting the results**
  - **"Most linear" is not "latest".** With only 6 time points, many windows fall on the early burst (0–30 or 0–60 min) rather than a classic post-lag linear phase.
  - **Negative lag times are clipped to 0.** A negative τ means the line crosses Q = 0 before t = 0, which is typical of an initial burst followed by a slower linear phase. A lag time can't be negative, so the script prints 0. In the SE table, the ± shown is still the SE of the unclipped estimate. In the SD table each cell is clipped before averaging, so a group whose cells straddle 0 can show a positive mean with a large SD.
  - **Short windows.** 3-point windows have 1 degree of freedom, so their SEs are rough. 2-point windows have R² = 1 by construction.
  - **Unweighted fits.** Each time point counts equally, even though the replicate spread grows over time.

---

## Outputs

| Block | Contents |
|---|---|
| `Q (cumulative amount) at 180 min` | Q ± SD for each geometry × condition |
| `Steady-state parameters across replicate Franz cells` | For each group: n cells, each cell's window, and mean ± SD of Jss, τ, Papp and R² |
| `Derived steady-state parameters from cumulative Q(t)` | Window, slope, Jss, τ, Papp and R², each ± its standard error (reported in the Supplementary Information) |

---

## Settings reference

The settings are at the top of the script, apart from `r2_thresh`, which is in the steady-state section.

| Variable | Value | Meaning |
|---|---|---|
| `Vs` | `1.0` mL | Aliquot volume removed (and replaced with PBS) at each sampling |
| `A_cm2` | `1.7549275` cm² | Diffusion area |
| `C0_ugmL` | `20000` µg/mL | Donor (hydrogel) BSA concentration, used for Papp |
| `r2_thresh` | `0.957` | Minimum R² for a steady-state window |
| `nnlsCsv` | `NNLS_Cam_Ready_results.csv` in the script's folder | Per-sample input (180 rows) |

---

## Data notes

All 180 samples are used as measured; none were adjusted.

- **Below detection.** A value of −2.12 µg/mL means no BSA was detected. It is what the calibration line, (0 − 0.001112) / 0.000524533, gives for a zero BSA signal. 
- **Non-monotonic readings.** Sampling replaces 1 mL of the ~10 mL receptor with fresh PBS, so readings can dip slightly between time points. A few cells show larger decreases at later time points. In the cumulative curves these show up as reduced or negative increments.
- **Baseline readings.** Some cells read slightly above zero at t = 0.
- **Receptor volumes.** Receptor volume depends on which adapter was used: 10.255033 mL (`Rec`) or 9.850898 mL (`Cir`). See [Receptor volume](#receptor-volume). The per-replicate SD table uses each cell's own volume. Q(180) and the SE table use the group median.

---

## Known issues and caveats

- **Volume choice.** The group-level cumulative curves (Q(180) and the SE table) use the **median** receptor volume of each group, not each sample's own volume.
- **Display name.** The script renames geometry `R` to `Plain`, for better readability. The CSV and file names keep `R`.
