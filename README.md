# Quantum-Oscillator Path-Integral Monte Carlo

This project uses Euclidean path-integral Monte Carlo (PIMC) to study three
one-dimensional quantum systems:

- the harmonic oscillator;
- a quartic anharmonic oscillator;
- a symmetric double-well oscillator.

The program samples periodic paths in imaginary time, estimates ground-state
energies and probability densities, extracts first excited-state energies from
time-displaced correlation functions, and measures thresholded Euclidean
instanton densities in the double wells.

Every reported energy is compared with an independent finite-difference
solution of the same Schrödinger Hamiltonian. The reference calculation is
used only for validation: it does not set the Monte Carlo result, the fit
window, or the initial energy-gap estimate.

The default study demonstrates:

- discretisation of the Euclidean action on a periodic time lattice;
- local Metropolis, free-particle staging, centroid, and reflection moves;
- independent reproducible random streams;
- virial energy estimation;
- autocorrelation-aware block errors;
- periodic correlation-function fits;
- block-bootstrap spectral uncertainties;
- ground-state projection and lattice-spacing checks;
- the distinction between Euclidean instantons and real-time dynamics.

## Repository Contents

```text
.
├── .gitignore
├── README.md
├── quantum_oscillator_pimc.py
└── requirements.txt
```

The numerical study is self-contained and requires no external dataset. The
original coursework script, lecturer-provided assignment sheet, poster
material, generated figures, caches, and verification outputs are deliberately
excluded.

## Requirements

The project has been verified using:

| Component | Tested version |
|---|---:|
| Python | 3.12.13 |
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| Matplotlib | 3.10.8 |

The compatible package ranges are recorded in `requirements.txt`.

## Installation

From the repository directory, create and activate a virtual environment.

On macOS or Linux:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

## Running the Study

Run the full publication configuration from the repository root:

```bash
python quantum_oscillator_pimc.py
```

This prints the production and convergence tables and displays five
Matplotlib figures.

For a non-interactive run:

```bash
python quantum_oscillator_pimc.py --no-show
```

To save the figures as PNG files without opening interactive windows:

```bash
python quantum_oscillator_pimc.py --no-show --save-dir figures
```

For a reduced diagnostic run:

```bash
python quantum_oscillator_pimc.py --quick --no-show
```

In the verified environment, the full deterministic run completed in
approximately two minutes and the quick run completed in approximately 24
seconds. Runtime depends on the processor and the installed NumPy and SciPy
builds. The results documented below come from the full configuration, not
from `--quick`.

## Physical Background

In units where $hbar=1$, the Hamiltonian of a one-dimensional particle is

$$
H
=
\frac{p^2}{2m}+V(x).
$$

The project considers three even potentials.

For the harmonic oscillator,

$$
V_{\mathrm{H}}(x)=\frac{1}{2}kx^2.
$$

For the quartic anharmonic oscillator,

$$
V_{\mathrm{A}}(x)
=
\frac{1}{2}kx^2+\lambda x^4,
\qquad \lambda>0.
$$

For the symmetric double well,

$$
V_{\mathrm{D}}(x)
=
ax^4+bx^2,
\qquad a>0,
\qquad b<0.
$$

The harmonic-oscillator angular frequency and exact spectrum are

$$
\omega=\sqrt{\frac{k}{m}},
$$

and

$$
E_n
=
\left(n+\frac{1}{2}\right)\omega.
$$

Its ground-state probability density is

$$
|\psi_0(x)|^2
=
\sqrt{\frac{m\omega}{\pi}}
e^{-m\omega x^2}.
$$

The anharmonic and double-well spectra do not have equally simple general
closed forms. The program therefore also constructs and diagonalises a
finite-difference Hamiltonian for all nine parameter sets.

## Euclidean Path Integral

At inverse temperature, or total Euclidean extent, $\beta$, the partition
function is

$$
Z(\beta)
=
\operatorname{Tr}\left(e^{-\beta H}\right).
$$

After continuing to imaginary time $\tau$, it can be written as an integral
over periodic paths:

$$
Z(\beta)
=
\int_{x(\beta)=x(0)}\mathcal{D}x\;e^{-S_E[x]}.
$$

The Euclidean action is

$$
S_E[x]
=
\int_0^\beta
\left[
\frac{m}{2}
\left(\frac{dx}{d\tau}\right)^2
+V(x)
\right]d\tau.
$$

The factor $e^{-S_E}$ is real and non-negative for the potentials used here,
so it can be sampled as a probability weight.

For sufficiently large $\beta$, excited-state contributions to
$e^{-\beta H}$ are exponentially suppressed. The ensemble then approaches
the quantum ground state rather than the classical distribution
$e^{-\beta V(x)}$.

## Lattice Discretisation

The interval $[0,\beta)$ is divided into $N$ sites with spacing

$$
a_\tau=\frac{\beta}{N}.
$$

With periodic indexing $x_N=x_0$, the primitive lattice action is

$$
S_E^{(N)}
=
\sum_{j=0}^{N-1}
\left[
\frac{m}{2a_\tau}(x_{j+1}-x_j)^2
+a_\tau V(x_j)
\right].
$$

The first term is essential. It couples neighbouring imaginary-time sites and
distinguishes a quantum path ensemble from independent scalar sampling.

The default production lattice uses

$$
\beta=16,
\qquad
N=256,
\qquad
a_\tau=0.0625.
$$

## Main Aim

The numerical study tests whether:

1. periodic-path sampling reproduces known harmonic-oscillator results;
2. PIMC energies agree with an independent Schrödinger solver for all three
   potential families;
3. sampled time-slice distributions approach $|\psi_0(x)|^2$;
4. the first spectral gap can be extracted from imaginary-time correlations;
5. autocorrelation-aware errors describe the observed agreement;
6. the default results remain stable under changes in lattice spacing and
   Euclidean extent;
7. double-well paths contain resolvable Euclidean instanton transitions.

## What the Code Does

For each of the nine production parameter sets, the program:

1. creates four statistically independent deterministic random streams;
2. constructs a periodic imaginary-time path for each chain;
3. performs 2,000 thermalisation sweeps;
4. adapts the local proposal width only during the early thermalisation stage;
5. freezes the proposal before measurements begin;
6. combines local, staging, centroid, and reflection updates;
7. records 2,000 configurations per chain, separated by four compound sweeps;
8. estimates the virial energy on every saved path;
9. computes time-origin-averaged periodic correlations;
10. estimates the integrated energy-autocorrelation time;
11. blocks the measurements and verifies that the blocks exceed the measured
    correlation scale;
12. fits the first energy gap with a periodic spectral model;
13. bootstraps the gap and first excited-state energy;
14. samples the ground-state probability density;
15. calculates a finite-difference Hamiltonian benchmark;
16. counts thresholded Euclidean instantons for double wells;
17. prints reproducible comparison tables;
18. performs separate lattice-spacing and Euclidean-extent checks;
19. creates five diagnostic figures.

## Monte Carlo Updates

### Local checkerboard moves

For one lattice site, the program proposes

$$
x_j'=x_j+\delta,
$$

where

$$
\delta\sim\mathcal{U}(-w,w).
$$

Only the two adjacent kinetic links and the local potential term change. The
proposal is accepted with probability

$$
P_{\mathrm{accept}}
=
\min\left(1,e^{-\Delta S_E}\right).
$$

Even sites and odd sites are updated in separate vectorised passes. Sites of
one parity do not share a direct link, so their proposals can be evaluated
together using the opposite parity as fixed neighbours.

### Proposal adaptation

The initial proposal width is scaled with the free kinetic fluctuation:

$$
w_0
=
1.8\sqrt{\frac{a_\tau}{m}}.
$$

During the early thermalisation sweeps, the width is adjusted every 100 sweeps
towards a target local acceptance of $0.50$. Adaptation stops before the final
thermalisation segment and remains frozen throughout measurement. The verified
production acceptance rates lie between $0.500$ and $0.502$.

### Free-particle staging moves

Purely local updates become inefficient as $a_\tau$ decreases because the
kinetic term strongly correlates neighbouring sites. The program therefore
also updates segments of 16 and 32 links using free-particle Brownian bridges.

For fixed segment endpoints, the proposed interior is drawn exactly from the
conditional kinetic distribution. The kinetic contribution cancels in the
Metropolis-Hastings ratio, leaving only the change in potential action:

$$
P_{\mathrm{stage}}
=
\min\left[
1,
\exp\left(
-a_\tau\sum_{j\in\mathrm{segment}}
\left[V(x_j')-V(x_j)\right]
\right)
\right].
$$

These non-local moves update long-wavelength structure more efficiently than
site moves alone.

### Centroid and reflection moves

A centroid proposal translates every time slice by the same random amount:

$$
x_j'=x_j+c.
$$

Because all differences $x_{j+1}-x_j$ are unchanged, only the potential part
of the action enters the acceptance decision.

All three potentials are even. The global reflection

$$
x_j'=-x_j
$$

therefore preserves the full action exactly and is accepted automatically.
It helps balance the two symmetry-related sectors without altering even
observables or creating artificial instantons within a path.

### Independent chains and reproducibility

The default root seed is

```python
seed = 2025
```

`numpy.random.SeedSequence` spawns independent child streams for every chain,
potential, convergence run, and bootstrap analysis. Reusing the complete
configuration and root seed reproduces the documented output exactly.

Changing the seed changes the finite Monte Carlo sample but should not change
the statistical conclusions.

## Ground-State Energy

The code uses the one-dimensional virial estimator. For a stationary state,

$$
\left\langle T\right\rangle
=
\frac{1}{2}
\left\langle xV'(x)\right\rangle.
$$

The total energy estimator on one periodic path is therefore

$$
E_{\mathrm{virial}}
=
\frac{1}{N}
\sum_{j=0}^{N-1}
\left[
V(x_j)+\frac{1}{2}x_jV'(x_j)
\right].
$$

For the harmonic potential this reduces to

$$
E_{\mathrm{H}}
=
k\left\langle x^2\right\rangle.
$$

For the anharmonic potential,

$$
E_{\mathrm{A}}
=
k\left\langle x^2\right\rangle
+3\lambda\left\langle x^4\right\rangle.
$$

For the double well,

$$
E_{\mathrm{D}}
=
2b\left\langle x^2\right\rangle
+3a\left\langle x^4\right\rangle.
$$

Negative ground-state energies are possible for a double-well potential whose
minima lie below zero. An additive shift in the potential would shift every
energy by the same constant without changing the paths' physical structure.

## Correlated Statistical Analysis

### Thermalisation and measurement spacing

No observable is retained during the first 2,000 production sweeps. After
thermalisation, four compound sweeps separate successive saved paths.

These choices reduce, but do not remove, Markov-chain dependence. The program
therefore diagnoses the remaining autocorrelation and analyses blocked rather
than individual measurements.

### Integrated autocorrelation time

For a sequence of virial-energy measurements $E_i$, the normalized Monte Carlo
autocorrelation is

$$
\rho(t)
=
\frac{
\left\langle(E_i-\bar E)(E_{i+t}-\bar E)\right\rangle
}{
\left\langle(E_i-\bar E)^2\right\rangle
}.
$$

The integrated time is estimated as

$$
\tau_{\mathrm{int}}
=
\frac{1}{2}+\sum_{t=1}^{W}\rho(t),
$$

where the summation window ends at the first non-positive adjacent-lag pair.
The autocorrelation itself is evaluated with an FFT.

The production values lie between $1.66$ and $3.52$ saved measurements.

### Blocked energy error

Every chain is divided into blocks of 50 measurements. Block means from all
four independent streams are then combined.

If $B$ is the total number of blocks and $\bar E_b$ is block $b$'s mean, the
reported standard error is

$$
\sigma_{E_0}
=
\frac{s_{\mathrm{block}}}{\sqrt{B}}.
$$

The code requires the block length to exceed five times the measured mean
$\tau_{\mathrm{int}}$. A custom configuration that violates this safeguard
raises a descriptive `RuntimeError` rather than reporting an inadequately
blocked uncertainty.

## Ground-State Probability Density

For a periodic finite-temperature path integral, the marginal distribution of
one time slice is proportional to the diagonal density matrix
$\rho(x,x;\beta)$. In the large-$\beta$ limit,

$$
\frac{\rho(x,x;\beta)}{Z(\beta)}
\longrightarrow
|\psi_0(x)|^2.
$$

Imaginary-time translation invariance allows all lattice sites to contribute
to the histogram. The implementation retains a regular subset of approximately
64 sites from every measured path to control memory use without privileging a
particular time origin.

The plotted histogram is compared with the ground-state eigenvector returned
by the independent Hamiltonian solver.

## First Excited State

### Periodic correlation function

For every measured path, the code calculates the circular, time-origin-averaged
correlation

$$
C(\tau_\ell)
=
\frac{1}{N}
\sum_{j=0}^{N-1}
x_jx_{j+\ell},
$$

with periodic site indices. The calculation uses the convolution theorem and
an FFT.

For an even potential, the position operator $x$ is odd under parity. Its
leading spectral contribution therefore connects the even ground state to the
odd first excited state. At sufficiently large imaginary-time separation,

$$
C(\tau)
\simeq
A
\left[
e^{-\Delta E\tau}
+e^{-\Delta E(\beta-\tau)}
\right],
$$

where

$$
\Delta E=E_1-E_0.
$$

The second exponential is the backward-propagating contribution required by
periodic Euclidean time.

### Fit-window selection

The fitter considers starting points between approximately
$\tau=0.25$ and $\tau=1.5$. The latest positive point more than two standard
errors from zero sets the preferred upper limit. The selector retains a
minimum viable window when necessary and does not fit beyond the first half of
the available $[0,\beta/2]$ range.

Candidate windows are compared for gap stability. The selector uses the
earliest window that is stable against later starts and has an acceptable
diagonal reduced-$\chi^2$ diagnostic. If no candidate meets both conditions,
it selects the valid window with the best diagnostic score.

The optimiser's initial gap comes from early measured log ratios:

$$
\Delta E_{\mathrm{initial}}
=
-\frac{1}{a_\tau}
\log\left[
\frac{C(\tau+a_\tau)}{C(\tau)}
\right].
$$

The numerical Hamiltonian benchmark is not used as an initial value or fit
target.

### Block-bootstrap uncertainty

The same 50-measurement blocks used for the energy are applied to the complete
correlation vector. The default analysis resamples these blocks 250 times.

Every resampled correlation is fitted on the selected window. The gap error is
the sample standard deviation of the successful bootstrap gaps. The excited
energy is calculated jointly in each resample:

$$
E_1^{(r)}
=
E_0^{(r)}+\Delta E^{(r)}.
$$

This preserves the covariance between the sampled energy and correlation
within each bootstrap replicate. All 250 fits succeeded for every production
parameter set.

## Independent Hamiltonian Benchmark

For every potential, the code discretises

$$
H
=
-\frac{1}{2m}\frac{d^2}{dx^2}+V(x)
$$

on a separate position grid. The central second derivative produces a real
symmetric tridiagonal matrix with diagonal elements

$$
H_{ii}
=
\frac{1}{m\Delta x^2}+V(x_i)
$$

and off-diagonal elements

$$
H_{i,i\pm1}
=
-\frac{1}{2m\Delta x^2}.
$$

`scipy.linalg.eigh_tridiagonal` returns the two lowest eigenvalues and the
ground-state eigenvector on a 4,000-point interior grid.

This solver provides a controlled numerical benchmark for the anharmonic and
double-well systems and reproduces the exact harmonic spectrum to the displayed
precision. It is methodologically independent of the PIMC path ensemble.

## Euclidean Instantons

For

$$
V(x)=ax^4+bx^2,
\qquad a>0,
\qquad b<0,
$$

the classical minima occur at

$$
x_{\min}
=
\sqrt{\frac{-b}{2a}}.
$$

The default well threshold is

$$
x_{\mathrm{threshold}}
=
\frac{1}{2}x_{\min}.
$$

Each time slice is assigned to the positive well when
$x\geq x_{\mathrm{threshold}}$, the negative well when
$x\leq-x_{\mathrm{threshold}}$, and the intervening barrier region otherwise.
After unresolved barrier slices are removed, changes between the two resolved
well labels are counted around the periodic path.

The reported diagnostic is

$$
r_{\mathrm{instanton}}
=
\frac{N_{\mathrm{transitions}}}{\beta}.
$$

This is an instanton density along imaginary time. It is not a directly
measured real-time transition frequency. Its quantitative value depends on the
threshold, lattice spacing, Euclidean extent, and treatment of rapid
recrossings. The independently fitted odd-state gap is the more direct spectral
measure of tunnelling-induced level splitting.

## Parameter Sets

The default production study uses:

| ID | Family | Parameters |
|---|---|---|
| H1 | Harmonic | $k=1$, $m=1$ |
| H2 | Harmonic | $k=2$, $m=1$ |
| H3 | Harmonic | $k=1$, $m=0.5$ |
| A1 | Anharmonic | $k=1$, $\lambda=0.1$, $m=1$ |
| A2 | Anharmonic | $k=1$, $\lambda=0.5$, $m=1$ |
| A3 | Anharmonic | $k=1$, $\lambda=0.1$, $m=0.5$ |
| D1 | Double well | $a=0.25$, $b=-1$, $m=1$ |
| D2 | Double well | $a=0.5$, $b=-1$, $m=1$ |
| D3 | Double well | $a=0.25$, $b=-1$, $m=0.5$ |

## Default Configuration

| Setting | Default | Description |
|---|---:|---|
| `beta` | `16.0` | Total Euclidean extent |
| `n_sites` | `256` | Periodic imaginary-time sites |
| `spacing` | `0.0625` | Derived lattice spacing $\beta/N$ |
| `thermalization_sweeps` | `2_000` | Discarded compound sweeps per chain |
| `measurements` | `2_000` | Saved configurations per chain |
| `sweeps_between` | `4` | Compound sweeps between measurements |
| `n_chains` | `4` | Independent streams per parameter set |
| `block_size` | `50` | Measurements per analysis block |
| `bootstrap_samples` | `250` | Block-bootstrap replicates |
| `target_acceptance` | `0.50` | Local proposal adaptation target |
| `seed` | `2025` | Reproducible root seed |

Custom studies can be run programmatically:

```python
import numpy as np

from quantum_oscillator_pimc import (
    SimulationConfig,
    anharmonic_potential,
    run_study,
)


config = SimulationConfig(
    beta=12.0,
    n_sites=192,
    thermalization_sweeps=2_000,
    measurements=2_000,
    n_chains=4,
    seed=1234,
)

potential = anharmonic_potential(k=1.0, quartic=0.2, mass=1.0)
seed = np.random.SeedSequence(config.seed)
result = run_study(potential, config, seed)

print(result.energy, result.energy_error)
print(result.excited_energy, result.excited_energy_error)
```

`n_sites` must be even and at least 32. The configuration also validates the
Euclidean extent, chain count, thermalisation, block count, bootstrap count,
measurement spacing, and target acceptance.

## Verified Default Results

### Ground and first excited states

| ID | PIMC $E_0$ | Reference $E_0$ | PIMC $E_1$ | Reference $E_1$ |
|---|---:|---:|---:|---:|
| H1 | $0.506324\pm0.005383$ | 0.500000 | $1.472180\pm0.026861$ | 1.499998 |
| H2 | $0.711434\pm0.005766$ | 0.707106 | $2.126104\pm0.033868$ | 2.121315 |
| H3 | $0.704913\pm0.005137$ | 0.707106 | $2.126178\pm0.032953$ | 2.121316 |
| A1 | $0.558154\pm0.005446$ | 0.559146 | $1.799281\pm0.029214$ | 1.769498 |
| A2 | $0.697162\pm0.006231$ | 0.696175 | $2.316026\pm0.034630$ | 2.324398 |
| A3 | $0.808701\pm0.005962$ | 0.817837 | $2.633387\pm0.035454$ | 2.617221 |
| D1 | $-0.312821\pm0.006372$ | -0.299522 | $0.043106\pm0.011407$ | 0.046370 |
| D2 | $0.058275\pm0.006209$ | 0.068893 | $0.854296\pm0.020010$ | 0.856512 |
| D3 | $-0.136181\pm0.007431$ | -0.130419 | $0.676137\pm0.022086$ | 0.661427 |

Every PIMC energy agrees with the corresponding independent benchmark within
three reported standard errors.

### First energy gaps

| ID | PIMC $\Delta E$ | Reference $\Delta E$ | Fit window | Diagonal $\chi^2/\mathrm{dof}$ |
|---|---:|---:|---:|---:|
| H1 | $0.965856\pm0.027696$ | 0.999998 | $[0.250,4.000]$ | 0.17 |
| H2 | $1.414670\pm0.034428$ | 1.414210 | $[0.250,4.000]$ | 0.14 |
| H3 | $1.421265\pm0.033927$ | 1.414210 | $[0.250,3.312]$ | 0.04 |
| A1 | $1.241127\pm0.029787$ | 1.210353 | $[0.250,3.562]$ | 0.14 |
| A2 | $1.618864\pm0.035102$ | 1.628223 | $[0.250,3.250]$ | 0.61 |
| A3 | $1.824687\pm0.036045$ | 1.799385 | $[0.250,2.625]$ | 0.08 |
| D1 | $0.355927\pm0.010928$ | 0.345892 | $[0.250,4.000]$ | 0.07 |
| D2 | $0.796020\pm0.019181$ | 0.787619 | $[0.250,4.000]$ | 0.15 |
| D3 | $0.812318\pm0.020487$ | 0.791847 | $[0.250,4.000]$ | 0.16 |

The displayed $\chi^2/\mathrm{dof}$ values use diagonal errors. Correlation
values at different imaginary-time separations are calculated from the same
paths and are therefore mutually correlated. These values are window
diagnostics, not complete covariance-aware goodness-of-fit probabilities.

### Sampling diagnostics

| ID | Local acceptance | Mean $\tau_{\mathrm{int},E}$ |
|---|---:|---:|
| H1 | 0.501 | 3.52 |
| H2 | 0.501 | 3.16 |
| H3 | 0.500 | 2.85 |
| A1 | 0.500 | 2.93 |
| A2 | 0.500 | 2.17 |
| A3 | 0.500 | 2.35 |
| D1 | 0.500 | 1.66 |
| D2 | 0.500 | 2.01 |
| D3 | 0.502 | 1.92 |

The 50-measurement analysis block is more than five times the largest mean
energy-autocorrelation time.

### Euclidean instanton densities

| ID | Transitions per unit Euclidean time |
|---|---:|
| D1 | $0.155000\pm0.003030$ |
| D2 | $0.327281\pm0.002860$ |
| D3 | $0.319453\pm0.002856$ |

The deeper, more widely separated D1 wells show the smallest thresholded
transition density and the smallest independently fitted odd-state gap.

## Convergence Checks

The default program performs two additional harmonic-oscillator studies with
exact ground-state energy $E_0=0.5$.

### Lattice-spacing stability

At fixed $\beta=8$:

| $a_\tau$ | PIMC $E_0$ |
|---:|---:|
| 0.2500 | $0.497204\pm0.003772$ |
| 0.1250 | $0.504065\pm0.004011$ |
| 0.0625 | $0.491876\pm0.005290$ |

All three estimates are statistically compatible with the exact result. The
central values need not vary monotonically because each point is an independent
finite Monte Carlo estimate.

### Euclidean-extent stability

At fixed $a_\tau=0.125$:

| $\beta$ | PIMC $E_0$ |
|---:|---:|
| 6 | $0.500104\pm0.004299$ |
| 10 | $0.498860\pm0.003727$ |
| 16 | $0.501295\pm0.003599$ |

The estimates remain stable as the Euclidean extent increases and are all
consistent with ground-state saturation.

## Terminal Output

For each parameter set, the program prints:

- the potential family and parameters;
- PIMC and reference $E_0$;
- PIMC and reference $E_1-E_0$;
- PIMC and reference $E_1$;
- binned or bootstrap uncertainties;
- the local acceptance rate;
- the integrated energy-autocorrelation time;
- the selected correlation-fit window;
- the diagonal reduced-$\chi^2$ diagnostic;
- the number of successful bootstrap fits;
- the Euclidean instanton density for double wells.

Separate tables report the lattice-spacing and Euclidean-extent checks.

## Plots

The program creates five Matplotlib figures.

### Energy comparison

Two panels compare PIMC and independently diagonalised $E_0$ and $E_1$ values
for all nine parameter sets. PIMC points include the reported uncertainties.

### Ground-state probability densities

A $3\times3$ grid compares the sampled time-slice histograms with the
finite-difference $|\psi_0(x)|^2$ curves.

### Correlation fits

Representative harmonic, anharmonic, and double-well correlations are plotted
on logarithmic axes. The automatically selected fit ranges and periodic
spectral curves are shown explicitly.

### Double-well tunnelling

The upper panel shows a representative periodic double-well path with the
positive and negative well thresholds. The lower panel compares the three
Euclidean instanton densities.

### Convergence diagnostics

Two panels show the harmonic ground-state energy under changes in lattice
spacing and total Euclidean extent.

Figures are displayed by default. They are saved only when `--save-dir` is
supplied, and generated files are ignored by Git.

## Returned Results

`run_default_analysis()` returns

```python
results, convergence = run_default_analysis(config)
```

`results` is a list of nine `StudyResult` objects. `convergence` is a dictionary
with `"lattice"` and `"extent"` entries.

Each `StudyResult` contains:

| Attribute | Description |
|---|---|
| `potential` | Validated `Potential` definition |
| `config` | Validated `SimulationConfig` |
| `energy`, `energy_error` | Virial $E_0$ estimate and blocked error |
| `gap`, `gap_error` | Correlation-fit $E_1-E_0$ and bootstrap error |
| `excited_energy`, `excited_energy_error` | Joint-bootstrap $E_1$ result |
| `acceptance` | Mean production local acceptance |
| `proposal_width` | Mean frozen local proposal width |
| `tau_energy` | Mean integrated energy-autocorrelation time |
| `reference` | Two eigenvalues, grid, and ground-state density |
| `density_samples` | Retained time-slice samples for histogramming |
| `correlation_time` | Imaginary-time separations |
| `correlation`, `correlation_error` | Mean periodic correlation and blocked errors |
| `fit` | Gap, uncertainty, window, diagnostic, and fitted curve |
| `energy_blocks` | Virial-energy block means |
| `correlation_blocks` | Complete correlation block means |
| `instanton_rate`, `instanton_rate_error` | Double-well Euclidean diagnostic |
| `representative_path` | One measured path for plotting |
| `instanton_threshold` | Applied double-well threshold |

Example:

```python
from quantum_oscillator_pimc import SimulationConfig, run_default_analysis


config = SimulationConfig(seed=2025)
results, convergence = run_default_analysis(config)

first_double_well = results[6]

print(first_double_well.energy, first_double_well.energy_error)
print(first_double_well.gap, first_double_well.gap_error)
print(first_double_well.instanton_rate)
print(convergence["lattice"])
```

## Main Functions

| Function or class | Description |
|---|---|
| `Potential` | Stores and validates one even polynomial potential |
| `SimulationConfig` | Stores and validates the lattice and analysis settings |
| `harmonic_potential()` | Constructs $kx^2/2$ |
| `anharmonic_potential()` | Constructs $kx^2/2+\lambda x^4$ |
| `double_well_potential()` | Constructs $ax^4+bx^2$ |
| `default_potentials()` | Returns the nine documented parameter sets |
| `local_action_change()` | Evaluates the local Euclidean action difference |
| `metropolis_sweep()` | Updates both checkerboard parities |
| `staging_move()` | Proposes a conditional free-particle Brownian bridge |
| `centroid_move()` | Translates the complete path |
| `compound_sweep()` | Combines all production update types |
| `circular_correlation()` | Calculates the periodic time-origin average by FFT |
| `integrated_autocorrelation_time()` | Estimates the positive-pair integrated time |
| `virial_energy()` | Evaluates the path virial estimator |
| `block_means()` | Blocks each chain independently |
| `reference_solution()` | Diagonalises the tridiagonal Hamiltonian |
| `select_fit_window()` | Selects a stable correlation-fit region |
| `bootstrap_fit()` | Bootstraps the gap and excited energy jointly |
| `count_instantons()` | Counts thresholded periodic well transitions |
| `run_study()` | Runs and analyses one potential |
| `run_convergence_study()` | Tests lattice spacing and Euclidean extent |
| `run_default_analysis()` | Runs all production and convergence ensembles |
| `create_figures()` | Creates the five diagnostic figures |
| `save_figures()` | Writes selected figures as PNG files |
| `print_results()` | Prints the production comparison table |
| `print_convergence()` | Prints the convergence tables |

## Repository Preparation and Verification

This repository develops the original coursework exploration into a
reproducible portfolio study. The central scientific correction is the
replacement of scalar Boltzmann sampling with complete periodic Euclidean
paths containing the required kinetic coupling.

The maintenance and analysis work includes:

- renaming the script to an importable filename;
- separating potentials, configuration, simulation, analysis, plotting, and
  command-line behaviour;
- implementing the discretised kinetic-plus-potential action;
- adding checkerboard, staging, centroid, and reflection moves;
- freezing adaptive proposals before measurement;
- adding independent deterministic random streams;
- replacing reshaped scalar samples with genuine time-lattice paths;
- adding the virial energy estimator;
- time-origin averaging the quantum correlation function;
- adding stable periodic fits and data-derived initial gaps;
- bootstrapping spectral results by complete measurement block;
- adding finite-difference Hamiltonian validation;
- replacing raw Markov-chain sign changes with thresholded Euclidean
  instanton counts;
- adding lattice-spacing and Euclidean-extent studies;
- declaring and validating the package dependencies.

Verification included:

- syntax and direct-entrypoint checks;
- exact deterministic reproduction under the default seed;
- validation failures for invalid masses, potentials, lattices, blocks, and
  acceptance targets;
- four independent production chains per parameter set;
- automatic rejection of inadequately short analysis blocks;
- local acceptance and energy-autocorrelation checks;
- 250 successful block-bootstrap fits per parameter set;
- PIMC agreement with every benchmark within three reported errors;
- harmonic analytical-spectrum agreement;
- lattice-spacing and Euclidean-extent stability;
- visual inspection of all five figures;
- a non-interactive fresh execution with exit code 0;
- credential, private-path, dependency, archive, and exclusion checks.

The original source and assignment brief remain unchanged and are not included
in this repository.

## Numerical Considerations

### Finite Euclidean extent

The ground-state interpretation is asymptotic in $\beta$. At insufficient
extent, the sampled density matrix includes thermally excited states. The
default $\beta=16$ and the separate extent scan support ground-state saturation
for the studied parameter sets, but custom potentials with smaller gaps may
require larger values.

### Primitive-action discretisation

The lattice action is a finite-$a_\tau$ primitive approximation. Results should
be checked as $a_\tau$ decreases while the physical Euclidean extent remains
fixed. Smaller spacing can increase autocorrelation and computational cost.

### Monte Carlo autocorrelation

Staging and centroid moves substantially reduce the long-wavelength slowing
seen with local updates alone. They do not make saved paths independent. The
block-size safeguard and reported $\tau_{\mathrm{int},E}$ diagnose only the
default energy sequence; custom observables can have different correlation
scales.

### Virial estimator

The virial estimator has favourable continuum behaviour for these smooth
potentials. It is not the only valid PIMC energy estimator. Thermodynamic and
combined estimators may provide useful cross-checks in extended work.

### Correlation fit window

Short imaginary times can contain contributions from higher odd states. Long
times become noise dominated. The automatic selector tests fit stability, but
the remaining window choice is a systematic consideration. The displayed
diagonal $\chi^2$ values do not incorporate the full covariance across
imaginary-time separations.

### Bootstrap interpretation

The bootstrap resamples measurement blocks and propagates their joint energy
and correlation fluctuations. It assumes the 50-measurement blocks are long
enough to behave approximately independently. The explicit block-to-$\tau$
safeguard supports this assumption for the default energy sequence.

### Finite-difference benchmark

The reference solver uses a finite spatial interval, finite grid spacing, and
Dirichlet boundary approximation. The interval is chosen so the displayed
ground-state density is negligible near its boundaries. More weakly confined
custom potentials should recheck box-size and grid convergence.

### Euclidean tunnelling diagnostic

The thresholded transition count identifies well-to-well structure along
imaginary time. It is sensitive to the chosen threshold and to lattice
resolution, and it should not be interpreted as a measured real-time rate.
The first odd-state gap provides an independent spectral diagnostic.

### Statistical agreement

Agreement within a chosen number of standard errors does not prove that every
systematic error is negligible. The independent Hamiltonian comparison,
continuum checks, extent checks, block tests, and visual density overlays are
used together rather than relying on one metric.

### Computational cost

The default study runs nine four-chain production ensembles and six additional
convergence ensembles. The `--quick` option is intended for installation and
workflow checks; its wider uncertainties should not replace the documented
production results.

## Expected Behaviour

A successful default run should show:

- local acceptance rates close to 0.50;
- energy autocorrelation times comfortably below the 50-measurement block;
- sampled densities following the independent $|\psi_0|^2$ curves;
- ground and first excited energies consistent with the benchmarks;
- harmonic energies close to the exact analytical spectrum;
- positive, stable first-gap fits in signal-dominated windows;
- all 250 bootstrap fits completing for every parameter set;
- a smaller D1 double-well gap and instanton density than D2 and D3;
- stable harmonic $E_0$ under the reported lattice and extent changes;
- identical output when the same seed and full configuration are reused.

## References

- Shikhar Mittal, Marise J. E. Westbroek, Peter R. King, and Dimitri D.
  Vvedensky, [*Path integral Monte Carlo method for the quantum anharmonic
  oscillator*](https://arxiv.org/abs/1811.04669).
- Wolfhard Janke and Tilman Sauer, [*Optimal Energy Estimation in Path-Integral
  Monte Carlo Simulations*](https://arxiv.org/abs/cond-mat/9710117).
- Colin Morningstar, [*Monte Carlo Methods in Quantum Mechanics*](https://arxiv.org/abs/hep-lat/0702020).

These references provide methodological context. The implementation and
documented results are produced by this repository.

## Coursework Context

The original program was written by Jack Turner for PH-353 Computational
Physics coursework at Swansea University in 2025. This repository presents the
oscillator study as an academic and portfolio project with a corrected quantum
path ensemble, additional statistical analysis, independent validation,
convergence diagnostics, and documentation.

The lecturer-provided assignment sheet and separate poster submission are not
distributed with the repository.

## Licence

No open-source licence has currently been applied.

Copyright remains with Jack Turner. The absence of a licence means that
permission to copy, modify, or redistribute the code should not be assumed.
