# Improvement Record: Quantum-Oscillator Path-Integral Monte Carlo — Gate 3

## Archive represented

This document describes the numerical rebuild contained in:

```text
Quantum-Oscillator-Path-Integral-Monte-Carlo-Gate3.zip
```

The accompanying `353 CA3 Notes.py` is the untouched original source used as the baseline. Its SHA-256 checksum is:

```text
51c19f633293ccf9a13a303dc7d3f98eb7c65b2946e20cece59a88266d3be609  353 CA3 Notes.py
```

Gate 3 contains the verified numerical rebuild from Gate 2 and adds its publication README. The simulation code and reference results were not changed during the documentation gate.

## Executive summary

The original script explored harmonic, anharmonic, and double-well potentials using a scalar Metropolis sampler. Although several functions referred to paths and excited states, the core sampler updated one coordinate according to a classical weight proportional to

$$
e^{-\beta V(x)}.
$$

That distribution does not include the kinetic coupling required for a quantum Euclidean path integral. It therefore could not, by itself, provide a controlled calculation of quantum ground-state densities, energies, spectral gaps, or tunnelling paths.

The rebuilt project implements an actual periodic Euclidean-time lattice. It samples complete paths using local and non-local Monte Carlo moves, estimates correlated uncertainties, extracts spectral gaps from periodic correlation functions, compares every parameter set with an independent finite-difference Hamiltonian, and performs lattice-spacing and Euclidean-extent checks.

This was a scientific-method rebuild, not only a style refactor.

## Original implementation limitations

- The main Metropolis routine sampled one scalar coordinate rather than a periodic path.
- The sampling weight omitted the kinetic term coupling neighbouring imaginary-time sites.
- The scalar density $e^{-\beta V(x)}$ was treated as a theoretical comparison even though it is a classical Boltzmann density, not generally $|\psi_0(x)|^2$.
- Ground- and excited-state energy functions were not derived from a controlled lattice estimator and spectral fit.
- Tunnelling was inferred without a well-defined periodic-path ensemble and threshold prescription.
- One implicit random stream was used without deterministic independent chains.
- The code had no burn-in adaptation policy, autocorrelation-aware block design, bootstrap spectral uncertainty, continuum check, ground-state projection check, or independent Schrödinger benchmark.
- Plotting and numerical work were coupled in top-level execution.
- There was no CLI, quick mode, figure-saving mode, validated configuration, or dependency declaration.

## Improvement summary

| Area | Improvement | Why it was made | Impact |
|---|---|---|---|
| Quantum formulation | Replaced scalar sampling with a periodic Euclidean-time lattice | Quantum mechanics requires the kinetic coupling between neighbouring time slices | The simulation now samples a legitimate discretised path integral |
| Local updates | Added even/odd checkerboard Metropolis sweeps | Neighbouring path sites are coupled | Independent parity updates allow efficient vectorised local proposals |
| Non-local updates | Added free-particle staging, centroid, and reflection moves | Local updates mix slowly on fine lattices and in double wells | Long-wavelength and between-well modes are sampled more effectively |
| Adaptation | Tuned local proposal width during early thermalisation and then froze it | Adapting during measurement would alter the stationary chain | Production acceptance is stable near 0.50 without compromising the target distribution |
| Random streams | Added deterministic independent generators for every chain | One implicit stream weakens reproducibility and multi-chain diagnostics | Default results reproduce exactly while chains remain independent |
| Energy estimator | Added a virial energy estimator appropriate to the lattice ensemble | The original energy calculation was not a controlled quantum estimator | Ground-state energies can be compared with independent Hamiltonian eigenvalues |
| Correlation analysis | Added periodic time-origin-averaged correlation functions | Excited-state gaps are encoded in imaginary-time decay | $E_1-E_0$ is extracted from the sampled paths rather than assumed |
| Fit selection | Added deterministic fit-window selection and periodic spectral models | Arbitrary windows can bias gap estimates | The chosen window and fit diagnostic are reported for every ensemble |
| Uncertainty | Added blocking, integrated autocorrelation times, and block bootstrap fits | Saved configurations and correlation points are correlated | Reported errors reflect the chain structure and nonlinear fit uncertainty |
| Benchmark | Added finite-difference Hamiltonian diagonalisation | Anharmonic and double-well spectra lack simple general closed forms | Every Monte Carlo energy has an independent numerical reference |
| Tunnelling | Added thresholded Euclidean instanton counting | Raw zero crossings are dominated by small fluctuations | Double-well transition densities have a stated operational definition |
| Convergence | Added lattice-spacing and Euclidean-extent studies | One lattice cannot establish discretisation or ground-state projection stability | The default harmonic result is checked against changes in $a_\tau$ and $\beta$ |
| Packaging | Added validated data classes, CLI modes, save support, requirements, and `.gitignore` | The original notes were not repository-ready | The study can be reproduced and reviewed as a standalone project |

## Euclidean path-integral correction

For a one-dimensional Hamiltonian

$$
H=\frac{p^2}{2m}+V(x),
$$

the Euclidean partition function is

$$
Z(\beta)
=
\int_{x(\beta)=x(0)}\mathcal D x\;e^{-S_E[x]},
$$

with action

$$
S_E[x]
=
\int_0^\beta
\left[
\frac{m}{2}\left(\frac{dx}{d\tau}\right)^2
+V(x)
\right]d\tau.
$$

The rebuild discretises $[0,\beta)$ into $N$ sites with

$$
a_\tau=\frac{\beta}{N}
$$

and periodic indexing $x_N=x_0$. The primitive lattice action is

$$
S_E^{(N)}
=
\sum_{j=0}^{N-1}
\left[
\frac{m}{2a_\tau}(x_{j+1}-x_j)^2
+a_\tau V(x_j)
\right].
$$

The kinetic term

$$
\frac{m}{2a_\tau}(x_{j+1}-x_j)^2
$$

is the decisive difference from the original scalar sampler. It correlates adjacent time slices and gives the ensemble its quantum character.

## Monte Carlo improvements

### Local checkerboard moves

For site $j$, the program proposes

$$
x_j'=x_j+\delta,
\qquad
\delta\sim\mathcal U(-w,w).
$$

The proposal is accepted with

$$
P_{\mathrm{accept}}
=
\min\left(1,e^{-\Delta S_E}\right).
$$

Even and odd sites are updated in separate vectorised passes. Sites with the same parity are not directly linked, so their local changes can be evaluated together.

### Thermalisation-only adaptation

The initial width scales as

$$
w_0=1.8\sqrt{\frac{a_\tau}{m}}.
$$

It is adjusted toward a 50% local acceptance rate during early thermalisation and then frozen before measurements. The impact is efficient sampling without an adaptive production kernel.

### Staging moves

Fine lattices strongly couple neighbouring sites. The rebuild proposes complete segments using free-particle Brownian bridges. Because the proposal samples the conditional kinetic action exactly, the acceptance ratio depends only on the potential change:

$$
P_{\mathrm{stage}}
=
\min\left[
1,
\exp\left(
-a_\tau\sum_{j\in\mathrm{segment}}
[V(x_j')-V(x_j)]
\right)
\right].
$$

This reduces the random-walk behaviour of purely local updates.

### Centroid and reflection moves

Centroid shifts move the entire path and improve exploration of long-wavelength modes. Reflection proposals help an even-potential path move between symmetry-related wells. These updates are especially important for double-well ensembles.

## Statistical improvements

### Autocorrelation-aware blocking

The integrated autocorrelation time is estimated from the energy series. Production measurements are divided into blocks longer than the observed correlation scale. If $\bar E_b$ is the mean of block $b$, the error is based on the distribution of block means rather than treating every saved path as independent.

The default block length of 50 measurements exceeds five times the largest reported mean energy autocorrelation time.

### Spectral-gap extraction

For an even ground state and the position operator, the leading periodic correlation model is

$$
C(\tau)
=
A\left[
e^{-\Delta E\tau}
+e^{-\Delta E(\beta-\tau)}
\right],
$$

where

$$
\Delta E=E_1-E_0.
$$

The program averages over all time origins, selects a reproducible fit window, fits the periodic model, and bootstraps blocks of paths. The first excited state is then

$$
E_1=E_0+\Delta E.
$$

This replaces the original heuristic excited-state calculation with an observable tied to the sampled Euclidean ensemble.

### Independent Hamiltonian benchmark

The Schrödinger operator

$$
-\frac{1}{2m}\frac{d^2}{dx^2}+V(x)
$$

is discretised on a separate spatial grid and diagonalised. These eigenvalues do not control the Monte Carlo simulation or fit window; they are used only for validation.

## Double-well tunnelling improvement

For

$$
V(x)=ax^4+bx^2,
\qquad
a>0,
\qquad
b<0,
$$

the minima lie at

$$
x_\mathrm{min}=\pm\sqrt{\frac{-b}{2a}}.
$$

The rebuild defines positive and negative well thresholds and counts transitions only after a path reaches the opposite threshold. This suppresses false instantons caused by rapid recrossings near $x=0$.

The reported quantity is a Euclidean transition density, not a real-time tunnelling rate.

## Verification and measured impact

The production suite contains nine parameter sets: three harmonic, three quartic anharmonic, and three double-well systems.

Every reported PIMC ground- and first-excited-state energy agreed with the independent finite-difference benchmark within three reported standard errors.

For the standard harmonic oscillator, the rebuild obtained

$$
E_0=0.506324\pm0.005383
$$

against the exact value

$$
E_0=0.5,
$$

and

$$
E_1=1.472180\pm0.026861
$$

against the reference value approximately $1.5$.

The separate convergence studies found harmonic ground-state estimates statistically compatible with $0.5$ across the tested lattice spacings and Euclidean extents.

The three double-well ensembles produced explicit instanton-density estimates. The deepest, most widely separated wells had the smallest transition density and the smallest independently fitted odd-state gap, matching the expected qualitative relationship.

## Gate 3 repository impact

The Gate 3 archive contains:

- `quantum_oscillator_pimc.py`, containing the complete rebuilt simulation;
- compatible NumPy, SciPy, and Matplotlib requirements;
- a `.gitignore` excluding caches and generated figures;
- `--quick`, `--no-show`, and `--save-dir` execution modes;
- structured return values and deterministic configuration classes;
- a publication-quality `README.md` explaining the physics, equations, algorithms, validation evidence, installation, commands, API, outputs, limitations, and reproducibility controls.

### Documentation impact

The Gate 3 README turns the verified implementation into a repository that can be evaluated without first reverse-engineering the source. It distinguishes the original classical scalar sampler from the rebuilt quantum path integral, derives the lattice quantities with properly rendered LaTeX, reports the verified energy and gap results, and explains how to reproduce quick or full runs. It also exposes the discretisation, fit-covariance, instanton-definition, and computational-cost limitations next to the claims they qualify. This gate improved clarity and publication readiness without changing the numerical method or its measured outputs.

## Remaining limitations

- The primitive lattice action has finite-$a_\tau$ discretisation error.
- Spectral fits use diagonal correlation errors rather than a full covariance matrix.
- Finite-difference benchmarks have finite spatial range and spacing.
- Instanton density depends on the chosen threshold definition.
- Full production runs are computationally more expensive than the original scalar sampler.
- No open-source licence was applied to the produced archive.
