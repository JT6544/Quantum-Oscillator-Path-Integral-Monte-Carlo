"""Euclidean path-integral Monte Carlo for one-dimensional quantum oscillators.

The program studies harmonic, quartic anharmonic, and symmetric double-well
potentials in units where hbar = 1. It samples periodic imaginary-time paths,
estimates the ground-state energy and density, extracts the first energy gap
from time-displaced correlation functions, and diagnoses Euclidean instantons
in double-well configurations.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import curve_fit


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Potential:
    """Even polynomial potential ``V(x) = quadratic*x^2 + quartic*x^4``."""

    label: str
    family: str
    mass: float
    quadratic: float
    quartic: float = 0.0

    def __post_init__(self) -> None:
        if self.mass <= 0:
            raise ValueError("mass must be positive")
        if self.quartic < 0:
            raise ValueError("quartic coefficient must be non-negative")
        if self.quartic == 0 and self.quadratic <= 0:
            raise ValueError("a purely quadratic potential must be confining")
        if self.family not in {"harmonic", "anharmonic", "double_well"}:
            raise ValueError(f"unknown potential family: {self.family}")

    def value(self, x: FloatArray) -> FloatArray:
        return self.quadratic * x**2 + self.quartic * x**4

    def derivative(self, x: FloatArray) -> FloatArray:
        return 2.0 * self.quadratic * x + 4.0 * self.quartic * x**3

    @property
    def minimum_position(self) -> float:
        if self.family != "double_well":
            return 0.0
        return float(np.sqrt(-self.quadratic / (2.0 * self.quartic)))

    @property
    def parameter_text(self) -> str:
        if self.family == "harmonic":
            return f"k={2.0 * self.quadratic:g}, m={self.mass:g}"
        if self.family == "anharmonic":
            return (
                f"k={2.0 * self.quadratic:g}, "
                f"lambda={self.quartic:g}, m={self.mass:g}"
            )
        return (
            f"a={self.quartic:g}, b={self.quadratic:g}, "
            f"m={self.mass:g}"
        )


@dataclass(frozen=True)
class SimulationConfig:
    """Configuration for a path-integral Monte Carlo ensemble."""

    beta: float = 16.0
    n_sites: int = 256
    thermalization_sweeps: int = 2_000
    measurements: int = 2_000
    sweeps_between: int = 4
    n_chains: int = 4
    block_size: int = 50
    bootstrap_samples: int = 250
    target_acceptance: float = 0.50
    seed: int = 2025

    def __post_init__(self) -> None:
        if self.beta <= 0:
            raise ValueError("beta must be positive")
        if self.n_sites < 32 or self.n_sites % 2:
            raise ValueError("n_sites must be an even integer of at least 32")
        if self.thermalization_sweeps < 100:
            raise ValueError("thermalization_sweeps must be at least 100")
        if self.measurements < 2 * self.block_size:
            raise ValueError("measurements must contain at least two blocks")
        if self.sweeps_between < 1:
            raise ValueError("sweeps_between must be positive")
        if self.n_chains < 2:
            raise ValueError("at least two independent chains are required")
        if self.block_size < 2:
            raise ValueError("block_size must be at least 2")
        if self.bootstrap_samples < 50:
            raise ValueError("bootstrap_samples must be at least 50")
        if not 0.2 <= self.target_acceptance <= 0.8:
            raise ValueError("target_acceptance must lie between 0.2 and 0.8")

    @property
    def spacing(self) -> float:
        return self.beta / self.n_sites


@dataclass(frozen=True)
class ReferenceSolution:
    energies: FloatArray
    x: FloatArray
    ground_density: FloatArray


@dataclass(frozen=True)
class FitResult:
    gap: float
    gap_error: float
    fit_start: float
    fit_end: float
    reduced_chi_squared: float
    successful_bootstraps: int
    curve: FloatArray


@dataclass
class StudyResult:
    potential: Potential
    config: SimulationConfig
    energy: float
    energy_error: float
    gap: float
    gap_error: float
    excited_energy: float
    excited_energy_error: float
    acceptance: float
    proposal_width: float
    tau_energy: float
    reference: ReferenceSolution
    density_samples: FloatArray
    correlation_time: FloatArray
    correlation: FloatArray
    correlation_error: FloatArray
    fit: FitResult
    energy_blocks: FloatArray
    correlation_blocks: FloatArray
    instanton_rate: float | None = None
    instanton_rate_error: float | None = None
    representative_path: FloatArray | None = None
    instanton_threshold: float | None = None


def harmonic_potential(k: float, mass: float) -> Potential:
    return Potential(
        label="Harmonic oscillator",
        family="harmonic",
        mass=mass,
        quadratic=0.5 * k,
    )


def anharmonic_potential(k: float, quartic: float, mass: float) -> Potential:
    return Potential(
        label="Anharmonic oscillator",
        family="anharmonic",
        mass=mass,
        quadratic=0.5 * k,
        quartic=quartic,
    )


def double_well_potential(a: float, b: float, mass: float) -> Potential:
    if a <= 0 or b >= 0:
        raise ValueError("a double well requires a > 0 and b < 0")
    return Potential(
        label="Double-well oscillator",
        family="double_well",
        mass=mass,
        quadratic=b,
        quartic=a,
    )


def default_potentials() -> list[Potential]:
    """Return the nine parameter sets used by the default study."""

    return [
        harmonic_potential(k=1.0, mass=1.0),
        harmonic_potential(k=2.0, mass=1.0),
        harmonic_potential(k=1.0, mass=0.5),
        anharmonic_potential(k=1.0, quartic=0.1, mass=1.0),
        anharmonic_potential(k=1.0, quartic=0.5, mass=1.0),
        anharmonic_potential(k=1.0, quartic=0.1, mass=0.5),
        double_well_potential(a=0.25, b=-1.0, mass=1.0),
        double_well_potential(a=0.5, b=-1.0, mass=1.0),
        double_well_potential(a=0.25, b=-1.0, mass=0.5),
    ]


def local_action_change(
    old: FloatArray,
    proposed: FloatArray,
    left: FloatArray,
    right: FloatArray,
    potential: Potential,
    spacing: float,
) -> FloatArray:
    """Return the action change for independent same-parity site proposals."""

    kinetic_scale = potential.mass / (2.0 * spacing)
    old_kinetic = (old - left) ** 2 + (right - old) ** 2
    new_kinetic = (proposed - left) ** 2 + (right - proposed) ** 2
    potential_change = spacing * (
        potential.value(proposed) - potential.value(old)
    )
    return kinetic_scale * (new_kinetic - old_kinetic) + potential_change


def metropolis_sweep(
    path: FloatArray,
    potential: Potential,
    spacing: float,
    proposal_width: float,
    rng: np.random.Generator,
) -> float:
    """Perform one checkerboard Metropolis sweep and return acceptance."""

    accepted = 0
    n_sites = path.size
    for parity in (0, 1):
        indices = np.arange(parity, n_sites, 2)
        old = path[indices].copy()
        proposed = old + rng.uniform(-proposal_width, proposal_width, old.size)
        left = path[(indices - 1) % n_sites]
        right = path[(indices + 1) % n_sites]
        delta_action = local_action_change(
            old, proposed, left, right, potential, spacing
        )
        accept = np.log(rng.random(old.size)) < -delta_action
        path[indices[accept]] = proposed[accept]
        accepted += int(np.count_nonzero(accept))
    return accepted / n_sites


def staging_move(
    path: FloatArray,
    potential: Potential,
    spacing: float,
    segment_length: int,
    rng: np.random.Generator,
) -> bool:
    """Propose a free-particle Brownian bridge with fixed endpoints.

    The proposal is drawn exactly from the kinetic part of the conditional
    path distribution, so the Metropolis decision contains only the change
    in the potential action.
    """

    n_sites = path.size
    segment_length = min(segment_length, n_sites - 2)
    if segment_length < 2:
        return False
    start = int(rng.integers(0, n_sites))
    indices = (start + np.arange(segment_length + 1)) % n_sites
    left_endpoint = path[indices[0]]
    right_endpoint = path[indices[-1]]
    increments = rng.normal(
        scale=np.sqrt(spacing / potential.mass), size=segment_length
    )
    increments -= np.mean(increments)
    increments += (right_endpoint - left_endpoint) / segment_length
    proposed_segment = left_endpoint + np.cumsum(increments)
    proposed_interior = proposed_segment[:-1]
    interior_indices = indices[1:-1]
    old_interior = path[interior_indices]
    potential_change = spacing * np.sum(
        potential.value(proposed_interior) - potential.value(old_interior)
    )
    if np.log(rng.random()) < -potential_change:
        path[interior_indices] = proposed_interior
        return True
    return False


def centroid_move(
    path: FloatArray,
    potential: Potential,
    spacing: float,
    width: float,
    rng: np.random.Generator,
) -> bool:
    """Translate the whole path; the periodic kinetic action is unchanged."""

    shift = float(rng.uniform(-width, width))
    proposed = path + shift
    potential_change = spacing * np.sum(
        potential.value(proposed) - potential.value(path)
    )
    if np.log(rng.random()) < -potential_change:
        path[:] = proposed
        return True
    return False


def compound_sweep(
    path: FloatArray,
    potential: Potential,
    config: SimulationConfig,
    proposal_width: float,
    rng: np.random.Generator,
) -> float:
    """Combine local, staging, centroid, and symmetry-preserving moves."""

    local_acceptance = metropolis_sweep(
        path, potential, config.spacing, proposal_width, rng
    )
    short_segment = min(16, config.n_sites // 4)
    long_segment = min(32, config.n_sites // 2)
    staging_move(path, potential, config.spacing, short_segment, rng)
    staging_move(path, potential, config.spacing, long_segment, rng)
    centroid_width = 0.8 / np.sqrt(config.beta)
    centroid_move(path, potential, config.spacing, centroid_width, rng)
    if rng.random() < 0.05:
        path *= -1.0
    return local_acceptance


def circular_correlation(path: FloatArray) -> FloatArray:
    """Time-origin-averaged periodic correlation for one path."""

    transform = np.fft.rfft(path)
    correlation = np.fft.irfft(transform * transform.conjugate(), n=path.size)
    return np.asarray(correlation[: path.size // 2 + 1].real / path.size)


def integrated_autocorrelation_time(values: FloatArray) -> float:
    """Estimate integrated autocorrelation time using positive pair sums."""

    centered = np.asarray(values, dtype=float) - np.mean(values)
    if centered.size < 4 or np.allclose(centered, 0.0):
        return 0.5
    n = centered.size
    padded = 1 << (2 * n - 1).bit_length()
    transform = np.fft.rfft(centered, n=padded)
    autocovariance = np.fft.irfft(transform * transform.conjugate(), n=padded)[:n]
    autocovariance /= np.arange(n, 0, -1)
    if autocovariance[0] <= 0:
        return 0.5
    autocorrelation = autocovariance / autocovariance[0]
    tau = 0.5
    for start in range(1, n - 1, 2):
        pair = autocorrelation[start] + autocorrelation[start + 1]
        if pair <= 0:
            break
        tau += float(pair)
    return max(0.5, tau)


def count_instantons(path: FloatArray, threshold: float) -> int:
    """Count thresholded well-to-well transitions on one periodic path."""

    states = np.zeros(path.size, dtype=np.int8)
    states[path >= threshold] = 1
    states[path <= -threshold] = -1
    resolved = states[states != 0]
    if resolved.size < 2:
        return 0
    changes = int(np.count_nonzero(resolved[1:] != resolved[:-1]))
    changes += int(resolved[-1] != resolved[0])
    return changes


def initial_path(
    potential: Potential,
    n_sites: int,
    chain_index: int,
    rng: np.random.Generator,
) -> FloatArray:
    """Construct a dispersed but finite starting path."""

    if potential.family == "double_well":
        sign = -1.0 if chain_index % 2 else 1.0
        center = sign * potential.minimum_position
        return center + rng.normal(scale=0.15, size=n_sites)
    omega = np.sqrt(2.0 * potential.quadratic / potential.mass)
    scale = np.sqrt(1.0 / (2.0 * potential.mass * omega))
    return rng.normal(scale=scale, size=n_sites)


def virial_energy(path: FloatArray, potential: Potential) -> float:
    """Compute the path-averaged virial energy estimator."""

    values = potential.value(path) + 0.5 * path * potential.derivative(path)
    return float(np.mean(values))


def block_means(values: FloatArray, block_size: int) -> FloatArray:
    """Block the measurement axis independently for every chain."""

    n_chains, n_measurements = values.shape[:2]
    n_blocks = n_measurements // block_size
    if n_blocks < 2:
        raise ValueError("too few complete blocks")
    trimmed = values[:, : n_blocks * block_size]
    new_shape = (n_chains, n_blocks, block_size) + values.shape[2:]
    blocked = trimmed.reshape(new_shape).mean(axis=2)
    return blocked.reshape((n_chains * n_blocks,) + values.shape[2:])


def run_chain(
    potential: Potential,
    config: SimulationConfig,
    chain_index: int,
    seed: np.random.SeedSequence,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, float, float, FloatArray]:
    """Run one independently seeded path chain."""

    rng = np.random.default_rng(seed)
    path = initial_path(potential, config.n_sites, chain_index, rng)
    spacing = config.spacing
    proposal_width = 1.8 * np.sqrt(spacing / potential.mass)

    adaptation_interval = 100
    interval_acceptance = []
    frozen_thermalization = max(adaptation_interval, config.thermalization_sweeps // 5)
    adaptation_end = config.thermalization_sweeps - frozen_thermalization
    for sweep in range(config.thermalization_sweeps):
        interval_acceptance.append(
            compound_sweep(path, potential, config, proposal_width, rng)
        )
        if (
            sweep + 1 <= adaptation_end
            and (sweep + 1) % adaptation_interval == 0
        ):
            observed = float(np.mean(interval_acceptance))
            proposal_width *= np.exp(0.8 * (observed - config.target_acceptance))
            proposal_width = float(
                np.clip(
                    proposal_width,
                    0.2 * np.sqrt(spacing / potential.mass),
                    5.0 * np.sqrt(spacing / potential.mass),
                )
            )
            interval_acceptance.clear()

    n_lags = config.n_sites // 2 + 1
    energies = np.empty(config.measurements)
    correlations = np.empty((config.measurements, n_lags))
    density_stride = max(1, config.n_sites // 64)
    density_samples = np.empty(
        (config.measurements, (config.n_sites + density_stride - 1) // density_stride)
    )
    instantons = np.full(config.measurements, np.nan)
    threshold = 0.5 * potential.minimum_position
    acceptance_values = []
    representative_path = path.copy()

    for measurement in range(config.measurements):
        for _ in range(config.sweeps_between):
            acceptance_values.append(
                compound_sweep(path, potential, config, proposal_width, rng)
            )
        energies[measurement] = virial_energy(path, potential)
        correlations[measurement] = circular_correlation(path)
        selected = path[::density_stride]
        density_samples[measurement, : selected.size] = selected
        if potential.family == "double_well":
            instantons[measurement] = count_instantons(path, threshold)
        if measurement == config.measurements // 2:
            representative_path = path.copy()

    return (
        energies,
        correlations,
        density_samples.ravel(),
        instantons,
        float(np.mean(acceptance_values)),
        proposal_width,
        representative_path,
    )


def reference_solution(
    potential: Potential,
    n_grid: int = 4_000,
) -> ReferenceSolution:
    """Diagonalise a finite-difference Hamiltonian for validation."""

    if potential.family == "double_well":
        scale = max(6.0, 3.5 * potential.minimum_position)
    else:
        scale = 8.0 if potential.mass >= 0.75 else 10.0
    x = np.linspace(-scale, scale, n_grid + 2)[1:-1]
    dx = x[1] - x[0]
    diagonal = 1.0 / (potential.mass * dx**2) + potential.value(x)
    off_diagonal = np.full(n_grid - 1, -0.5 / (potential.mass * dx**2))
    energies, eigenvectors = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        select="i",
        select_range=(0, 1),
        check_finite=False,
    )
    ground_density = eigenvectors[:, 0] ** 2 / dx
    return ReferenceSolution(
        energies=np.asarray(energies),
        x=x,
        ground_density=ground_density,
    )


def periodic_correlation_model(
    time: FloatArray,
    amplitude: float,
    gap: float,
    beta: float,
) -> FloatArray:
    return amplitude * (np.exp(-gap * time) + np.exp(-gap * (beta - time)))


def fit_correlation_window(
    time: FloatArray,
    correlation: FloatArray,
    error: FloatArray,
    beta: float,
    start_index: int,
    end_index: int,
    initial_gap: float,
) -> tuple[FloatArray, FloatArray, float]:
    selected_time = time[start_index : end_index + 1]
    selected_correlation = correlation[start_index : end_index + 1]
    selected_error = np.maximum(
        error[start_index : end_index + 1],
        np.max(error[start_index : end_index + 1]) * 1e-3,
    )

    def model(t: FloatArray, amplitude: float, gap: float) -> FloatArray:
        return periodic_correlation_model(t, amplitude, gap, beta)

    parameters, covariance = curve_fit(
        model,
        selected_time,
        selected_correlation,
        p0=(max(correlation[0], 1e-8), initial_gap),
        sigma=selected_error,
        absolute_sigma=True,
        bounds=([0.0, 0.02], [np.inf, 10.0]),
        maxfev=20_000,
    )
    residual = (selected_correlation - model(selected_time, *parameters)) / selected_error
    degrees_of_freedom = max(1, selected_time.size - parameters.size)
    reduced_chi_squared = float(np.sum(residual**2) / degrees_of_freedom)
    return parameters, covariance, reduced_chi_squared


def select_fit_window(
    time: FloatArray,
    correlation: FloatArray,
    error: FloatArray,
    beta: float,
    initial_gap: float,
) -> tuple[int, int, FloatArray, FloatArray, float]:
    """Select the earliest statistically acceptable stable cosh window."""

    positive = np.flatnonzero(
        (correlation > 0.0) & (correlation > 2.0 * np.maximum(error, 1e-15))
    )
    maximum_end = min(len(time) - 2, int(0.5 * len(time)))
    if positive.size:
        maximum_end = min(maximum_end, int(positive[-1]))
    minimum_points = 8
    minimum_end = max(minimum_points + 2, int(np.ceil(1.5 / (time[1] - time[0]))))
    end_index = max(minimum_end, maximum_end)
    end_index = min(end_index, len(time) - 2)

    first_start = max(1, int(np.ceil(0.25 / (time[1] - time[0]))))
    last_start = min(
        int(np.ceil(1.5 / (time[1] - time[0]))),
        end_index - minimum_points + 1,
    )
    candidates = []
    for start_index in range(first_start, last_start + 1, 2):
        try:
            parameters, covariance, reduced_chi_squared = fit_correlation_window(
                time,
                correlation,
                error,
                beta,
                start_index,
                end_index,
                initial_gap,
            )
        except (RuntimeError, ValueError, FloatingPointError):
            continue
        gap_error = float(np.sqrt(max(covariance[1, 1], 0.0)))
        candidates.append(
            (
                start_index,
                parameters,
                covariance,
                reduced_chi_squared,
                gap_error,
            )
        )
    if not candidates:
        raise RuntimeError("no valid correlation fit window was found")

    for index, candidate in enumerate(candidates[:-2]):
        _, parameters, _, reduced_chi_squared, gap_error = candidate
        later = candidates[index + 1 : index + 3]
        stable = all(
            abs(parameters[1] - item[1][1])
            <= 2.0 * np.hypot(gap_error, item[4])
            for item in later
        )
        if stable and reduced_chi_squared <= 2.5:
            return (
                candidate[0],
                end_index,
                candidate[1],
                candidate[2],
                candidate[3],
            )

    best = min(
        candidates,
        key=lambda item: abs(np.log(max(item[3], 1e-8))) + 0.02 * time[item[0]],
    )
    return best[0], end_index, best[1], best[2], best[3]


def correlation_initial_gap(time: FloatArray, correlation: FloatArray) -> float:
    """Estimate a fit starting value from early positive log ratios."""

    spacing = time[1] - time[0]
    last_index = min(len(correlation) - 1, max(4, int(np.ceil(1.0 / spacing))))
    ratios = correlation[1 : last_index + 1] / correlation[:last_index]
    valid = np.isfinite(ratios) & (ratios > 0.0) & (ratios < 1.0)
    if not np.any(valid):
        return 1.0
    estimates = -np.log(ratios[valid]) / spacing
    return float(np.clip(np.median(estimates), 0.05, 5.0))


def bootstrap_fit(
    energy_blocks: FloatArray,
    correlation_blocks: FloatArray,
    time: FloatArray,
    beta: float,
    start_index: int,
    end_index: int,
    initial_gap: float,
    samples: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    """Bootstrap the joint excited-state energy and gap uncertainties."""

    n_blocks = energy_blocks.size
    gaps = []
    excited_energies = []
    correlation_scale = np.std(correlation_blocks, axis=0, ddof=1)
    correlation_scale /= np.sqrt(n_blocks)
    for _ in range(samples):
        indices = rng.integers(0, n_blocks, n_blocks)
        boot_correlation = np.mean(correlation_blocks[indices], axis=0)
        try:
            parameters, _, _ = fit_correlation_window(
                time,
                boot_correlation,
                correlation_scale,
                beta,
                start_index,
                end_index,
                initial_gap,
            )
        except (RuntimeError, ValueError, FloatingPointError):
            continue
        gap = float(parameters[1])
        gaps.append(gap)
        excited_energies.append(float(np.mean(energy_blocks[indices])) + gap)
    required = max(40, int(0.8 * samples))
    if len(gaps) < required:
        raise RuntimeError(
            f"only {len(gaps)} of {samples} bootstrap correlation fits succeeded"
        )
    return (
        float(np.std(gaps, ddof=1)),
        float(np.std(excited_energies, ddof=1)),
        len(gaps),
    )


def analyse_ensemble(
    potential: Potential,
    config: SimulationConfig,
    chain_data: list[
        tuple[FloatArray, FloatArray, FloatArray, FloatArray, float, float, FloatArray]
    ],
    bootstrap_seed: np.random.SeedSequence,
) -> StudyResult:
    energies = np.stack([item[0] for item in chain_data])
    correlations = np.stack([item[1] for item in chain_data])
    density_samples = np.concatenate([item[2] for item in chain_data])
    instantons = np.stack([item[3] for item in chain_data])
    acceptance = float(np.mean([item[4] for item in chain_data]))
    proposal_width = float(np.mean([item[5] for item in chain_data]))

    energy_blocks = block_means(energies, config.block_size)
    correlation_blocks = block_means(correlations, config.block_size)
    energy = float(np.mean(energy_blocks))
    energy_error = float(np.std(energy_blocks, ddof=1) / np.sqrt(energy_blocks.size))
    tau_energy = float(
        np.mean([integrated_autocorrelation_time(chain) for chain in energies])
    )
    if config.block_size < 5.0 * tau_energy:
        raise RuntimeError(
            f"block size {config.block_size} is too short for tau={tau_energy:.2f}"
        )

    correlation = np.mean(correlation_blocks, axis=0)
    correlation_error = np.std(correlation_blocks, axis=0, ddof=1)
    correlation_error /= np.sqrt(correlation_blocks.shape[0])
    time = np.arange(correlation.size) * config.spacing

    reference = reference_solution(potential)
    initial_gap = correlation_initial_gap(time, correlation)
    start_index, end_index, parameters, _, reduced_chi_squared = select_fit_window(
        time,
        correlation,
        correlation_error,
        config.beta,
        initial_gap,
    )
    gap = float(parameters[1])
    rng = np.random.default_rng(bootstrap_seed)
    gap_error, excited_error, successful_bootstraps = bootstrap_fit(
        energy_blocks,
        correlation_blocks,
        time,
        config.beta,
        start_index,
        end_index,
        gap,
        config.bootstrap_samples,
        rng,
    )
    excited_energy = energy + gap
    curve = periodic_correlation_model(time, parameters[0], gap, config.beta)
    fit = FitResult(
        gap=gap,
        gap_error=gap_error,
        fit_start=float(time[start_index]),
        fit_end=float(time[end_index]),
        reduced_chi_squared=reduced_chi_squared,
        successful_bootstraps=successful_bootstraps,
        curve=curve,
    )

    instanton_rate = None
    instanton_rate_error = None
    threshold = None
    if potential.family == "double_well":
        instanton_blocks = block_means(instantons, config.block_size) / config.beta
        instanton_rate = float(np.mean(instanton_blocks))
        instanton_rate_error = float(
            np.std(instanton_blocks, ddof=1) / np.sqrt(instanton_blocks.size)
        )
        threshold = 0.5 * potential.minimum_position

    return StudyResult(
        potential=potential,
        config=config,
        energy=energy,
        energy_error=energy_error,
        gap=gap,
        gap_error=gap_error,
        excited_energy=excited_energy,
        excited_energy_error=excited_error,
        acceptance=acceptance,
        proposal_width=proposal_width,
        tau_energy=tau_energy,
        reference=reference,
        density_samples=density_samples,
        correlation_time=time,
        correlation=correlation,
        correlation_error=correlation_error,
        fit=fit,
        energy_blocks=energy_blocks,
        correlation_blocks=correlation_blocks,
        instanton_rate=instanton_rate,
        instanton_rate_error=instanton_rate_error,
        representative_path=chain_data[0][6],
        instanton_threshold=threshold,
    )


def run_study(
    potential: Potential,
    config: SimulationConfig,
    seed: np.random.SeedSequence,
) -> StudyResult:
    """Run and analyse all independent chains for one potential."""

    children = seed.spawn(config.n_chains + 1)
    chain_data = [
        run_chain(potential, config, index, children[index])
        for index in range(config.n_chains)
    ]
    return analyse_ensemble(potential, config, chain_data, children[-1])


def concise_config(config: SimulationConfig) -> SimulationConfig:
    """Return a smaller configuration used for convergence diagnostics."""

    if config.measurements < 1_000:
        return replace(
            config,
            thermalization_sweeps=600,
            measurements=600,
            n_chains=2,
            block_size=30,
            bootstrap_samples=80,
        )
    measurements = max(1_600, config.measurements * 4 // 5)
    block_size = 50
    return replace(
        config,
        thermalization_sweeps=max(1_200, config.thermalization_sweeps * 3 // 5),
        measurements=measurements,
        n_chains=4,
        block_size=block_size,
        bootstrap_samples=max(100, config.bootstrap_samples // 2),
    )


def run_convergence_study(
    config: SimulationConfig,
    seed: np.random.SeedSequence,
) -> dict[str, list[tuple[float, float, float]]]:
    """Check lattice-spacing and Euclidean-extent convergence for the HO."""

    potential = harmonic_potential(k=1.0, mass=1.0)
    base = concise_config(config)
    lattice_settings = [(0.25, 32), (0.125, 64), (0.0625, 128)]
    beta_for_lattice = 8.0
    beta_settings = [6.0, 10.0, 16.0]
    spacing_for_beta = 0.125
    children = seed.spawn(len(lattice_settings) + len(beta_settings))

    lattice_results = []
    for index, (spacing, n_sites) in enumerate(lattice_settings):
        current = replace(base, beta=beta_for_lattice, n_sites=n_sites)
        result = run_study(potential, current, children[index])
        lattice_results.append((spacing, result.energy, result.energy_error))

    beta_results = []
    offset = len(lattice_settings)
    for index, beta in enumerate(beta_settings):
        n_sites = int(round(beta / spacing_for_beta))
        if n_sites % 2:
            n_sites += 1
        current = replace(base, beta=beta, n_sites=n_sites)
        result = run_study(potential, current, children[offset + index])
        beta_results.append((beta, result.energy, result.energy_error))

    return {"lattice": lattice_results, "extent": beta_results}


def format_estimate(value: float, error: float) -> str:
    return f"{value:.6f} +/- {error:.6f}"


def print_results(results: Iterable[StudyResult]) -> None:
    print("\nQuantum oscillator PIMC results")
    print(
        "model | parameters | E0 PIMC | E0 ref | gap PIMC | gap ref | "
        "E1 PIMC | E1 ref | accept | tau_E"
    )
    print("-" * 170)
    for result in results:
        reference_gap = result.reference.energies[1] - result.reference.energies[0]
        print(
            f"{result.potential.family:11s} | "
            f"{result.potential.parameter_text:27s} | "
            f"{format_estimate(result.energy, result.energy_error):23s} | "
            f"{result.reference.energies[0]:.6f} | "
            f"{format_estimate(result.gap, result.gap_error):23s} | "
            f"{reference_gap:.6f} | "
            f"{format_estimate(result.excited_energy, result.excited_energy_error):23s} | "
            f"{result.reference.energies[1]:.6f} | "
            f"{result.acceptance:.3f} | {result.tau_energy:.2f}"
        )
        print(
            f"  fit: tau=[{result.fit.fit_start:.3f}, {result.fit.fit_end:.3f}], "
            f"chi2/dof={result.fit.reduced_chi_squared:.2f}, "
            f"bootstraps={result.fit.successful_bootstraps}"
        )
        if result.instanton_rate is not None:
            print(
                "  Euclidean instanton rate: "
                f"{format_estimate(result.instanton_rate, result.instanton_rate_error or 0.0)}"
            )


def print_convergence(convergence: dict[str, list[tuple[float, float, float]]]) -> None:
    print("\nHarmonic-oscillator convergence checks")
    print("Lattice spacing at beta=8")
    for spacing, energy, error in convergence["lattice"]:
        print(f"  a={spacing:7.4f}: E0={format_estimate(energy, error)}")
    print("Euclidean extent at a=0.125")
    for beta, energy, error in convergence["extent"]:
        print(f"  beta={beta:5.1f}: E0={format_estimate(energy, error)}")


def plot_energy_comparison(results: list[StudyResult]) -> plt.Figure:
    labels = ["H1", "H2", "H3", "A1", "A2", "A3", "D1", "D2", "D3"]
    positions = np.arange(len(results))
    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for axis, level in zip(axes, (0, 1), strict=True):
        estimates = [r.energy if level == 0 else r.excited_energy for r in results]
        errors = [r.energy_error if level == 0 else r.excited_energy_error for r in results]
        references = [r.reference.energies[level] for r in results]
        axis.errorbar(
            positions - 0.08,
            estimates,
            yerr=errors,
            fmt="o",
            capsize=3,
            label="PIMC",
        )
        axis.scatter(
            positions + 0.08,
            references,
            marker="x",
            s=55,
            label="Hamiltonian benchmark",
        )
        axis.set_ylabel(f"E{level}")
        axis.grid(alpha=0.25)
        axis.legend()
    axes[-1].set_xticks(positions, labels)
    axes[-1].set_xlabel("Parameter set (H: harmonic, A: anharmonic, D: double well)")
    figure.suptitle("PIMC energies and independent Hamiltonian benchmarks")
    figure.tight_layout()
    return figure


def plot_ground_densities(results: list[StudyResult]) -> plt.Figure:
    figure, axes = plt.subplots(3, 3, figsize=(13, 10))
    for axis, result in zip(axes.ravel(), results, strict=True):
        reference = result.reference
        lower, upper = np.quantile(result.density_samples, [0.002, 0.998])
        axis.hist(
            result.density_samples,
            bins=70,
            range=(lower, upper),
            density=True,
            alpha=0.55,
            color="tab:blue",
            label="PIMC",
        )
        mask = (reference.x >= lower) & (reference.x <= upper)
        axis.plot(
            reference.x[mask],
            reference.ground_density[mask],
            color="tab:red",
            linewidth=1.8,
            label=r"$|\psi_0|^2$",
        )
        axis.set_title(result.potential.parameter_text, fontsize=9)
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    for axis in axes[-1]:
        axis.set_xlabel("x")
    for axis in axes[:, 0]:
        axis.set_ylabel("Probability density")
    figure.suptitle("Ground-state probability densities")
    figure.tight_layout()
    return figure


def representative_results(results: list[StudyResult]) -> list[StudyResult]:
    return [results[0], results[3], results[6]]


def plot_correlations(results: list[StudyResult]) -> plt.Figure:
    selected = representative_results(results)
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for axis, result in zip(axes, selected, strict=True):
        time = result.correlation_time
        normalized = result.correlation / result.correlation[0]
        normalized_error = result.correlation_error / result.correlation[0]
        fit_curve = result.fit.curve / result.correlation[0]
        axis.errorbar(
            time,
            normalized,
            yerr=normalized_error,
            fmt=".",
            markersize=3,
            alpha=0.65,
            label="PIMC correlation",
        )
        fit_mask = (time >= result.fit.fit_start) & (time <= result.fit.fit_end)
        axis.plot(time[fit_mask], fit_curve[fit_mask], color="tab:red", label="periodic fit")
        axis.axvspan(result.fit.fit_start, result.fit.fit_end, alpha=0.08, color="tab:red")
        axis.set_yscale("log")
        axis.set_ylim(bottom=max(1e-4, np.min(normalized[normalized > 0]) * 0.5))
        axis.set_title(result.potential.family.replace("_", " ").title())
        axis.set_xlabel(r"$\tau$")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel(r"$C(\tau)/C(0)$")
    axes[0].legend(fontsize=8)
    figure.suptitle("Imaginary-time correlations and first-gap fits")
    figure.tight_layout()
    return figure


def plot_tunnelling(results: list[StudyResult]) -> plt.Figure:
    double_results = [r for r in results if r.potential.family == "double_well"]
    figure, axes = plt.subplots(2, 1, figsize=(11, 7))
    example = double_results[0]
    time = np.arange(example.config.n_sites) * example.config.spacing
    axes[0].plot(time, example.representative_path, linewidth=1.2)
    threshold = example.instanton_threshold or 0.0
    axes[0].axhline(threshold, color="tab:red", linestyle="--", linewidth=1)
    axes[0].axhline(-threshold, color="tab:red", linestyle="--", linewidth=1)
    axes[0].set_ylabel("x")
    axes[0].set_xlabel(r"$\tau$")
    axes[0].set_title("Representative periodic double-well path")
    axes[0].grid(alpha=0.25)

    positions = np.arange(len(double_results))
    rates = [r.instanton_rate for r in double_results]
    errors = [r.instanton_rate_error for r in double_results]
    labels = [r.potential.parameter_text for r in double_results]
    axes[1].bar(positions, rates, yerr=errors, capsize=4, color="tab:purple", alpha=0.75)
    axes[1].set_xticks(positions, labels, rotation=12, ha="right")
    axes[1].set_ylabel("Transitions per Euclidean time")
    axes[1].set_title("Thresholded Euclidean instanton density")
    axes[1].grid(axis="y", alpha=0.25)
    figure.tight_layout()
    return figure


def plot_convergence(
    convergence: dict[str, list[tuple[float, float, float]]]
) -> plt.Figure:
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    spacing, energy, error = map(np.asarray, zip(*convergence["lattice"], strict=True))
    axes[0].errorbar(spacing, energy, yerr=error, fmt="o-", capsize=4)
    axes[0].axhline(0.5, color="black", linestyle="--", label="exact")
    axes[0].set_xlabel("Lattice spacing a")
    axes[0].set_ylabel("E0")
    axes[0].set_title(r"Lattice-spacing stability at $\beta=8$")
    axes[0].legend()
    axes[0].grid(alpha=0.25)

    beta, energy, error = map(np.asarray, zip(*convergence["extent"], strict=True))
    axes[1].errorbar(beta, energy, yerr=error, fmt="o-", capsize=4)
    axes[1].axhline(0.5, color="black", linestyle="--", label="exact")
    axes[1].set_xlabel(r"Euclidean extent $\beta$")
    axes[1].set_ylabel("E0")
    axes[1].set_title("Euclidean-extent stability")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    figure.suptitle("Harmonic-oscillator convergence diagnostics")
    figure.tight_layout()
    return figure


def create_figures(
    results: list[StudyResult],
    convergence: dict[str, list[tuple[float, float, float]]],
) -> dict[str, plt.Figure]:
    return {
        "energy_comparison": plot_energy_comparison(results),
        "ground_state_densities": plot_ground_densities(results),
        "correlation_fits": plot_correlations(results),
        "double_well_tunnelling": plot_tunnelling(results),
        "convergence": plot_convergence(convergence),
    }


def save_figures(figures: dict[str, plt.Figure], directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name, figure in figures.items():
        figure.savefig(directory / f"{name}.png", dpi=180, bbox_inches="tight")


def quick_config(config: SimulationConfig) -> SimulationConfig:
    return replace(
        config,
        thermalization_sweeps=600,
        measurements=600,
        n_chains=2,
        block_size=30,
        bootstrap_samples=80,
    )


def run_default_analysis(
    config: SimulationConfig,
) -> tuple[list[StudyResult], dict[str, list[tuple[float, float, float]]]]:
    root_seed = np.random.SeedSequence(config.seed)
    children = root_seed.spawn(len(default_potentials()) + 1)
    results = []
    for index, potential in enumerate(default_potentials()):
        print(f"Running {potential.family}: {potential.parameter_text}")
        results.append(run_study(potential, config, children[index]))
    convergence = run_convergence_study(config, children[-1])
    return results, convergence


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2025, help="root random seed")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="use a reduced diagnostic run instead of publication settings",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="create figures without opening interactive windows",
    )
    parser.add_argument(
        "--save-dir",
        type=Path,
        help="optional directory in which to save PNG figures",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = SimulationConfig(seed=args.seed)
    if args.quick:
        config = quick_config(config)
    print("Quantum-oscillator Euclidean path-integral Monte Carlo")
    print(f"Seed: {config.seed}")
    print(
        f"beta={config.beta:g}, sites={config.n_sites}, "
        f"a={config.spacing:.5f}, chains={config.n_chains}"
    )
    print(
        f"thermalization={config.thermalization_sweeps}, "
        f"measurements={config.measurements}, "
        f"sweeps/measurement={config.sweeps_between}"
    )
    results, convergence = run_default_analysis(config)
    print_results(results)
    print_convergence(convergence)
    figures = create_figures(results, convergence)
    if args.save_dir is not None:
        save_figures(figures, args.save_dir)
        print(f"\nSaved {len(figures)} figures to {args.save_dir}")
    if args.no_show:
        plt.close("all")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
