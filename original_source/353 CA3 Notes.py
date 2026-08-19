import numpy as np
import matplotlib.pyplot as plt

def MetroGauss(V, mass=1, x0=0, epsilon=1, beta=1, T=100000, Therm=5000):    
    xs = np.zeros(T)  
    x = x0  
    accepted = 0  
    for i in range(T):  
        delta = 2 * epsilon * (np.random.rand() - 0.5)  
        xprime = x + delta 
        dH = V(xprime) - V(x)  
        P = min(1, np.exp(-beta * dH))  
        if np.random.rand() <= P:  
            x = xprime  
            accepted += 1  
        xs[i] = x  
    acceptance_rate = accepted / T
    print(f"Acceptance Rate: {acceptance_rate:.3f}")
    samples = xs[Therm:]
    return samples

def BinAnalyse(samples, binsize=20):
    L = len(samples) // binsize  
    binList = [np.mean(samples[i * binsize:(i + 1) * binsize]) for i in range(L)] 
    mean_value = np.mean(samples)  
    uncertainty = np.std(binList) / np.sqrt(L - 1)
    return mean_value, uncertainty

def AutoCorrelation(samples, max_lag=100):
    mean_x = np.mean(samples)  
    var_x = np.var(samples)  
    autocorr = np.correlate(samples - mean_x, samples - mean_x, mode='full')    
    autocorr = autocorr[len(samples)-1:] / (var_x * len(samples))  
    tau = 1 + 2 * sum(autocorr[1:max_lag])
    print(f"Estimated autocorrelation time: {tau:.2f}")
    
    return autocorr[:max_lag], tau

def plot_autocorrelation_comparison():
    V_harm = lambda x: HarmonicOscillator(x, k=1.0, mass=1.0)
    V_anharm = lambda x: AnharmonicOscillator(x, k=1.0, lambda_=0.1, mass=1.0)
    V_double = lambda x: DoubleWellPotential(x, a=0.25, b=-1.0, mass=1.0)
    
    samples_harm = MetroGauss(V_harm, mass=1.0, epsilon=0.5, beta=1.0, T=50000, Therm=5000)
    samples_anharm = MetroGauss(V_anharm, mass=1.0, epsilon=0.5, beta=1.0, T=50000, Therm=5000)
    samples_double = MetroGauss(V_double, mass=1.0, epsilon=0.5, beta=1.0, T=50000, Therm=5000)
    
    autocorr_harm, tau_harm = AutoCorrelation(samples_harm, max_lag=100)
    autocorr_anharm, tau_anharm = AutoCorrelation(samples_anharm, max_lag=100)
    autocorr_double, tau_double = AutoCorrelation(samples_double, max_lag=100)
    
    plt.figure(figsize=(10, 6))
    plt.plot(autocorr_harm, 'r-', label=f'Harmonic (τ={tau_harm:.2f})')
    plt.plot(autocorr_anharm, 'g-', label=f'Anharmonic (τ={tau_anharm:.2f})')
    plt.plot(autocorr_double, 'b-', label=f'Double Well (τ={tau_double:.2f})')
    plt.xlabel("Lag")
    plt.ylabel("Autocorrelation")
    plt.title("Autocorrelation Comparison Between Potentials")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return tau_harm, tau_anharm, tau_double

def CalculateEnergy(samples, V, potential_type="harmonic", k=1.0, lambda_=0.1, mass=1.0): 
    
    x_squared = samples**2     
    x_squared_mean, x_squared_err = BinAnalyse(x_squared)
    
    x_fourth = samples**4
    x_fourth_mean, x_fourth_err = BinAnalyse(x_fourth)
    
    if potential_type == "harmonic":
        omega = np.sqrt(k/mass)
        E0 = 0.5 * k * x_squared_mean
        dE0 = x_squared_err * 0.5 * k 
        
    elif potential_type == "anharmonic":
        E0 = 0.5 * k * x_squared_mean + 3 * lambda_ * x_fourth_mean
        dE0 = np.sqrt((0.5*k*x_squared_err)**2 + (3*lambda_*x_fourth_err)**2)  
        
    elif potential_type == "double_well": 
        potential_energies = np.array([V(x) for x in samples])
        V_mean, V_err = BinAnalyse(potential_energies) 
        E0 = V_mean + 0.5 * mass * x_squared_mean # Vmean will be negative and dominate. This is the case when potentials have regions below zero, so E0 MC is allowed to be negative in this case.
        dE0 = np.sqrt(V_err**2 + (0.5*mass*x_squared_err)**2) 
    
    return E0, dE0

def ComputeExcitedStateEnergy(paths, dtau=0.1, potential_type="harmonic", omega=1.0, E0=0.5): 
    
    N = paths.shape[1] 
    C_tau = np.zeros(N // 2) 

    for tau in range(N // 2):
        C_tau[tau] = np.mean(paths[:, 0] * paths[:, tau]) 
    
    C_tau = C_tau / C_tau[0] 
    
    # Using the improved formula for energy gap calculation:
    # dE = -1/a log(<x(0)x(tau+a)> / <x(0)x(tau)>)
    log_C = np.log(np.abs(C_tau) + 1e-10)
    start_idx = min(5, N // 4)
    end_idx = min(N // 3, len(log_C) - 1)
    
    # Calculate energy gaps using the log ratio method
    dE_values = []
    for tau in range(start_idx, end_idx-1):
        if C_tau[tau+1] != 0 and C_tau[tau] != 0:
            dE = -np.log(abs(C_tau[tau+1] / C_tau[tau])) / dtau
            dE_values.append(dE)
    
    if len(dE_values) > 0:
        E1_minus_E0 = np.mean(dE_values)
        dE1 = np.std(dE_values) / np.sqrt(len(dE_values))
        E1 = E0 + E1_minus_E0
    else:
        # Fallback if log ratio method fails
        slopes = -(log_C[start_idx+1:end_idx] - log_C[start_idx:end_idx-1]) / dtau
        E1_minus_E0 = np.mean(slopes)
        dE1 = np.std(slopes) / np.sqrt(len(slopes))
        E1 = E0 + E1_minus_E0
    
    if potential_type == "harmonic":
        theoretical_gap = omega
        print(f"Theoretical Gap: {theoretical_gap:.4f}, Measured Gap: {E1_minus_E0:.4f} ± {dE1:.4f}")
        # Calculate relative difference between theoretical and measured gap
        gap_difference = abs(theoretical_gap - E1_minus_E0) / theoretical_gap * 100
        print(f"Gap Difference: {gap_difference:.2f}%")
    
    print(f"Energy Gap (E1-E0): {E1_minus_E0:.4f} ± {dE1:.4f}")
    print(f"First Excited State Energy: E1 = {E1:.4f} ± {dE1:.4f}")
    
    return E1, dE1
    
def QuantumTunneling(paths): 
    
    sign_changes = np.sum(np.diff(np.sign(paths), axis=1) != 0, axis=1)
    num_tunneling = np.mean(sign_changes) 
    uncertainty = np.std(sign_changes) / np.sqrt(len(sign_changes)) 
    alpha = 1 - np.exp(-num_tunneling * 0.1)  
    print(f"Tunneling Events: {num_tunneling:.2f} ± {uncertainty:.2f}")
    print(f"Tunneling Parameter (alpha): {alpha:.4f}")
    
    return num_tunneling, uncertainty, alpha

def HarmonicOscillator(x, k=1, mass=1):
    return 0.5 * k * x**2

def AnharmonicOscillator(x, k=1, lambda_=0.1, mass=1):
    return 0.5 * k * x**2 + lambda_ * x**4

def DoubleWellPotential(x, a=1, b=-1, mass=1):
    return a * x**4 + b * x**2  

def theoretical_harmonic_pdf(x, mass=1.0, k=1.0, beta=1.0):
    omega = np.sqrt(k/mass)
    return np.sqrt(mass*omega/(np.pi)) * np.exp(-mass*omega*x**2)

def theoretical_anharmonic_pdf(x, mass=1.0, k=1.0, lambda_=0.1, beta=1.0):
    
    V = 0.5 * k * x**2 + lambda_ * x**4
    unnormalized_pdf = np.exp(-beta * V)
    x_range = np.linspace(-10, 10, 1000)
    V_range = 0.5 * k * x_range**2 + lambda_ * x_range**4
    Z = np.trapz(np.exp(-beta * V_range), x_range)
    
    return unnormalized_pdf / Z

def theoretical_double_well_pdf(x, mass=1.0, a=1.0, b=-1.0, beta=1.0):
    
    V = a * x**4 + b * x**2
    unnormalized_pdf = np.exp(-beta * V)
    x_range = np.linspace(-10, 10, 1000)
    V_range = a * x_range**4 + b * x_range**2
    Z = np.trapz(np.exp(-beta * V_range), x_range)
    
    return unnormalized_pdf / Z

def calculate_normalized_pdf(V, x_range, beta=1.0):
    
    potentials = np.array([V(x) for x in x_range])
    unnormalized_pdf = np.exp(-beta * potentials)
    Z = np.trapz(unnormalized_pdf, x_range)
    
    return unnormalized_pdf / Z

def plot_histogram(samples, V, title, x_range=(-3, 3), bins=50, theoretical_pdf=None, beta=1.0):

    fig, ax = plt.subplots(figsize=(10, 5))
    hist, bin_edges, _ = ax.hist(samples, bins=bins, density=True, alpha=0.7, color='blue', edgecolor='black')
    x_vals = np.linspace(x_range[0], x_range[1], 1000)
    pdf_vals = [theoretical_pdf(x) for x in x_vals]
    ax.plot(x_vals, pdf_vals, 'r-', linewidth=2, label='Theoretical PDF')

    ax.legend()
    ax.set_xlabel("Position (x)")
    ax.set_ylabel("Probability Density")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# Theoretical values calculation functions
def theoretical_harmonic_E0(k, mass):
    omega = np.sqrt(k/mass)
    return 0.5 * omega

def theoretical_harmonic_E1(k, mass, dtau=0.1):
    # Using the formula E1 = E0 + w_0
    # where w_0 = mu/sqrt(m) * (1 - (mu^2 * a^2)/(24*m))
    # with mu = sqrt(k), a = dtau
    mu = np.sqrt(k)
    omega0 = mu / np.sqrt(mass) * (1 - (mu**2 * dtau**2) / (24 * mass))
    E0 = theoretical_harmonic_E0(k, mass)
    return E0 + omega0

def theoretical_anharmonic_E0(k, lambda_, mass, beta=1.0):
    # First-order perturbation theory
    omega = np.sqrt(k/mass)
    E0_harmonic = 0.5 * omega
    # The anharmonic correction at first order
    correction = 3 * lambda_ / (2 * mass**2 * omega**2)
    return E0_harmonic + correction

def theoretical_anharmonic_E1(k, lambda_, mass, dtau=0.1, beta=1.0):
    # Using the formula E1 = E0 + w_0
    # For anharmonic, adjust the frequency with perturbation correction
    mu = np.sqrt(k)
    # Base harmonic frequency
    omega = mu / np.sqrt(mass)
    # Anharmonic correction to frequency
    delta_omega = -15 * lambda_ / (4 * mass**2 * omega**2)
    # Corrected frequency with time discretization effect
    omega0 = (omega + delta_omega) * (1 - ((omega + delta_omega)**2 * dtau**2) / (24 * mass))
    
    E0 = theoretical_anharmonic_E0(k, lambda_, mass, beta)
    return E0 + omega0

def theoretical_double_well_E0(a, b, mass, beta=1.0):
    # For small tunneling, approximate as harmonic near one of the wells
    x_min = np.sqrt(-b/(2*a))  # Position of minimum
    # Second derivative at minimum gives effective spring constant
    k_eff = -2*b + 12*a*x_min**2
    omega_eff = np.sqrt(k_eff/mass)
    # Ground state is approximately harmonic oscillator at one minimum
    return 0.5 * omega_eff + a*x_min**4 + b*x_min**2

def theoretical_double_well_E1(a, b, mass, tunneling_parameter=None, dtau=0.1, beta=1.0):
    # For double well with tunneling
    x_min = np.sqrt(-b/(2*a))  # Position of minimum
    k_eff = -2*b + 12*a*x_min**2
    mu = np.sqrt(k_eff)
    
    # Basic harmonic approximation with discretization correction
    omega0 = mu / np.sqrt(mass) * (1 - (mu**2 * dtau**2) / (24 * mass))
    
    E0 = theoretical_double_well_E0(a, b, mass, beta)
    
    # If tunneling parameter is provided, use it for a better estimate
    if tunneling_parameter is not None:
        # Tunneling splitting approximation
        delta_E = -0.1 * np.log(1 - tunneling_parameter)
        return E0 + delta_E
    else:
        # Otherwise use corrected harmonic approximation
        return E0 + omega0

def explore_parameters(potential_func, param_values, name, beta=1.0, epsilon=0.5, T=100000, Therm=5000, dtau=0.1):
    
    results = []  
    
    for params in param_values:
        
        if name == "Harmonic Oscillator":
            k, mass = params
            V = lambda x: potential_func(x, k=k, mass=mass)
            omega = np.sqrt(k/mass)
            title = f"{name}: mass={mass}, k={k}, ω={omega:.3f}"
            theoretical_pdf = lambda x: theoretical_harmonic_pdf(x, mass=mass, k=k, beta=beta)
            potential_type = "harmonic"
            theory_E0 = theoretical_harmonic_E0(k, mass)
            theory_E1 = theoretical_harmonic_E1(k, mass, dtau)
            
        elif name == "Anharmonic Oscillator":
            k, lambda_, mass = params
            V = lambda x: potential_func(x, k=k, lambda_=lambda_, mass=mass)
            title = f"{name}: mass={mass}, k={k}, λ={lambda_}"
            theoretical_pdf = lambda x: theoretical_anharmonic_pdf(x, mass=mass, k=k, lambda_=lambda_, beta=beta)
            potential_type = "anharmonic"
            theory_E0 = theoretical_anharmonic_E0(k, lambda_, mass, beta)
            theory_E1 = theoretical_anharmonic_E1(k, lambda_, mass, dtau, beta)
            
        elif name == "Double Well Potential":
            a, b, mass = params
            V = lambda x: potential_func(x, a=a, b=b, mass=mass)
            title = f"{name}: mass={mass}, a={a}, b={b}"
            theoretical_pdf = lambda x: theoretical_double_well_pdf(x, mass=mass, a=a, b=b, beta=beta)
            potential_type = "double_well"
            theory_E0 = theoretical_double_well_E0(a, b, mass, beta)
            # For E1, we'll update later after calculating tunneling
            theory_E1 = None
        
        print(f"\n{title}")
        
        samples = MetroGauss(V, mass=mass, epsilon=epsilon, beta=beta, T=T, Therm=Therm)
        
        if potential_type == "harmonic":
            E0, dE0 = CalculateEnergy(samples, V, potential_type=potential_type, k=k, mass=mass)
            
        elif potential_type == "anharmonic":
            E0, dE0 = CalculateEnergy(samples, V, potential_type=potential_type, k=k, lambda_=lambda_, mass=mass)
            
        else:
            E0, dE0 = CalculateEnergy(samples, V, potential_type=potential_type, k=abs(b), mass=mass)
            
        print(f"Ground State Energy: E0 = {E0:.4f} ± {dE0:.4f}")
        print(f"Theoretical Ground State: E0_theory = {theory_E0:.4f}")
        diff_percent_E0 = 100 * (E0 - theory_E0) / theory_E0
        print(f"E0 Difference: {diff_percent_E0:.2f}%")

        plot_histogram(samples, V, title, x_range=(-3,3), theoretical_pdf=theoretical_pdf, beta=beta)
        paths = np.reshape(samples, (-1, 100))        

        if potential_type == "harmonic":
            omega = np.sqrt(k/mass)
            E1, dE1 = ComputeExcitedStateEnergy(paths, dtau=dtau, potential_type=potential_type, omega=omega, E0=E0)
            
        else:
            E1, dE1 = ComputeExcitedStateEnergy(paths, dtau=dtau, potential_type=potential_type, E0=E0)
            
        if name == "Double Well Potential":
            tunneling, tunneling_err, alpha = QuantumTunneling(paths)
            E1_corrected = E0 - (0.1) * np.log(1-alpha)  
            print(f"Corrected E1 (using tunneling): {E1_corrected:.4f}")
            # Update the theoretical E1 with tunneling info
            theory_E1 = theoretical_double_well_E1(a, b, mass, tunneling_parameter=alpha, dtau=dtau, beta=beta)
            print(f"Theoretical Excited State: E1_theory = {theory_E1:.4f}")
            diff_percent_E1 = 100 * (E1_corrected - theory_E1) / theory_E1
            print(f"E1 Difference: {diff_percent_E1:.2f}%")
            results.append((params, E0, dE0, E1_corrected, dE1, tunneling, tunneling_err, 
                            theory_E0, theory_E1, diff_percent_E0, diff_percent_E1))
            
        else:
            print(f"Theoretical Excited State: E1_theory = {theory_E1:.4f}")
            diff_percent_E1 = 100 * (E1 - theory_E1) / theory_E1
            print(f"E1 Difference: {diff_percent_E1:.2f}%")
            results.append((params, E0, dE0, E1, dE1, theory_E0, theory_E1, diff_percent_E0, diff_percent_E1))
    
    return results

# Set the time discretization parameter
dtau = 0.1

print("\n=== Harmonic Oscillator Study ===")
harmonic_params = [
    (1.0, 1.0),    
    (2.0, 1.0),   
    (1.0, 0.5)     
]    
harmonic_results = explore_parameters(HarmonicOscillator, harmonic_params, "Harmonic Oscillator", beta=1.0, dtau=dtau)

print("\n=== Harmonic Oscillator Results ===")
print(f"{'k':^8} | {'mass':^8} | {'E0 (MC)':^15} | {'E0 (Theory)':^15} | {'Diff %':^10} | {'E1 (MC)':^15} | {'E1 (Theory)':^15} | {'Diff %':^10}")
print("-" * 110)
for result in harmonic_results:
    params, E0, dE0, E1, dE1, theory_E0, theory_E1, diff_E0, diff_E1 = result
    k, mass = params
    print(f"{k:^8.2f} | {mass:^8.2f} | {E0:^6.4f} ± {dE0:^6.4f} | {theory_E0:^15.4f} | {diff_E0:^10.2f} | {E1:^6.4f} ± {dE1:^6.4f} | {theory_E1:^15.4f} | {diff_E1:^10.2f}")

print("\n=== Anharmonic Oscillator Study ===")
anharmonic_params = [
    (1.0, 0.1, 1.0),    
    (1.0, 0.5, 1.0),    
    (1.0, 0.1, 0.5)     
]
anharmonic_results = explore_parameters(AnharmonicOscillator, anharmonic_params, "Anharmonic Oscillator", beta=1.0, dtau=dtau)

print("\n=== Anharmonic Oscillator Results ===")
print(f"{'k':^8} | {'λ':^8} | {'mass':^8} | {'E0 (MC)':^15} | {'E0 (Theory)':^15} | {'Diff %':^10} | {'E1 (MC)':^15} | {'E1 (Theory)':^15} | {'Diff %':^10}")
print("-" * 125)
for result in anharmonic_results:
    params, E0, dE0, E1, dE1, theory_E0, theory_E1, diff_E0, diff_E1 = result
    k, lambda_, mass = params
    print(f"{k:^8.2f} | {lambda_:^8.2f} | {mass:^8.2f} | {E0:^6.4f} ± {dE0:^6.4f} | {theory_E0:^15.4f} | {diff_E0:^10.2f} | {E1:^6.4f} ± {dE1:^6.4f} | {theory_E1:^15.4f} | {diff_E1:^10.2f}")

print("\n=== Double Well Potential Study ===")
doublewell_params = [
    (0.25, -1.0, 1.0),
    (0.5, -1.0, 1.0),     
    (0.25, -1.0, 0.5)     
]
doublewell_results = explore_parameters(DoubleWellPotential, doublewell_params,"Double Well Potential", beta=1.0, dtau=dtau)

print("\n=== Double Well Potential Results ===")
print(f"{'a':^8} | {'b':^8} | {'mass':^8} | {'E0 (MC)':^15} | {'E0 (Theory)':^15} | {'Diff %':^10} | {'E1 (MC)':^15} | {'E1 (Theory)':^15} | {'Diff %':^10} | {'Tunnel':^10}")
print("-" * 140)
for result in doublewell_results:
    params, E0, dE0, E1, dE1, tunneling, tunneling_err, theory_E0, theory_E1, diff_E0, diff_E1 = result
    a, b, mass = params
    print(f"{a:^8.2f} | {b:^8.2f} | {mass:^8.2f} | {E0:^6.4f} ± {dE0:^6.4f} | {theory_E0:^15.4f} | {diff_E0:^10.2f} | {E1:^6.4f} ± {dE1:^6.4f} | {theory_E1:^15.4f} | {diff_E1:^10.2f} | {tunneling:^10.2f}")

plot_autocorrelation_comparison()

