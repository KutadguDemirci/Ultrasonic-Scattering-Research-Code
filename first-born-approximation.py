import numpy as np

def first_born_approximation(frequency, rho_med, rho_scat, c, a, 
                              E_med, nu_med, E_scat, nu_scat, theta):
    
    # Computes the first Born approximation scattering amplitude f(theta)
    # and its squared magnitude |A|^2 for elastic scattering.

    # Parameters:
    #   frequency : Testing frequency (Hz). Example: 1e6 for 1 MHz.
    #   rho_med : Density of the medium (kg/m³). Example: 2700 for aluminum.
    #   rho_scat : Density of the scatterer (kg/m³). Example: 2800.
    #   c : Speed of sound in the medium (m/s). Example: 6420 m/s.
    #   a : Radius of the spherical scatterer (m). Example: 0.001 for 1 mm.
    #   E_med : Young's modulus of the medium (Pa). Example: 70e9 for aluminum.
    #   nu_med : Poisson's ratio of the medium. Example: 0.33.
    #   E_scat : Young's modulus of the scatterer (Pa). Example: 77e9 (if 10% stiffer than medium).
    #   nu_scat : Poisson's ratio of the scatterer. Example: 0.33.
    #   theta : Scattering angle (radians). For backscattering, use theta = np.pi.

    # Returns:
    #   f : Scattering amplitude f(theta).
    #   A2 : Squared magnitude |f(theta)|^2.
          
    # Angular frequency (omega) in rad/s:
    omega = 2 * np.pi * frequency  # e.g., for 1 MHz, omega = 2*pi*1e6

    # Effective stiffness of the medium, C0 (for longitudinal waves):
    # C0 = E_med*(1 - nu_med) / ((1 + nu_med)*(1 - 2*nu_med))
    C0 = E_med * (1 - nu_med) / ((1 + nu_med) * (1 - 2 * nu_med))
    
    # Effective stiffness of the scatterer:
    C_scat = E_scat * (1 - nu_scat) / ((1 + nu_scat) * (1 - 2 * nu_scat))
    
    # Contrasts:
    Delta_rho = rho_scat - rho_med       # Density contrast (kg/m³)
    Delta_C = C_scat - C0                # Stiffness contrast (Pa)
    
    # Wave number in the medium: k = omega / c (m⁻¹)
    k = omega / c
    
    # Momentum transfer: q = 2k * sin(theta/2) (m⁻¹)
    q = 2 * k * np.sin(theta / 2)
    
    # Compute the normalized spherical form factor F(q)
    qa = q * a
    if np.isclose(qa, 0):
        F_q = 1.0  # For very small qa, the limit is 1.
    else:
        F_q = (3 * (np.sin(qa) - qa * np.cos(qa))) / (qa)**3

    # Volume of the spherical scatterer: V = (4/3) * pi * a^3 (m³)
    V = (4 / 3) * np.pi * a**3

    # Define the effective contrast function U, assumed constant inside the scatterer.
    # U = omega^2 * [Delta_rho + (Delta_C / C0) * (1 - cos^2(theta))]
    # (1 - cos^2(theta)) accounts for directional dependence in stiffness contrast.
    U = omega**2 * (Delta_rho + (Delta_C / C0) * (1 - np.cos(theta)**2))
    
    # Compute the scattering amplitude f(theta) using the first Born approximation:
    # f(theta) = -[U * V * F(q)] / (4*pi)
    f = - (U * V * F_q) / (4 * np.pi)
    
    # Compute the squared magnitude |f(theta)|^2 (backscattering intensity proxy)
    A2 = np.abs(f)**2
    
    return f, A2

if __name__ == "__main__":
    # ============================
    # Example Input Values:
    # ----------------------------
    # Testing frequency (Hz)
    frequency = 1e6        # 1 MHz
    
    # Medium properties:
    rho_med = 2300         # kg/m³, e.g., aluminum density.
    c = 1000               # m/s, speed of sound in aluminum.
    E_med = 70e9           # Pa, Young's modulus of aluminum.
    nu_med = 0.33          # Poisson's ratio for aluminum.
    
    # Scatterer properties:
    rho_scat = 3000        # kg/m³, e.g., a slightly different density.
    E_scat = 77e9          # Pa, e.g., 10% higher than the medium.
    nu_scat = 0.33         # Poisson's ratio for scatterer (assumed same here).
    
    # Grain (scatterer) size:
    a = 0.0173/2000              # m, i.e., 1 mm radius.
    
    # Scattering angle: for backscattering, set theta = pi (180°)
    theta = np.pi          # radians
    
    # ============================
    # Compute the scattering amplitude and squared amplitude:
    f, A2 = first_born_approximation(frequency, rho_med, rho_scat, c, a, 
                                     E_med, nu_med, E_scat, nu_scat, theta)
    
    # Print the results:
    print("Scattering amplitude f(theta):", f)
    print("Squared amplitude |A|^2:", A2)