import numpy as np

def F(q, a):
    """
    Computes the normalized spherical form factor F(q).

    Parameters:
      q : float
          Momentum transfer (1/m).
      a : float
          Radius of the scatterer (m).

    Returns:
      float : Value of F(q).
    """
    qa = q * a
    # Handle the limit qa -> 0 to avoid division by zero
    if np.isclose(qa, 0):
        return 1.0
    return (3 * (np.sin(qa) - qa * np.cos(qa))) / (qa)**3

def scattering_amplitude(frequency,  rho_medium, rho_scatterer,
                         lambda_medium, lambda_scatterer,
                         mu_medium, mu_scatterer,
                         c, a, theta):
    
    # Computes the scattering amplitude f(θ) using the first Born approximation.

    # Parameters:
    #   frequency : Testing frequency (Hz). [Example: 1e6 for 1 MHz]
    #   rho_medium : Density of the medium (kg/m³). [Example: 2700 kg/m³ for aluminum]
    #   rho_scatterer : Density of the scatterer (kg/m³). [Example: 2800 kg/m³]
    #   lambda_medium : First Lamé parameter of the medium (Pa). [Units: Pascals]
    #   lambda_scatterer : First Lamé parameter of the scatterer (Pa). [Units: Pascals]
    #   mu_medium : Shear modulus (second Lamé parameter) of the medium (Pa). [Units: Pascals]
    #   mu_scatterer : Shear modulus (second Lamé parameter) of the scatterer (Pa). [Units: Pascals]
    #   c : Speed of the wave in the medium (m/s). [Example: 6420 m/s]
    #   a : Radius of the spherical scatterer (m). [Example: 0.001 m for 1 mm radius]
    #   theta : Scattering angle (radians). [For backscattering, use theta = π]

    # Returns:
    #   float : Scattering amplitude
    
    
    # Angular frequency (omega) in rad/s
    omega = 2 * np.pi * frequency

    # Compute contrasts (difference between scatterer and medium properties)
    delta_rho = rho_scatterer - rho_medium                     # Density contrast (kg/m³)
    delta_lambda = lambda_scatterer - lambda_medium              # Lamé parameter contrast (Pa)
    delta_mu = mu_scatterer - mu_medium                          # Shear modulus contrast (Pa)

    # Volume of the spherical scatterer, V = (4/3)πa³ (m³)
    V = (4 / 3) * np.pi * a**3

    # Compute wave number in the medium: k = omega / c (1/m)
    k = omega / c

    # Compute momentum transfer: q = 2k sin(θ/2) (1/m)
    q = 2 * k * np.sin(theta / 2)

    # Compute the normalized spherical form factor
    F_q = F(q, a)

    # Denominator term common in the elastic moduli (Pa)
    denom = lambda_medium + 2 * mu_medium

    # Compute the scattering amplitude f(θ) using the provided formula
    # f(θ) = -ω²/(4πc²) * [ (δρ/ρ_medium) cosθ + (δλ/denom) - (δμ/denom) cos²θ ] * V * F(q)
    term1 = (delta_rho / rho_medium) * np.cos(theta)
    term2 = delta_lambda / denom
    term3 = (delta_mu / denom) * (np.cos(theta) ** 2)
    
    f_theta = - (omega**2 / (4 * np.pi * c**2)) * (term1 + term2 - term3) * V * F_q

    return f_theta

# Example usage:
if __name__ == "__main__":
    # Input material properties:
    # Frequency
    frequency = 3.5e6            # Hz

    # Medium properties:
    rho_medium = 2300            # kg/m³
    lambda_medium = 68e9         # Pa (example value)
    mu_medium = 79e9             # Pa (example value)
    c = 1000                     # m/s

    # Scatterer properties:
    rho_scatterer = 3000         # kg/m³
    lambda_scatterer = 150e9     # Pa (example value)
    mu_scatterer = 60e9          # Pa (example value)

    # Scatterer size and scattering angle:
    a = 0.0173/2000              # m (1 mm radius)
    theta = np.pi                # radians (backscattering)

    f_theta = scattering_amplitude(frequency, rho_medium, rho_scatterer,
                                   lambda_medium, lambda_scatterer,
                                   mu_medium, mu_scatterer,
                                   c, a, theta)
    print("Scattering amplitude f(θ):", f_theta)
    print("Scattering amplitude |A|^2:", f_theta**2)
    print("Backscattering Coefficient :", 1e13*(f_theta**2))
    
    
    