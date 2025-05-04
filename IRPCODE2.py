# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 20:50:12 2025

@author: kutadgu
"""

import numpy as np
import math

def F(q, a):
    
    # Computes the normalized spherical form factor F(q).

    # Parameters:
    #   q : Momentum transfer (1/m).
    #   a : Radius of the scatterer (m).

    # Returns: Value of F(q).
    
    qa = q * a
    # Handle the limit qa -> 0 to avoid division by zero
    
    if np.isclose(qa, 0):
        return 1.0
    return (3 * (np.sin(qa) - qa * np.cos(qa))) / (qa)**3

def scattering_amplitude(frequency, rho_medium, rho_scatterer,
                         lambda_medium, lambda_scatterer,
                         mu_medium, mu_scatterer,
                         c, a, theta):
    #
    # Computes the scattering amplitude using the first Born approximation, the First-Born Approximation
    # can be found on the Wu and Aki paper.

    # Parameters:
    #   frequency : Testing frequency (Hz). [Example: 1e6 for 1 MHz]
    #   rho_medium : Density of the medium (kg/m³).
    #   rho_scatterer : Density of the scatterer (kg/m³).
    #   lambda_medium : First Lamé parameter of the medium (Pa).
    #   lambda_scatterer : First Lamé parameter of the scatterer (Pa).
    #   mu_medium : Shear modulus of the medium (Pa).
    #   mu_scatterer : Shear modulus of the scatterer (Pa).
    #   c : Speed of the wave in the medium (m/s).
    #   a : Radius of the spherical scatterer (m).
    #   theta : Scattering angle (radians).

    # Returns:
    #   float : Scattering amplitude f(θ)


    # Angular frequency (omega) in rad/s
    omega = 2 * np.pi * frequency

    # Compute contrasts (difference between scatterer and medium properties)
    delta_rho = rho_scatterer - rho_medium        # Density contrast (kg/m³)
    delta_lambda = lambda_scatterer - lambda_medium # Lamé parameter contrast (Pa)
    delta_mu = mu_scatterer - mu_medium             # Shear modulus contrast (Pa)

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

    # Compute the scattering amplitude f(θ)
    term1 = (delta_rho / rho_medium) * np.cos(theta)
    term2 = delta_lambda / denom
    term3 = (delta_mu / denom) * (np.cos(theta) ** 2)
    
    f_theta = - (omega**2 / (4 * np.pi * c**2)) * (term1 + term2 - term3) * V * F_q

    return f_theta

def calculate_scatterer_mass_and_volume(
    target_n,           # Target number of scatterers per m³ (n)
    radius,             # Radius of a single scatterer
    density,            # Density of scatterer material
    radius_unit='m',    # 'm', 'cm', 'mm', or 'micron'
    density_unit='kg/m3', # 'kg/m3' or 'g/cm3'
    cup_volume_ml=75    # Testing volume in mL
):
    # Convert radius to meters
    if radius_unit == 'cm':
        radius_m = radius * 0.01
    elif radius_unit == 'mm':
        radius_m = radius * 0.001
    elif radius_unit == 'micron':
        radius_m = radius * 1e-6
    else:
        radius_m = radius  # Assume meters
    
    # Convert density to kg/m³
    if density_unit == 'g/cm3':
        density_kgm3 = density * 1000  # 1 g/cm³ = 1000 kg/m³
    else:
        density_kgm3 = density
    
    # Calculate volume of one scatterer (sphere)
    volume_per_scatterer = (4/3) * math.pi * (radius_m ** 3)
    
    # Mass of one scatterer (in kg)
    mass_per_scatterer_kg = density_kgm3 * volume_per_scatterer
    
    # Total number of scatterers in the test volume
    cup_volume_m3 = cup_volume_ml * 1e-6  # mL to m³ conversion
    total_scatterers = target_n * cup_volume_m3
    
    # Total mass required (in kg) and convert to grams
    total_mass_kg = total_scatterers * mass_per_scatterer_kg
    total_mass_grams = total_mass_kg * 1000
    
    # Calculate the volume of the scatterer material (in m³ and then mL)
    total_material_volume_m3 = total_mass_kg / density_kgm3
    total_material_volume_ml = total_material_volume_m3 * 1e6
    
    return total_mass_grams, total_material_volume_ml

if __name__ == "__main__":
    # ---- Born Approximation Parameters ----
    # Input material properties:
    frequency = 2.5e6          # Hz
    rho_medium = 2300          # kg/m³
    lambda_medium = 68e9       # Pa
    mu_medium = 79e9           # Pa
    c = 1000                 # m/s (given)
    
    # Scatterer properties:
    rho_scatterer = 3000       # kg/m³
    lambda_scatterer = 150e9   # Pa
    mu_scatterer = 60e9        # Pa
    
    # Scatterer size and scattering angle:
    # Here we use a grain diameter of 0.0173 mm so that radius is half of that in meters.
    # (0.0173 mm = 0.0173/1000 m; therefore, radius a = (0.0173/2)/1000 m)
    a = (0.0173 / 2) / 1000    # m
    theta = np.pi              # Backscattering

    # Calculate scattering amplitude using the Born approximation:
    f_theta = scattering_amplitude(frequency, rho_medium, rho_scatterer,
                                   lambda_medium, lambda_scatterer,
                                   mu_medium, mu_scatterer,
                                   c, a, theta)
    print("Scattering amplitude f(θ):", f_theta)
    A_squared = f_theta**2
    print("Scattering amplitude |A|^2:", A_squared)
    print("Backscattering Coefficient:", 1e13 * A_squared)

    # ---- Particle Number and Mass/Volume Calculation ----
    # We want to achieve: n * |A|^2 between 1e-3 and 5e-3.
    prod_min = 1e-3
    prod_max = 5e-3

    # Compute the range for n:
    n_min = prod_min / A_squared
    n_max = prod_max / A_squared
    n_baseline = (n_min + n_max) / 2  # A baseline value in between

    print("\nRequired number density n to achieve n*|A|^2 in the desired range:")
    print("Minimum n =", n_min, "m⁻³")
    print("Maximum n =", n_max, "m⁻³")
    print("Baseline n =", n_baseline, "m⁻³")

    # Now, using the baseline number density, compute the total mass and volume
    # For the scatterer mass calculation, we use the scatterer density as 3000 kg/m³.
    # Since the calculate_scatterer_mass_and_volume function accepts the radius unit,
    # we pass our computed radius (in meters) and specify 'm' as the radius unit.
    mass_grams, volume_ml = calculate_scatterer_mass_and_volume(
        target_n=n_baseline,
        radius=a,             # in meters
        density=3000,         # kg/m³ for scatterer material
        radius_unit='m',
        density_unit='kg/m3',
        cup_volume_ml=75      # test volume in mL
    )

    print("\nFor the baseline number density:")
    print("Total scatterer mass =", mass_grams, "grams")
    print("Total scatterer volume =", volume_ml, "mL")
