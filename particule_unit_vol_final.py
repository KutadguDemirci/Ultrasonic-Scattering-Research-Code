# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 18:23:49 2025

@author: kutadgu
"""

import math

def calculate_scatterer_mass_and_volume(
    target_n,           # Target number of scatterers per m³ (e.g., 1e6)
    radius,             # Radius of a single scatterer (e.g., 1e-6 for 1 micron)
    density,            # Density of scatterer material (e.g., 3200 for SiC in kg/m³)
    radius_unit='m',    # 'm', 'cm', 'mm', or 'micron'
    density_unit='kg/m3', # 'kg/m3' or 'g/cm3'
    cup_volume_ml=75    # Default 75 mL
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
    
    # Total number of scatterers needed in the cup
    cup_volume_m3 = cup_volume_ml * 1e-6  # Convert mL to m³ (1 mL = 1e-6 m³)
    total_scatterers = target_n * cup_volume_m3
    
    # Total mass required (in kg)
    total_mass_kg = total_scatterers * mass_per_scatterer_kg
    total_mass_grams = total_mass_kg * 1000  # Convert kg to grams
    
    # Calculate volume of the material corresponding to the total mass
    # Volume (in m³) = mass (in kg) / density (in kg/m³)
    total_material_volume_m3 = total_mass_kg / density_kgm3
    # Convert m³ to mL (1 m³ = 1e6 mL)
    total_material_volume_ml = total_material_volume_m3 * 1e6
    
    return total_mass_grams, total_material_volume_ml

# Example usage for SiC scatterers:
target_n = 4.6535e12  # 1e15 scatterers per m³
radius = 17.3/2     # 9.3 microns (example value)
density = 3      # 3 g/cm³ for SiC

mass_grams, volume_ml = calculate_scatterer_mass_and_volume(
    target_n, 
    radius, 
    density, 
    radius_unit='micron', 
    density_unit='g/cm3'
)

print(f"Required mass for {target_n:.1e} scatterers/m³: {mass_grams:.4f} grams")
print(f"This mass of material corresponds to a volume of: {volume_ml:.4f} mL")
