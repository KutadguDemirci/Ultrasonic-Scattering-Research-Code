import math

def calculate_number_density(total_mass, sample_volume, grain_radius, grain_density):
    
    # Calculate the number density of scatterers (number per unit volume).

    # Parameters:
    #   total_mass: Total mass of the scatterers (kg)
    #   sample_volume: Volume of the sample in which scatterers are distributed (m³)
    #   grain_radius: Radius of one scatterer grain (m)
    #   grain_density: Density of the scatterer material (kg/m³)

    # Returns:
    #   number_density: Number of scatterers per unit volume (number/m³)
    #   total_scatterers: Total number of scatterers (dimensionless)
    
    
    # Calculate the volume of a single perfectly spherical grain
    grain_volume = (4/3) * math.pi * grain_radius**3

    # Calculate the mass of a single grain
    mass_per_grain = grain_density * grain_volume

    # Total number of scatterers
    total_scatterers = total_mass / mass_per_grain

    # Number density (number per unit volume)
    number_density = total_scatterers / sample_volume

    return number_density, total_scatterers

# Example usage:
if __name__ == "__main__":
    # Input parameters (you can change these values as needed)
    total_mass = 0.1          # Total mass of scatterers in kg
    sample_volume = 0.001     # Sample volume in m³ (e.g., 1 liter = 0.001 m³)
    grain_radius = 0.0005     # Grain radius in m (0.5 mm)
    grain_density = 3200      # Density of silicone carbide in kg/m³

    n, N_total = calculate_number_density(total_mass, sample_volume, grain_radius, grain_density)
    print("Number density of scatterers: {:.3e} scatterers/m³".format(n))
    print("Total number of scatterers: {:.3e}".format(N_total))
