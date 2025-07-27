import math
import re
# Constants
BOLTZMANN_CONSTANT = 1.38e-23  
PLANCKS_CONSTANT = 6.63e-34   
PI = 3.1416    # Pi
AVOGADRO_NUMBER = 6.023e23 

def calculate_atomic_mass(formula_str):
    atom_symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt"
    ]
    atom_masses = [
        1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897, 24.305,
        26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961,
        54.938, 55.845, 58.9332, 58.6934, 63.546, 65.39, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.8,
        85.4678, 87.62, 88.9059, 91.224, 92.9064, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411,
        114.818, 118.71, 121.76, 127.6, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.24,
        145, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967, 178.49,
        180.948, 183.84, 186.207, 190.23, 192.217, 195.078, 196.966, 200.59, 204.383, 207.2, 208.98, 209,
        210, 222, 223, 226, 227, 232.038, 231.036, 238.029, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259,
        262, 261, 262, 266, 264, 277, 268
    ]
    atoms = dict(zip(atom_symbols, atom_masses))
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula_str)
    total_mass = 0.0
    for symbol, count in matches:
        if symbol not in atoms:
            raise ValueError(f"Unknown atom symbol: {symbol}")
        n = int(count) if count else 1
        total_mass += atoms[symbol] * n
    return total_mass

def calculate_tetragonal_shear_modulus(c11, c12):
    return (c11 - c12) / 2

def calculate_cauchy_pressure(c12, c44):
    return (c12 - c44)

def calculate_kleixman_parameter(c11, c12):
    return (c11 + 8 * c12) / (7* c11 + 2 * c12)

def calculate_pughs_ratio(bulk_modulus, shear_modulus):
    return shear_modulus / bulk_modulus

def calculate_machinability(bulk_modulus, c44):
    return bulk_modulus / c44

def calculate_a1(c11, c13, c33, c44):
    return (4 * c44) / (c11 + c33 - 2 * c13)

def calculate_a2(c22, c23, c33, c55):
    return (4 * c55) / (c22 + c33 - 2 * c23)

def calculate_a3(a1, a2):
    return a1*a2

def calculate_zener_anisotropy_factor(c11, c12, c44):
    return (2 * c44) / (c11 - c12)

def calculate_c44_reuss(c11, c12, c44):
    upper_term = c44*(c11 - c12)
    lower_term = 3 * (c11 - c12) + 4 * c44
    return 5/3 * (upper_term / lower_term)

def calculate_c44_voigt(c44_reuss, c11, c12, c44):
    upper_term = (c11 - c12 - 2 * c44)
    lower_term = 3*(c11 - c12) + 4 * c44
    return c44_reuss + 3/5 * (upper_term**2 / lower_term)

def calculate_universal_log_euclidean_index(c44_voigt, c44_reuss, bulk_modulus_voigt, bulk_modulus_reuss):
    log_bulk_modulus_ratio = math.log(bulk_modulus_voigt/bulk_modulus_reuss)
    second_term = math.log(c44_voigt/c44_reuss)
    summed_term = log_bulk_modulus_ratio ** 2 + 5 * second_term ** 2
    return math.sqrt(summed_term)

def calculate_universal_anisotropy_index(bulk_modulus_voigt, bulk_modulus_reuss, shear_modulus_voigt, shear_modulus_reuss):
    bulk_modulus_ratio = (bulk_modulus_voigt/bulk_modulus_reuss)
    shear_modulus_ratio = (shear_modulus_voigt/shear_modulus_reuss)
    return 5 * shear_modulus_ratio + bulk_modulus_ratio - 6

def calculate_anisotropy_in_compressibility(bulk_modulus_voigt, bulk_modulus_reuss):
    upper_term = bulk_modulus_voigt - bulk_modulus_reuss
    lower_term = bulk_modulus_voigt + bulk_modulus_reuss
    ratio = upper_term / lower_term
    # Returns the value in percentage
    return ratio * 100

def calculate_anisotropy_in_shear(shear_modulus_voigt, shear_modulus_reuss, shear_modulus):
    return (shear_modulus_voigt - shear_modulus_reuss) / (2 * shear_modulus)

def calculate_equivalent_zener_anisotropy(universal_anisotropy_index):
    first_term = 1 + 5/12 * universal_anisotropy_index
    second_term = math.sqrt(first_term**2 - 1)
    return first_term + second_term

def calculate_uniaxial_bulk_modulus_Ba(c11, c12, c13, c22, c23, c33):
    uniaxial_alpha = calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33)
    uniaxial_beta = calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33)
    uniaxial_lambda = calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta)
    return uniaxial_lambda/(1 + uniaxial_alpha + uniaxial_beta)

def calculate_uniaxial_bulk_modulus_Bb(c11, c12, c13, c22, c23, c33):
    uniaxial_alpha = calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33)
    uniaxial_beta = calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33)
    uniaxial_lambda = calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta)
    uniaxial_bulk_modulus_Ba = uniaxial_lambda/(1 + uniaxial_alpha + uniaxial_beta)
    return uniaxial_bulk_modulus_Ba/uniaxial_alpha

def calculate_uniaxial_bulk_modulus_Bc(c11, c12, c13, c22, c23, c33):
    uniaxial_alpha = calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33)
    uniaxial_beta = calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33)
    uniaxial_lambda = calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta)
    uniaxial_bulk_modulus_Ba = uniaxial_lambda/(1 + uniaxial_alpha + uniaxial_beta)
    return uniaxial_bulk_modulus_Ba/uniaxial_beta

def calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33):
    upper_term = (c11 - c12) * (c33 -c13) - (c23 - c13) * (c11 - c13)
    lower_term = (c33 - c13) * (c22 - c12) - (c13 - c23) * (c12 - c23)
    return upper_term / lower_term

def calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33):
    upper_term = (c22 - c12) * (c11 - c13) - (c11 - c12) * (c23 - c12)
    lower_term = (c22 - c12) * (c33 - c13) - (c12 - c23) * (c13 - c23)
    return upper_term / lower_term

def calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta):
    return c11 + 2 * c12 * uniaxial_alpha + c22 * uniaxial_alpha**2 + 2 * c13 * uniaxial_beta + c33 * uniaxial_beta**2 + 2 * c33 * uniaxial_alpha * uniaxial_beta

def calculate_uniaxial_beta_relax(c11, c12, c13, c22, c23, c33):
    uniaxial_alpha = calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33)
    uniaxial_beta = calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33)
    uniaxial_lambda = calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta)
    return uniaxial_lambda/(1 + uniaxial_alpha + uniaxial_beta) ** 2

def calculate_anisotropic_bulk_modulus_ABa(uniaxial_bulk_modulus_Ba, uniaxial_bulk_modulus_Bc):
    return uniaxial_bulk_modulus_Ba/uniaxial_bulk_modulus_Bc

def calculate_anisotropic_bulk_modulus_ABc(uniaxial_bulk_modulus_Bb, uniaxial_bulk_modulus_Bc):
    return uniaxial_bulk_modulus_Bc/uniaxial_bulk_modulus_Bb

def calculate_linear_compressibility_beta_a(c11, c12, c13, c33):
    upper_term = c33 - c13
    lower_term = (c11 + c12) * c33 - 2 * (c13**2)
    return upper_term / lower_term

def calculate_linear_compressibility_beta_c(c11, c12, c13, c33):
    upper_term = c11 + c12 -  2 * c13
    lower_term = (c11 + c12) * c33 - 2 * (c13**2)
    return upper_term / lower_term

def calculate_ratio_of_linear_compressibility(linear_compressibility_beta_a, linear_compressibility_beta_c):
    return linear_compressibility_beta_c/linear_compressibility_beta_a

def calculate_transverse_velocity(shear_modulus, density):
    return math.sqrt(shear_modulus*1e9/density)

def calculate_longitudinal_velocity(bulk_modulus, shear_modulus, density):
    return math.sqrt((3* bulk_modulus * 1e9 + 4 * shear_modulus * 1e9)/(3 * density))

def calculate_average_velocity(longitudinal_velocity, transverse_velocity):
    summed_term = 2 / transverse_velocity**3 + 1/longitudinal_velocity**3
    multiplied_term =1/3 * summed_term
    average_velocity = pow(multiplied_term, -1/3)
    return average_velocity

def calculate_density(formula_factor, molar_mass, cell_volume):
    density = (formula_factor * molar_mass ) / ( AVOGADRO_NUMBER * cell_volume )
    return density

def calculate_acoustic_impedance(density, shear_modulus):
    return math.sqrt(density * shear_modulus*1e9)

def calculate_radiation_factor_on_intensity_of_sound(density, shear_modulus):
    return math.sqrt(shear_modulus * 1e9/(density**3))

def calculate_gruneisen_parameter(poission_ratio):
    upper_term = 3* ( 1 + poission_ratio)
    lower_term = 2 * (2 - 3 * poission_ratio)
    return upper_term / lower_term

def calculate_debye_temperature_with_cell_volume(total_number_of_ions_in_cell, cell_volume, average_velocity):
    fraction_term = 3 * total_number_of_ions_in_cell / (4 * PI * cell_volume)
    second_term = fraction_term**(1/3)
    multiplied_term = PLANCKS_CONSTANT/BOLTZMANN_CONSTANT * second_term * average_velocity
    return multiplied_term

def calculate_debye_temperature_with_molar_mass(total_number_of_ions_in_cell, molar_mass, density, average_velocity):
    first_fraction_term = 3 * total_number_of_ions_in_cell / (4 * PI)
    second_fraction_term = AVOGADRO_NUMBER * density / (molar_mass)
    fraction_term = first_fraction_term * second_fraction_term 
    second_term = fraction_term**(1/3)
    multiplied_term = PLANCKS_CONSTANT/BOLTZMANN_CONSTANT * second_term * average_velocity
    return multiplied_term

def calculate_melting_temperature(c11, c33):
    melting_temperature = 354 + 4.5/3 * (2 * c11 + c33)
    return melting_temperature

def calculate_thermal_expansion_coefficient_using_shear_modulus(shear_modulus):
    return 1.6e-3/shear_modulus
def calculate_thermal_expansion_coefficient_using_melting_temperature(melting_temp):
    return 0.02/melting_temp

def calculate_heat_capacity_per_unit_volume(total_number_of_ions_in_cell, cell_volume):
    number_of_atoms_per_unit_volume = total_number_of_ions_in_cell / cell_volume
    heat_capacity = 3 * BOLTZMANN_CONSTANT * number_of_atoms_per_unit_volume
    return heat_capacity

def calculate_minimum_thermal_conductivity_clarke(average_velocity, number_of_ions_in_cell, cell_volume):
    average_volume_per_atom = cell_volume / number_of_ions_in_cell
    return BOLTZMANN_CONSTANT*average_velocity*average_volume_per_atom**(-2/3)

def calculate_minimum_thermal_conductivity_cahill(longitudinal_velocity, transverse_velocity, number_of_ions_in_cell, cell_volume):
    last_term = (longitudinal_velocity + 2 * transverse_velocity)
    atoms_per_unit_volume = number_of_ions_in_cell/cell_volume
    return BOLTZMANN_CONSTANT/2.48 * atoms_per_unit_volume**(2/3) * last_term

def calculate_minimum_thermal_conductivity_long(possion_ratio, molar_mass, number_of_ions, density, young_modulus):
    first_part = (2 + 2 * possion_ratio)**(3/2)
    second_part = ((1/1-possion_ratio)- possion_ratio)**(3/2)
    third_part = (molar_mass/(number_of_ions*density*AVOGADRO_NUMBER))**(-2/3)
    fourth_part = (young_modulus/density)**(1/2)
    first_term = 1/3 * ( 2* first_part + second_part)
    final = first_term**(-1/3) * BOLTZMANN_CONSTANT * third_part * fourth_part
    return final

def calculate_a_gamma(grunisen_parameter):
    upper_term = 5.72e7 * 0.849 
    lower_term = 2 * (1 - 0.514/grunisen_parameter + 0.224/(grunisen_parameter)**2)
    return upper_term/lower_term

def calculate_lattice_thermal_conductivity(grunisen_parameter, debye_temperature, molar_mass, number_of_ions_in_cell, cell_volume):
    delta = pow((cell_volume/number_of_ions_in_cell), 1/3)
    M_av = molar_mass/number_of_ions_in_cell
    a_gamma = calculate_a_gamma(grunisen_parameter)
    upper_term = a_gamma * M_av * debye_temperature**3 * delta
    lower_term = grunisen_parameter**2 * number_of_ions_in_cell**(2/3) * 300
    return upper_term/lower_term

def calculate_H_micro(poission_ratio):
    upper_term = 1 - 2 * poission_ratio
    lower_term = 6 * (1 + poission_ratio)
    return upper_term / lower_term

def calculate_H_macro(bulk_modulus, shear_modulus):
    inner_term = (shear_modulus/ bulk_modulus)**2 * shear_modulus
    return 2 * inner_term**0.585 -3 

def calculate_H_v_tian(bulk_modulus, shear_modulus):
    inner_term = (shear_modulus/ bulk_modulus)**1.137 * shear_modulus**0.708
    return 0.92 * inner_term

def calculate_H_v_teter(shear_modulus):
    return 0.151 * shear_modulus

def calculate_H_v_mazhnik(poission_ratio, young_modulus):
    upper_term = 1 - 8.5 * poission_ratio + 19.5 * poission_ratio**2
    lower_term = 1 - 7.6 * poission_ratio + 12.2 * poission_ratio**2 + 19.6 * poission_ratio**3
    x_sigma = upper_term / lower_term
    return 0.096 * x_sigma * young_modulus

def calculate_wavelength_of_dominant_phonon_at_300k(average_velocity):
    upper_term = 12.566 * average_velocity * 1e-12
    return upper_term/300


    
if __name__ == "__main__":
    # print(calculate_lattice_thermal_conductivity(0.862,464.18,177.46,4,116.40e-30))
    # print('this file is running')
    print(calculate_atomic_mass("Mg8Si3P2"))
