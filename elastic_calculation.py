from elastic_formula import *
from data_extraction import *

def calculate(scf_in_path, scf_out_path, compound_name):
    scf_out_path = scf_out_path
    scf_in_path = scf_in_path
    formula = compound_name
        
    c11, c12, c13, c14, c15, c16 = extract_elastic_constants_from_file(scf_out_path)[:6]
    c21, c22, c23, c24, c25, c26 = extract_elastic_constants_from_file(scf_out_path)[6:12]
    c31, c32, c33, c34, c35, c36 = extract_elastic_constants_from_file(scf_out_path)[12:18]
    c41, c42, c43, c44, c45, c46 = extract_elastic_constants_from_file(scf_out_path)[18:24]
    c51, c52, c53, c54, c55, c56 = extract_elastic_constants_from_file(scf_out_path)[24:30]
    c61, c62, c63, c64, c65, c66 = extract_elastic_constants_from_file(scf_out_path)[30:36]
    
    # Check if other constants extraction returns valid data
    other_constants = extract_other_constants_from_file(scf_out_path)
    if isinstance(other_constants, str):
        print(f"Error extracting other constants: {other_constants}")
        return
    
    other_constants = extract_other_constants_from_file(scf_out_path)
    bulk_modulus_voigt, bulk_modulus_reuss = other_constants['bulk_modulus_voigt'], other_constants['bulk_modulus_reuss']
    bulk_modulus = other_constants['bulk_modulus']
    shear_modulus_voigt, shear_modulus_reuss = other_constants['shear_modulus_voigt'], other_constants['shear_modulus_reuss']
    shear_modulus = other_constants['shear_modulus']
    young_modulus_voigt, young_modulus_resuss = other_constants['young_modulus_voigt'], other_constants['young_modulus_reuss']
    young_modulus = other_constants['young_modulus']
    poisson_ratio_voigt, possion_ratio_resuss = other_constants['poisson_ratio_voigt'], other_constants['poisson_ratio_reuss'] 
    poisson_ratio = other_constants['poisson_ratio']

    number_of_ions = extract_number_of_ions_from_file(scf_in_path)
    number_of_species = extract_number_of_species_from_file(scf_in_path)
    formula_factor = extract_max_number_of_species_formula_unit_from_file(scf_in_path)
    final_cell_volume = extract_cell_volume_from_file(scf_out_path)

    molar_mass = calculate_atomic_mass(formula) / 1000

    lable_width = 47
    value_width = 12
    line_width = 71
    unit_width = 8
    print("=" * line_width)
    print(f"{'Elastic Properties':^71}")
    print("=" * line_width)
    print(f"{'Final cell volume':<{lable_width}} : {final_cell_volume:>{value_width}.5e} {'m3':>{unit_width}}")
    print(f"{'Bulk modulus':<{lable_width}} : {bulk_modulus:>{value_width}.5f} {'GPa':>{unit_width}}")
    print(f"{'Shear modulus':<{lable_width}} : {shear_modulus:>{value_width}.5f} {'GPa':>{unit_width}}")
    print(f"{'Young modulus':<{lable_width}} : {young_modulus:>{value_width}.5f} {'GPa':>{unit_width}}")
    print(f"{'Poisson ratio':<{lable_width}} : {poisson_ratio:>{value_width}.5f}")
    
    Pughs_Ratio = calculate_pughs_ratio(bulk_modulus, shear_modulus)
    print(f"{'Pughs Ratio':<{lable_width}} : {Pughs_Ratio:>{value_width}.5f}")

    tryetragonal_Shear_modulus = calculate_tetragonal_shear_modulus(c11, c12)
    print(f"{'Tetragonal Shear modulus':<{lable_width}} : {tryetragonal_Shear_modulus:>{value_width}.5f} {'GPa':>{unit_width}}")

    cauchy_pressure = calculate_cauchy_pressure(c12, c44)
    print(f"{'Cauchy pressure':<{lable_width}} : {cauchy_pressure:>{value_width}.5f} {'GPa':>{unit_width}}")

    kleixman_parameter = calculate_kleixman_parameter(c11, c12)
    print(f"{'Kleixman parameter':<{lable_width}} : {kleixman_parameter:>{value_width}.5f}")

    machinability = calculate_machinability(bulk_modulus, c44)
    print(f"{'Machinability':<{lable_width}} : {machinability:>{value_width}.5f}")

    print(f"{'(c11 - c12)':<{lable_width}} : {c11 - c12:>{value_width}.5f} {'GPa':>{unit_width}}")

    compressibility = 1. / bulk_modulus
    print(f"{'Compressibility':<{lable_width}} : {compressibility:>{value_width}.5f} {'1/GPa':>{unit_width}}")

    elastic_debye_temperature = extract_elastic_debye_temperature_from_file(scf_out_path)
    print(f"{'Elastic Debye temperature':<{lable_width}} : {elastic_debye_temperature:>{value_width}.5f} {'K':>{unit_width}}")

    averaged_debye_sound_velocity = extract_average_debye_sound_velocity_from_file(scf_out_path)
    print(f"{'Average Debye sound velocity':<{lable_width}} : {averaged_debye_sound_velocity:>{value_width}.5f} {'m/s':>{unit_width}}")
    print()

    print("=" * line_width)
    print(f"{'Elastic Anisotropy (Hexagonal)':^71}")
    print("=" * line_width)
    a1 = calculate_a1(c11, c13, c33, c44)
    print(f"{'A_1':<{lable_width}} : {a1:>{value_width}.5f}")
    
    a2 = calculate_a2(c22, c23, c33, c55)
    print(f"{'A_2':<{lable_width}} : {a2:>{value_width}.5f}")
    
    a3 = calculate_a3(a1, a2)
    print(f"{'A_3':<{lable_width}} : {a3:>{value_width}.5f}")
    
    zener_anisotropy_factor = calculate_zener_anisotropy_factor(c11, c12, c44)
    print(f"{'Zener anisotropy factor (A)':<{lable_width}} : {zener_anisotropy_factor:>{value_width}.5f}")

    c44_reuss = calculate_c44_reuss(c11, c12, c44)
    print(f"{'C44 Reuss':<{lable_width}} : {c44_reuss:>{value_width}.5f}")

    c44_voigt = calculate_c44_voigt(c44_reuss, c11, c12, c44)
    print(f"{'C44 Voigt':<{lable_width}} : {c44_voigt:>{value_width}.5f}")

    universal_log_euclidean_index = calculate_universal_log_euclidean_index(c44_voigt, c44_reuss, bulk_modulus_voigt, bulk_modulus_reuss)
    print(f"{'Universal log euclidean index (A_L)':<{lable_width}} : {universal_log_euclidean_index:>{value_width}.5f}")

    universal_anisotropy_index = calculate_universal_anisotropy_index(bulk_modulus_voigt, bulk_modulus_reuss, shear_modulus_voigt, shear_modulus_reuss)
    print(f"{'Universal anisotropy index (A_u)':<{lable_width}} : {universal_anisotropy_index:>{value_width}.5f} >= 0")

    anisotropy_in_compressibility = calculate_anisotropy_in_compressibility(bulk_modulus_voigt, bulk_modulus_reuss)
    print(f"{'Anisotropy in compressibility (A_B)':<{lable_width}} : {anisotropy_in_compressibility:>{value_width}.5} %")

    anisotropy_in_shear = calculate_anisotropy_in_shear(shear_modulus_voigt, shear_modulus_reuss, shear_modulus)
    print(f"{'Anisotropy in shear (A_G)':<{lable_width}} : {anisotropy_in_shear:>{value_width}.5f}")

    equivalent_zener_anisotropy = calculate_equivalent_zener_anisotropy(universal_log_euclidean_index)
    print(f"{'Equivalent zener anisotropy (A_eq)':<{lable_width}} : {equivalent_zener_anisotropy:>{value_width}.5f}")
    print()

    print("=" * line_width)
    print(f"{'Universal Bulk modulus along a, b & c axis':^71}")
    print("=" * line_width)
    uniaxial_bulk_modulus_Ba = calculate_uniaxial_bulk_modulus_Ba(c11, c12, c13, c22, c23, c33)
    print(f"{'Uniaxial bulk modulus Ba':<{lable_width}} : {uniaxial_bulk_modulus_Ba:>{value_width}.5f}")

    uniaxial_bulk_modulus_Bb = calculate_uniaxial_bulk_modulus_Bb(c11, c12, c13, c22, c23, c33)
    print(f"{'Uniaxial bulk modulus Bb':<{lable_width}} : {uniaxial_bulk_modulus_Bb:>{value_width}.5f}")

    uniaxial_bulk_modulus_Bc = calculate_uniaxial_bulk_modulus_Bc(c11, c12, c13, c22, c23, c33)
    print(f"{'Uniaxial bulk modulus Bc':<{lable_width}} : {uniaxial_bulk_modulus_Bc:>{value_width}.5f}")

    uniaxial_alpha = calculate_uniaxial_alpha(c11, c12, c13, c22, c23, c33)
    print(f"{'α':<{lable_width}} : {uniaxial_alpha:>{value_width}.5f}")

    uniaxial_beta = calculate_uniaxial_beta(c11, c12, c13, c22, c23, c33)
    print(f"{'β':<{lable_width}} : {uniaxial_beta:>{value_width}.5f}")

    uniaxial_lambda = calculate_uniaxial_lambda(c11, c12, c13, c22, c33, uniaxial_alpha, uniaxial_beta)
    print(f"{'λ':<{lable_width}} : {uniaxial_lambda:>{value_width}.5f}")

    uniaxial_beta_relax = calculate_uniaxial_beta_relax(c11, c12, c13, c22, c23, c33)
    print(f"{'β_relax':<{lable_width}} : {uniaxial_beta_relax:>{value_width}.5f}")
    print()


    print("=" * line_width)
    print(f"{'Anisotropic Bulk modulus along a & c axis':^71}")
    print("=" * line_width)
    anisotropic_bulk_modulus_ABa = calculate_anisotropic_bulk_modulus_ABa(uniaxial_bulk_modulus_Ba, uniaxial_bulk_modulus_Bc)
    print(f"{'Anisotropic bulk modulus ABa':<{lable_width}} : {anisotropic_bulk_modulus_ABa:>{value_width}.5f}")

    anisotropic_bulk_modulus_ABc = calculate_anisotropic_bulk_modulus_ABc(uniaxial_bulk_modulus_Ba, uniaxial_bulk_modulus_Bc)
    print(f"{'Anisotropic bulk modulus ABb':<{lable_width}} : {anisotropic_bulk_modulus_ABc:>{value_width}.5f}")
    print()

    print("=" * line_width)
    print(f"{'Linear compressibility along a & c axis':^71}")
    print("=" * line_width)
    linear_compressibility_beta_a = calculate_linear_compressibility_beta_a(c11, c12, c13, c33)
    print(f"{'Linear compressibility βa':<{lable_width}} : {linear_compressibility_beta_a:>{value_width}.5e}")

    linear_compressibility_beta_c = calculate_linear_compressibility_beta_c(c11, c12, c13, c33)
    print(f"{'Linear compressibility βc':<{lable_width}} : {linear_compressibility_beta_c:>{value_width}.5e}")

    ratio_of_linear_compressibility = calculate_ratio_of_linear_compressibility(linear_compressibility_beta_a, linear_compressibility_beta_c)
    print(f"{'Ratio of linear compressibility':<{lable_width}} : {ratio_of_linear_compressibility:>{value_width}.5f}")
    print()

    print("=" * line_width)
    print(f"{'Hardness (GPa)':^71}")
    print("=" * line_width)
    H_micro = calculate_H_micro(poisson_ratio)
    print(f"{'H_micro':<{lable_width}} : {H_micro:>{value_width}.5f}")

    H_macro = calculate_H_macro(bulk_modulus, shear_modulus)
    print(f"{'H_macro':<{lable_width}} : {H_macro:>{value_width}.5f}")

    H_v_tian = calculate_H_v_tian(bulk_modulus, shear_modulus)
    print(f"{'(Hv) tian':<{lable_width}} : {H_v_tian:>{value_width}.5f}")

    H_v_teter = calculate_H_v_teter(shear_modulus)
    print(f"{'(Hv) teter':<{lable_width}} : {H_v_teter:>{value_width}.5f}")

    H_v_mazhnik = calculate_H_v_mazhnik(poisson_ratio, young_modulus)
    print(f"{'(Hv) mazhnik':<{lable_width}} : {H_v_mazhnik:>{value_width}.5f}")
    print()

    print("=" * line_width)
    print(f"{'Acoustic Properties':^71}")
    print("=" * line_width)
    density_calculated = calculate_density(formula_factor, molar_mass, final_cell_volume)
    print(f"{'Density (calculated)':<{lable_width}} : {density_calculated:>{value_width}.5f} {'kg/m3':>{unit_width}}")

    # density = extract_density_from_file(castep_file_path)
    # if density is not None:
    #     print(f"{'Density (from file)':<{lable_width}} : {density:>{value_width}.5f} {'kg/m3':>{unit_width}}")
    # else:
    #     print("**warning: You need to have later version of material studio to extract density from file**")

    longitudinal_velocity = calculate_longitudinal_velocity(bulk_modulus, shear_modulus, density_calculated)
    print(f"{'Longitudinal velocity':<{lable_width}} : {longitudinal_velocity:>{value_width}.5f} {'m/s':>{unit_width}}")

    transverse_velocity = calculate_transverse_velocity(shear_modulus, density_calculated)
    print(f"{'Transverse velocity':<{lable_width}} : {transverse_velocity:>{value_width}.5f} {'m/s':>{unit_width}}")

    averaged_velocity_calculated = calculate_average_velocity(longitudinal_velocity, transverse_velocity)
    print(f"{'Averaged velocity (calculated)':<{lable_width}} : {averaged_velocity_calculated:>{value_width}.5f} {'m/s':>{unit_width}}")

    acoustic_impedance = calculate_acoustic_impedance(density_calculated, shear_modulus)
    print(f"{'Acoustic impedance':<{lable_width}} : {acoustic_impedance:>{value_width}.5e} {'rayl':>{unit_width}}")

    radiation_factor = calculate_radiation_factor_on_intensity_of_sound(density_calculated, shear_modulus)
    print(f"{'Radiation factor on intensity of sound':<{lable_width}} : {radiation_factor:>{value_width}.5f}")
    print()

    print("=" * line_width)
    print(f"{'Thermo Physical Properties':^71}")
    print("=" * line_width)
    debye_temperature_with_cell_volume = calculate_debye_temperature_with_cell_volume(number_of_ions, final_cell_volume, averaged_velocity_calculated)
    print(f"{'Debye temperature with cell volume':<{lable_width}} : {debye_temperature_with_cell_volume:>{value_width}.5f} {'k':>{unit_width}}")

    debye_temperature_with_molar_mass = calculate_debye_temperature_with_molar_mass(number_of_ions, molar_mass, density_calculated, averaged_velocity_calculated)
    print(f"{'Debye temperature with molar mass':<{lable_width}} : {debye_temperature_with_molar_mass:>{value_width}.5f} {'K':>{unit_width}}")

    melting_temperature = calculate_melting_temperature(c11, c33)
    print(f"{'Melting temperature':<{lable_width}} : {melting_temperature:>{value_width}.5f} +- 300K")

    thermal_expansion_coefficient_using_melting_temperature = calculate_thermal_expansion_coefficient_using_melting_temperature(melting_temperature)
    print(f"{'Thermal expansion coeff. with Tm':<{lable_width}} : {thermal_expansion_coefficient_using_melting_temperature:>{value_width}.5e}")

    thermal_expansion_coefficient_using_shear_modulus = calculate_thermal_expansion_coefficient_using_shear_modulus(shear_modulus)
    print(f"{'Thermal expansion coeff. with G':<{lable_width}} : {thermal_expansion_coefficient_using_shear_modulus:>{value_width}.5e}")

    heat_capacity = calculate_heat_capacity_per_unit_volume(number_of_ions, final_cell_volume)
    print(f"{'Heat capacity per unit volume':<{lable_width}} : {heat_capacity:>{value_width}.5e} {'J/molK':>{unit_width}}")

    minimum_thermal_conductivity_cahill = calculate_minimum_thermal_conductivity_cahill(longitudinal_velocity, transverse_velocity, number_of_ions, final_cell_volume)
    print(f"{'Minimum thermal conductivity Cahill':<{lable_width}} : {minimum_thermal_conductivity_cahill:>{value_width}.5e} {'W/mK':>{unit_width}}")

    minimum_thermal_conductivity_clarke = calculate_minimum_thermal_conductivity_clarke(averaged_velocity_calculated, number_of_ions, final_cell_volume)
    print(f"{'Minimum thermal conductivity Clarke':<{lable_width}} : {minimum_thermal_conductivity_clarke:>{value_width}.5e} {'W/mK':>{unit_width}}")

    minimum_thermal_conductivity_long = calculate_minimum_thermal_conductivity_long(poisson_ratio,molar_mass,number_of_ions,density_calculated, young_modulus)
    print(f"{'Minimum thermal conductivity Long (K_min_long)':<{lable_width}} : {minimum_thermal_conductivity_long:>{value_width}.5e} {'W/mK':>{unit_width}}")

    gruneisen_parameter = calculate_gruneisen_parameter(poisson_ratio)
    print(f"{'Gruneisen parameter':<{lable_width}} : {gruneisen_parameter:>{value_width}.5f}")

    debye_temperature_from_file = extract_elastic_debye_temperature_from_file(scf_out_path)
    print(f"{'Debye temperature (from file)':<{lable_width}} : {debye_temperature_from_file:>{value_width}.5f} {'K':>{unit_width}}")

    lattice_thermal_conductivity = calculate_lattice_thermal_conductivity(gruneisen_parameter, debye_temperature_from_file, molar_mass, number_of_ions, final_cell_volume)
    print(f"{'Lattice thermal conductivity (kph)':<{lable_width}} : {lattice_thermal_conductivity:>{value_width}.5e} {'W/mK':>{unit_width}}")

    wavelength_of_dominant_phonon_at_300k = calculate_wavelength_of_dominant_phonon_at_300k(averaged_velocity_calculated)
    print(f"{'Wavelength of dominant phonon at 300K':<{lable_width}} : {wavelength_of_dominant_phonon_at_300k:>{value_width}.5e} {'m':>{unit_width}}")

if __name__ == "__main__":
    scf_in_path = "MgB2.Elastic/MgB2.scf.in"
    scf_out_path = "MgB2.Elastic/MgB2.scf.out"
    compound_name = "MgB2"

    calculate(scf_in_path, scf_out_path, compound_name)
