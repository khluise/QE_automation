def extract_elastic_constants_from_file(file_path):
    pattern = "Elastic constants C_ij (kbar)"
    with open(file_path, 'r') as file:
        lines = file.readlines()
        start_line = None

        for i, line in enumerate(lines):
            line = line.strip()
            if pattern in line:
                start_line = i + 2
                break
        if start_line is not None:
            extracted_data = lines[start_line:start_line + 6]
            data_list = []
            for line in extracted_data:
                data_list.extend(float(number) for number in line.split())
            
            if len(data_list) == 42:
                c11, c12, c13 = data_list[1], data_list[2], data_list[3]
                c14, c15, c16 = data_list[4], data_list[5], data_list[6]
                c21, c22, c23 = data_list[8], data_list[9], data_list[10]
                c24, c25, c26 = data_list[11], data_list[12], data_list[13]
                c31, c32, c33 = data_list[15], data_list[16], data_list[17]
                c34, c35, c36 = data_list[18], data_list[19], data_list[20]
                c41, c42, c43 = data_list[22], data_list[23], data_list[24]
                c44, c45, c46 = data_list[25], data_list[26], data_list[27]
                c51, c52, c53 = data_list[29], data_list[30], data_list[31]
                c54, c55, c56 = data_list[32], data_list[33], data_list[34]
                c61, c62, c63 = data_list[36], data_list[37], data_list[38]
                c64, c65, c66 = data_list[39], data_list[40], data_list[41]
                
                #convert kbar to GPa
                c11, c12, c13 = c11 * 0.1, c12 * 0.1, c13 * 0.1
                c14, c15, c16 = c14 * 0.1, c15 * 0.1, c16 * 0.1
                c21, c22, c23 = c21 * 0.1, c22 * 0.1, c23 * 0.1
                c24, c25, c26 = c24 * 0.1, c25 * 0.1, c26 * 0.1
                c31, c32, c33 = c31 * 0.1, c32 * 0.1, c33 * 0.1
                c34, c35, c36 = c34 * 0.1, c35 * 0.1, c36 * 0.1
                c41, c42, c43 = c41 * 0.1, c42 * 0.1, c43 * 0.1
                c44, c45, c46 = c44 * 0.1, c45 * 0.1, c46 * 0.1
                c51, c52, c53 = c51 * 0.1, c52 * 0.1, c53 * 0.1
                c54, c55, c56 = c54 * 0.1, c55 * 0.1, c56 * 0.1
                c61, c62, c63 = c61 * 0.1, c62 * 0.1, c63 * 0.1
                c64, c65, c66 = c64 * 0.1, c65 * 0.1, c66 * 0.1
                return (
                    c11, c12, c13, c14, c15, c16,
                    c21, c22, c23, c24, c25, c26,
                    c31, c32, c33, c34, c35, c36,
                    c41, c42, c43, c44, c45, c46,
                    c51, c52, c53, c54, c55, c56,
                    c61, c62, c63, c64, c65, c66
                )

            else:
                return "Data list and constants list lengths do not match."
        else:
            return "'Elastic constants C_ij (kbar)' not found in the file."

def extract_other_constants_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        start_line = None
        pattern_1 = "Voigt approximation:"
        for i, line in enumerate(lines):
                line = line.strip()
                if pattern_1 in line:
                    start_line = i + 1
                    break
        if start_line is not None:
            extracted_data = lines[start_line:start_line + 20]
            data_list = []
            for line in extracted_data:
                words = line.split('=')
                if len(words) == 2:
                    value = words[1].strip().split()[0]
                    try:
                        data_list.append(float(value))
                    except ValueError:
                        continue
            # bulk_modulus_voigt, bulk_modulus_reuss, bulk_modulus = data_list[0], data_list[1], data_list[2]
            # shear_modulus_voigt, shear_modulus_reuss, shear_modulus = data_list[3], data_list[4], data_list[5]
            # lame_lambda_voigt, lame_lambda_reuss, lame_lambda = data_list[6], data_list[7], data_list[8]
            # young_modulus_voigt, young_modulus_resuss, young_modulus = data[9], data[10], data[11]
            # poisson_ratio_voigt, possion_ratio_resuss, poisson_ratio = data[12], data[13], data[14]
            constant_dictonary = {'bulk_modulus_voigt': data_list[0] * 0.1, 
                                   'bulk_modulus_reuss': data_list[5] * 0.1,
                                   'bulk_modulus': data_list[10] * 0.1,
                                   'shear_modulus_voigt': data_list[2] * 0.1, 
                                   'shear_modulus_reuss': data_list[7] * 0.1, 
                                   'shear_modulus': data_list[12] * 0.1,
                                   'young_modulus_voigt': data_list[1] * 0.1, 
                                   'young_modulus_reuss': data_list[5] * 0.1,
                                   'young_modulus' : data_list[11] * 0.1, 
                                   'poisson_ratio_voigt': data_list[3], 
                                   'poisson_ratio_reuss': data_list[8],
                                   'poisson_ratio': data_list[14]}
            return constant_dictonary
        else:
            return "Error in extracting other constants"

def extract_cell_volume_from_file(file_path):
    pattern = "unit-cell volume"
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            cell_volume = float(last_match.split('=')[1].strip().split(' ')[0])
            return cell_volume * 1.48185e-31
        else:
            print("No match found")

def extract_density_from_file(file_path):
    pattern = 'density ='
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            density = float(last_match.split('=')[1].strip().split(' ')[0])
            return density * 1660.54
        else:
            print("No match found")

def extract_elastic_debye_temperature_from_file(file_path):
    pattern = 'Debye temperature ='
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            elastic_debye_temperature = float(last_match.split('=')[1].strip().split(' ')[0])
            return elastic_debye_temperature
        else:
            print("No match found")

def extract_average_debye_sound_velocity_from_file(file_path):
    pattern = 'Average Debye sound velocity'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            averaged_sound_velocity = float(last_match.split('=')[1].strip().split(' ')[0])
            return averaged_sound_velocity
        else:
            print("No match found")

def extract_number_of_ions_from_file(file_path):
    pattern = 'nat'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            # Check if '=' exists in the line
            if '=' in last_match:
                split_parts = last_match.split('=')
                if len(split_parts) > 1:
                    number_of_ions = float(split_parts[1].strip())
                    return number_of_ions
                else:
                    print("No value found after '=' in line:", last_match)
                    return None
            else:
                print("No '=' found in line:", last_match)
                return None
        else:
            print("No match found")
            return None

def extract_number_of_species_from_file(file_path):
    pattern = 'ntyp'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            # Check if '=' exists in the line
            if '=' in last_match:
                split_parts = last_match.split('=')
                if len(split_parts) > 1:
                    number_of_species = float(split_parts[1].strip())
                    return number_of_species
                else:
                    print("No value found after '=' in line:", last_match)
                    return None
            else:
                print("No '=' found in line:", last_match)
                return None
        else:
            print("No match found")
            return None
            
def extract_max_number_of_species_formula_unit_from_file(file_path):
    pattern_start = 'ATOMIC_POSITIONS'
    pattern_end = 'K_POINTS'
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    start_line = None
    end_line = None
    
    # Find start and end lines
    for i, line in enumerate(lines):
        if pattern_start in line:
            start_line = i + 1  # Skip the ATOMIC_POSITIONS line itself
        elif pattern_end in line and start_line is not None:
            end_line = i
            break
    
    if start_line is not None and end_line is not None:
        atomic_lines = lines[start_line:end_line]
        species_count = {}
        
        # Count each species
        for line in atomic_lines:
            line = line.strip()
            if line:  # Skip empty lines
                species = line.split()[0]  # First column is the species
                species_count[species] = species_count.get(species, 0) + 1
        
        if species_count:
            max_count = max(species_count.values())
            return max_count
        else:
            print("No atomic positions found")
            return None
    else:
        print("ATOMIC_POSITIONS or K_POINTS section not found")
        return None
def print_non_zero_value_of_Cij(elastic_constants):
    count = 0
    for i in range(6):
        for j in range(6):
            if elastic_constants[count] != 0:
                print(f"c{i+1}{j+1} {' ' * (20 - len(f'c{i+1}{j+1}'))}: {elastic_constants[count]}")
            count+=1

if __name__ == "__main__":
    scf_out_file = "MgB2.Elastic/MgB2.scf.out"
    scf_in_file = "MgB2.Elastic/MgB2.scf.in"
    elastic_constants = extract_elastic_constants_from_file(scf_out_file)
    if isinstance(elastic_constants, tuple):
        print_non_zero_value_of_Cij(elastic_constants)
    else:
        print(elastic_constants)
        
    other_constants = extract_other_constants_from_file(scf_out_file)
    if isinstance(other_constants, dict):
        for key, value in other_constants.items():
            print(f"{key}: {type(value)}")

    cell_volume = extract_cell_volume_from_file(scf_out_file)
    if cell_volume is not None:
        print(f"Cell volume: {cell_volume} m^3")
        
    elastic_debye_temperature = extract_elastic_debye_temperature_from_file(scf_out_file)
    if elastic_debye_temperature is not None:
        print(f"Elastic Debye temperature: {elastic_debye_temperature} K")
    number_of_ions = extract_number_of_ions_from_file(scf_in_file)
    if number_of_ions is not None:
        print(f"Number of ions: {number_of_ions}")
    number_of_species = extract_number_of_species_from_file(scf_in_file)
    if number_of_species is not None:
        print(f"Number of species: {number_of_species}")
    max_number_of_species_formula_unit = extract_max_number_of_species_formula_unit_from_file(scf_in_file)
    if max_number_of_species_formula_unit is not None:      
        print(f"Max number of species in formula unit: {max_number_of_species_formula_unit}")
