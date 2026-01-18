import os

def filter_mgf_with_nitrogen_rule(input_file, output_file):
    """
    Filters an MGF file based on specific chemical rules:
    1. Carbon count (C) between 18 and 25.
    2. Elements restricted to C, H, and O only.
    3. Mass tolerance of 0.005 Da.
    4. Ionization forms: [M+H]+ or [M+Na]+.
    5. Nitrogen Rule: For pure CHO neutral molecules, the number of H must be even.
    """
    # Exact atomic masses
    MASS_C = 12.000000
    MASS_H = 1.007825
    MASS_O = 15.994915
    MASS_PROTON = 1.007276  # Mass of H+
    MASS_NA = 22.989218      # Mass of Na+

    TOLERANCE = 0.005
    C_RANGE = range(18, 26)  # 18 to 25 inclusive

    match_count = 0

    def is_valid_cho_formula(pepmass):
        # Iterate through possible ionization adducts
        for adduct_mass in [MASS_PROTON, MASS_NA]:
            neutral_mass = pepmass - adduct_mass
            if neutral_mass <= 0:
                continue

            # Iterate through carbon range
            for c in C_RANGE:
                c_mass = c * MASS_C
                if c_mass > neutral_mass + TOLERANCE:
                    continue

                # Iterate through possible oxygen counts
                max_o = int((neutral_mass - c_mass + TOLERANCE) / MASS_O)
                for o in range(0, max_o + 1):
                    rem_h_mass = neutral_mass - (c_mass + o * MASS_O)

                    # Calculate required H count (rounded)
                    h = round(rem_h_mass / MASS_H)

                    # --- Core Validation Logic ---
                    # 1. H count must be non-negative
                    # 2. Nitrogen Rule: Neutral CHO molecules must have an even number of hydrogens
                    # 3. Valence Limit: H cannot exceed (2C + 2) (saturated alkane limit)
                    if h >= 0 and h % 2 == 0 and h <= (2 * c + 2):
                        calc_mass = c_mass + h * MASS_H + o * MASS_O
                        if abs(calc_mass - neutral_mass) <= TOLERANCE:
                            return True
        return False

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        return

    print(f"Processing: {input_file}...")

    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        current_ion = []
        keep_ion = False

        for line in f:
            current_ion.append(line)
            if line.startswith('PEPMASS='):
                try:
                    # Extract precursor mass
                    mz = float(line.split('=')[1].strip())
                    if is_valid_cho_formula(mz):
                        keep_ion = True
                except (ValueError, IndexError):
                    pass

            if line.strip() == 'END IONS':
                if keep_ion:
                    out.writelines(current_ion)
                    match_count += 1

                # Reset state for next ion
                current_ion = []
                keep_ion = False

    print("-" * 30)
    print(f"Processing Complete!")
    print(f"Spectra matching CHO rules (even H): {match_count}")
    print(f"Filtered results saved to: {output_file}")
    print("-" * 30)


# Execute the script
if __name__ == "__main__":
    input_filename = 'family_network.mgf'  # Ensure this file exists
    output_filename = 'precursor_formula_filter.mgf'
    filter_mgf_with_nitrogen_rule(input_filename, output_filename)