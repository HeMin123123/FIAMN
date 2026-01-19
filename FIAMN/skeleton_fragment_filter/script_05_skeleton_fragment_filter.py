import re
from typing import List, Tuple, Dict

# Define exact atomic masses (Unit: Da)
ATOMIC_MASSES = {
    'C': 12.0000000,      # Exact atomic mass of Carbon
    'H': 1.00782503207,   # Exact atomic mass of Hydrogen
    'O': 15.99491461956   # Exact atomic mass of Oxygen (Reference, not used in fragment matching)
}

# Mass tolerance (Da)
MASS_TOLERANCE = 0.005


def calculate_ch_mass(c_count: int, h_count: int) -> float:
    """
    Calculates the exact mass of a molecular formula consisting only of C and H.
    :param c_count: Number of Carbon atoms
    :param h_count: Number of Hydrogen atoms
    :return: Exact mass
    """
    return c_count * ATOMIC_MASSES['C'] + h_count * ATOMIC_MASSES['H']


def generate_ch_formulas(max_c: int = 50, max_h: int = 100) -> Dict[float, str]:
    """
    Pre-generates formulas for C/H combinations and their theoretical masses.
    :param max_c: Maximum number of Carbon atoms
    :param max_h: Maximum number of Hydrogen atoms
    :return: Dictionary where key is theoretical mass and value is the formula string
    """
    ch_mass_formula = {}
    for c in range(1, max_c + 1):
        # Reasonable range for Hydrogen: from CnH(2n+2) (alkanes) to CnH (alkynes/radicals)
        min_h = max(1, c)  # At least 1 Hydrogen
        max_h_for_c = 2 * c + 2  # Maximum H for a saturated hydrocarbon
        for h in range(min_h, min(max_h_for_c + 1, max_h + 1)):
            mass = calculate_ch_mass(c, h)
            ch_mass_formula[mass] = f"C{c}H{h}"
    return ch_mass_formula


def parse_mgf(file_path: str) -> List[Dict]:
    """
    Parses an MGF file to extract information for each mass spectrum entry.
    :param file_path: Path to the MGF file
    :return: List of dictionaries containing entry information
    """
    mgf_entries = []
    current_entry = None
    peaks_data = []

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Start of a new entry
            if line.startswith('BEGIN IONS'):
                if current_entry is not None and peaks_data:
                    current_entry['peaks'] = peaks_data
                    mgf_entries.append(current_entry)
                current_entry = {}
                peaks_data = []
            # End of an entry
            elif line.startswith('END IONS'):
                if current_entry is not None and peaks_data:
                    current_entry['peaks'] = peaks_data
                    mgf_entries.append(current_entry)
                current_entry = None
                peaks_data = []
            # Metadata key-value pairs (e.g., TITLE, PEPMASS)
            elif '=' in line and current_entry is not None:
                key, value = line.split('=', 1)
                current_entry[key.strip()] = value.strip()
            # Peak data (m/z intensity)
            elif current_entry is not None and re.match(r'^\d+\.?\d*\s+\d+\.?\d*$', line):
                mz_str, intensity_str = line.split()
                try:
                    mz = float(mz_str)
                    intensity = float(intensity_str)
                    peaks_data.append((mz, intensity))
                except ValueError:
                    # Skip malformed peak data
                    continue

    return mgf_entries


def filter_ch_peaks(mgf_entries: List[Dict], ch_mass_formula: Dict[float, str]) -> List[Dict]:
    """
    Filters fragment peaks to keep only those matching C/H formulas.
    :param mgf_entries: Parsed original MGF entries
    :param ch_mass_formula: Dictionary of C/H formulas and their masses
    :return: Filtered MGF entries
    """
    filtered_entries = []

    for entry in mgf_entries:
        filtered_peaks = []
        # Iterate through each peak in the entry
        for mz, intensity in entry['peaks']:
            # Match against C/H formulas within mass tolerance
            matched = False
            for theory_mass, formula in ch_mass_formula.items():
                if abs(mz - theory_mass) <= MASS_TOLERANCE:
                    # Optional: Store the matched formula info if needed
                    # filtered_peaks.append((mz, intensity, formula))
                    matched = True
                    break  # Stop searching once a match is found

            if matched:
                filtered_peaks.append((mz, intensity))

        # Only keep entries that still have valid peaks after filtering
        if filtered_peaks:
            new_entry = entry.copy()
            new_entry['peaks'] = filtered_peaks
            filtered_entries.append(new_entry)

    return filtered_entries


def write_mgf(file_path: str, filtered_entries: List[Dict]):
    """
    Writes filtered results to a new MGF file.
    :param file_path: Output file path
    :param filtered_entries: List of filtered MGF entries
    """
    with open(file_path, 'w', encoding='utf-8') as f:
        for entry in filtered_entries:
            f.write("BEGIN IONS\n")
            # Write metadata
            for key, value in entry.items():
                if key != 'peaks':
                    f.write(f"{key}={value}\n")
            # Write peak data (m/z and intensity only)
            for peak in entry['peaks']:
                if isinstance(peak, tuple):
                    mz, intensity = peak[:2]
                    f.write(f"{mz:.6f} {intensity:.2f}\n")
            f.write("END IONS\n\n")


def main(input_mgf: str, output_mgf: str):
    """
    Main function: Process MGF file, filter for C/H fragment peaks, and output results.
    :param input_mgf: Path to input MGF file
    :param output_mgf: Path to output MGF file
    """
    # 1. Pre-generate C/H formula-mass library
    print("Generating C/H formula mass library...")
    ch_mass_formula = generate_ch_formulas(max_c=50, max_h=100)

    # 2. Parse input MGF file
    print(f"Parsing input file: {input_mgf}")
    mgf_entries = parse_mgf(input_mgf)
    if not mgf_entries:
        print("No valid mass spectrometry data found!")
        return

    # 3. Filter for C/H fragment peaks
    print("Filtering for C/H fragment peaks...")
    filtered_entries = filter_ch_peaks(mgf_entries, ch_mass_formula)

    # 4. Write output MGF file
    print(f"Writing to output file: {output_mgf}")
    write_mgf(output_mgf, filtered_entries)

    print(f"Processing complete! Filtered {len(filtered_entries)} valid spectral entries.")


if __name__ == "__main__":
    # Please modify the following input and output file paths
    INPUT_MGF_PATH = "precursor_formula_filter.mgf"   # Input path
    OUTPUT_MGF_PATH = "skeleton_fragment_filter.mgf" # Output path

    main(INPUT_MGF_PATH, OUTPUT_MGF_PATH)