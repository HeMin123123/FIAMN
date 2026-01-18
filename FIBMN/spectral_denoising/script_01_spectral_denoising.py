import numpy as np
from scipy import stats
import os


import os
import numpy as np
from scipy import stats

def filter_mgf_by_nonzero_mode(input_path, output_path):
    """
    Filters MGF file spectra by removing peaks with intensity equal to 0
    or equal to the most common non-zero intensity value (noise floor).
    """
    if not os.path.exists(input_path):
        print(f"Error: Input file '{input_path}' not found.")
        return

    print(f"Processing: {input_path} ...")

    with open(input_path, 'r') as f:
        lines = f.readlines()

    filtered_mgf = []
    spectrum_data = []
    in_spectrum = False
    header_info = []

    total_ions = 0
    removed_ions = 0

    for line in lines:
        line = line.strip()
        if not line: continue

        if line == "BEGIN IONS":
            in_spectrum = True
            spectrum_data = []
            header_info = []
            filtered_mgf.append(line)
            continue

        if line == "END IONS":
            in_spectrum = False
            if spectrum_data:
                # 1. Convert to numpy arrays
                mz_array = np.array([x[0] for x in spectrum_data])
                intensity_array = np.array([x[1] for x in spectrum_data])

                # 2. Pre-filter: Remove peaks with 0 intensity
                valid_mask = intensity_array > 0
                mz_v = mz_array[valid_mask]
                int_v = intensity_array[valid_mask]

                if len(int_v) > 0:
                    # 3. Calculate mode of non-zero intensities
                    # (Rounding helps identify the baseline noise level)
                    rounded_int = np.round(int_v)
                    mode_result = stats.mode(rounded_int, keepdims=True)
                    nonzero_mode_val = mode_result.mode[0]

                    # 4. Final filter: Intensity must be > the non-zero mode
                    final_mask = int_v > nonzero_mode_val
                    final_mz = mz_v[final_mask]
                    final_int = int_v[final_mask]
                else:
                    final_mz = []
                    final_int = []

                # Statistics tracking
                total_ions += len(intensity_array)
                removed_ions += (len(intensity_array) - len(final_int))

                # 5. Write headers and filtered ions back to list
                filtered_mgf.extend(header_info)
                for m, i in zip(final_mz, final_int):
                    filtered_mgf.append(f"{m:.5f} {i:.1f}")

            filtered_mgf.append("END IONS")
            continue

        if in_spectrum:
            if "=" in line:
                header_info.append(line)
            else:
                try:
                    parts = line.split()
                    if len(parts) >= 2:
                        spectrum_data.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
        else:
            filtered_mgf.append(line)

    with open(output_path, 'w') as f:
        f.write("\n".join(filtered_mgf) + "\n")

    print(f"--- Processing Complete ---")
    print(f"Output File: {output_path}")
    print(f"Total Ions: {total_ions}")
    print(f"Removed Ions (0s and Mode): {removed_ions}")
    if total_ions > 0:
        print(f"Removal Rate: {removed_ions / total_ions:.1%}")


if __name__ == '__main__':
    # ==========================================
    # Configuration
    # ==========================================
    input_file = "Text.mgf"       # Replace with your filename
    output_file = "denoising.mgf"  # Output filename

    filter_mgf_by_nonzero_mode(input_file, output_file)
    pass

