def filter_mgf_by_total_intensity(input_file, output_file, min_total_intensity=2000):
    """
    Filters spectra in an MGF file whose total ion intensity is below a specified threshold.

    Args:
        input_file (str): Path to the input MGF file.
        output_file (str): Path to the output MGF file.
        min_total_intensity (float): Minimum total ion intensity threshold, default is 2000.
    """
    # Initialize variables
    current_spectrum = []  # Stores all lines of the current spectrum
    in_spectrum = False    # Flag to indicate if currently reading spectrum content
    total_intensity = 0.0  # Total ion intensity of the current spectrum

    try:
        with open(input_file, 'r', encoding='utf-8') as infile, \
                open(output_file, 'w', encoding='utf-8') as outfile:

            # Process the file line by line
            for line in infile:
                stripped_line = line.strip()

                # Start processing the spectrum when 'BEGIN IONS' is encountered
                if stripped_line == 'BEGIN IONS':
                    in_spectrum = True
                    current_spectrum = [line]
                    total_intensity = 0.0
                    continue

                # End processing, calculate total intensity, and decide whether to save
                elif stripped_line == 'END IONS':
                    in_spectrum = False
                    current_spectrum.append(line)

                    # Write to file only if total ion intensity meets the threshold
                    if total_intensity >= min_total_intensity:
                        outfile.writelines(current_spectrum)

                    continue

                # Handle content inside the spectrum block
                if in_spectrum:
                    current_spectrum.append(line)

                    # Parse peak data (Format: m/z intensity [additional info])
                    parts = stripped_line.split()
                    if len(parts) >= 2:
                        try:
                            # The second column is typically the intensity
                            intensity = float(parts[1])
                            total_intensity += intensity
                        except ValueError:
                            # Skip non-numeric lines (e.g., metadata/header lines like TITLE=)
                            continue

            print(f"Processing complete!")
            print(f"Input file: {input_file}")
            print(f"Output file: {output_file}")
            print(f"Filtering threshold: Total Ion Intensity >= {min_total_intensity}")

    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
    except PermissionError:
        print(f"Error: Permission denied for {input_file} or {output_file}.")
    except Exception as e:
        print(f"An error occurred during processing: {str(e)}")


# Example usage
if __name__ == "__main__":
    # Update these paths to your actual file locations
    input_mgf_path = "denoising.mgf"      # Input MGF file
    output_mgf_path = "filtered_TIC.mgf"  # Output filtered MGF file

    # Execute the filtering process
    filter_mgf_by_total_intensity(input_mgf_path, output_mgf_path, min_total_intensity=2000)