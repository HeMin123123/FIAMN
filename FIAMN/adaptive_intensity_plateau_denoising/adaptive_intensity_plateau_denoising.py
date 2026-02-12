import os
import numpy as np
from tqdm import tqdm

# Note: These modules must exist in your local environment
from finger_id.g_config import config, update_config_on_main
from finger_id.mgf import get_mgf, write_mgf
from finger_id.mgf_similarity import id_of_spectrum, copy_of_spectrum, np_array_of_spectrum


def adaptive_plateau_denoising(spectrum):
    """
    Performs adaptive denoising based on the longest plateau in the sorted intensity derivative.
    """
    spectrum = copy_of_spectrum(spectrum)
    mz_ori = spectrum['mz_array']
    intensity_ori = spectrum['intensity_array']

    # Pair m/z and intensity, filter out zero-intensity peaks
    mz_itst = list(zip(mz_ori, intensity_ori))
    mz_itst = list(filter(lambda mi: mi[1] > 0, mz_itst))
    
    # Sort peaks by intensity in ascending order
    mz_itst.sort(key=lambda mi: mi[1])
    print('Peak Count:', 'Original:', len(intensity_ori), 'Non-zero:', len(mz_itst))
    
    intensities = np.array([mi[1] for mi in mz_itst])

    # Calculate the first derivative (difference between adjacent sorted intensities)
    shifted_intensities = np.concatenate(((intensities[0],), intensities))[:-1]
    intensity_derivative = intensities - shifted_intensities
    
    print('Derivative non-zero count:', len([val for val in intensity_derivative if val > 0]))

    # Identify indices where the derivative is near zero (the plateau/flat regions)
    plateau_start_indices = []
    plateau_end_indices = []
    
    for i in range(1, len(intensity_derivative) - 1):
        if abs(intensity_derivative[i]) < 1e-6:
            # Start of a plateau
            if abs(intensity_derivative[i-1]) > 1e-6:
                plateau_start_indices.append(i)
            # End of a plateau
            if abs(intensity_derivative[i+1]) > 1e-6:
                plateau_end_indices.append(i)

    plateau_start_indices = np.array(plateau_start_indices)
    plateau_end_indices = np.array(plateau_end_indices)

    # Align start and end indices to ensure valid pairs
    if len(plateau_start_indices) > 0 and len(plateau_end_indices) > 0:
        if plateau_start_indices[0] > plateau_end_indices[0]:
            plateau_end_indices = plateau_end_indices[1:]
            plateau_start_indices = plateau_start_indices[0:len(plateau_end_indices)]
        
        plateau_lengths = plateau_end_indices - plateau_start_indices
    else:
        # Fallback if no plateau is found
        return 0, spectrum

    # Find the longest plateau (the region with the most stable noise)
    longest_plateau_idx = np.argmax(plateau_lengths)
    end_idx_of_longest_plateau = plateau_end_indices[longest_plateau_idx]
    
    # The threshold is defined as the intensity at the end of the longest plateau
    dynamic_threshold = intensities[end_idx_of_longest_plateau]

    # Filter: retain peaks where intensity > dynamic_threshold
    mz_arr, intensity_arr = np_array_of_spectrum(spectrum)
    valid_indices = np.where(intensity_arr > dynamic_threshold)
    
    filtered_intensity = intensity_arr[valid_indices]
    filtered_mz = mz_arr[valid_indices]

    # Store old data and update spectrum with denoised data
    spectrum['mz_array_old_filter_mi'] = spectrum['mz_array']
    spectrum['intensity_array_old_filter_mi'] = spectrum['intensity_array']
    spectrum['mz_array'] = filtered_mz
    spectrum['intensity_array'] = filtered_intensity
    
    return dynamic_threshold, spectrum


def process_mgf_file(mgf_fpath):
    print('Processing file:', mgf_fpath)
    spectrums = get_mgf(mgf_fpath, use_index=config['use_index'], show_progress=True, count=None)
    print('Total spectra in file:', len(spectrums))

    result_dir = config['mgf_result_dpath']
    base_name = os.path.basename(mgf_fpath)
    
    spectrums_filtered = []
    
    for i in tqdm(range(len(spectrums))):
        spec = spectrums[i]
        spec_id = id_of_spectrum(spec)
        print(f'Index {i} | Spectrum ID: {spec_id} | Original Peaks: {len(spec["intensity_array"])}')

        # Apply adaptive denoising
        threshold, denoised_spec = adaptive_plateau_denoising(spec)
        print(f'Denoised | Dynamic Threshold: {threshold:.4f}')

        # Save individual denoised spectrum
        output_path = os.path.join(result_dir, f'{base_name}_denoised_{spec_id}.mgf')
        with open(output_path, 'w') as f_out:
            write_mgf(f_out, denoised_spec)

        spectrums_filtered.append(spec)
    
    print(f'Completed processing: {len(spectrums_filtered)} spectra.')


def main():
    """
    Main entry point for adaptive MS2 denoising.
    """
    print('Starting Adaptive Denoising Workflow...')
    mgf_dir = config['mgf_dpath']
    
    # Locate all MGF files in the target directory
    mgf_files = [os.path.join(mgf_dir, f) for f in os.listdir(mgf_dir) if f.endswith('.mgf')]
    
    for file_path in mgf_files:
        process_mgf_file(file_path)


if __name__ == '__main__':
    # Configuration Setup
    cfg = {
        'mgf_dpath': 'data/denoising/input',
        'mgf_result_dpath': 'data/denoising/output',
        'filter_intensity_min_rate': 0.005,
        'filter_intensity_type': 'mode-',  
        'norm_intensity_type': 'max',  
        'merge_mz_delta': 0.01,
        'filter_2_pepmass_sp1_values': (287.2369, 289.2526,),
        'filter_2_pepmass_sp1_diff': 0.005,
        'filter_2_pepmass_sp1_values_enable': True,
    }
    update_config_on_main(cfg)

    # Ensure output directory exists
    if not os.path.exists(config['mgf_result_dpath']):
        os.makedirs(config['mgf_result_dpath'])

    main()
