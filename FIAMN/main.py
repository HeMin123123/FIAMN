import os
import time

# Importing core functions from individual modules
from spectral_denoising.script_01_spectral_denoising import filter_mgf_by_nonzero_mode
from TIC_filtering.script_02_TIC_filtering import filter_mgf_by_total_intensity
from family_Network_Cosine.script_03_cosine_similarity_calc import run_final_restoration_v2
from precursor_formula_filter.script_04_precursor_formula_filter import filter_mgf_with_nitrogen_rule
from skeleton_fragment_filter.script_05_skeleton_fragment_filter import main as run_skeleton_fragment_filter
from skeleton_similarity_network.script_06_skeleton_sim import run_similarity_network_pipeline


def run_fiamn_pipeline():
    print("=" * 50)
    print("Starting FIAMN Molecular Networking Analysis Pipeline")
    print("=" * 50)
    start_time = time.time()

    # --- Configuration Parameters ---
    raw_input = "spectral_denoising/Text.mgf"  # Raw input file
    target_mz = 289.253430034396  # Target seed ion m/z (Step 3)
    target_rt = 1201.59  # Target seed ion RT (Step 3)

    # Intermediate file paths
    denoised_mgf = "01_denoising.mgf"
    tic_filtered_mgf = "02_filtered_TIC.mgf"
    family_network_mgf = "03_family_network.mgf"
    precursor_filtered_mgf = "04_precursor_formula_filter.mgf"
    skeleton_filtered_mgf = "05_skeleton_fragment_filter.mgf"
    final_xgmml = "06_skeleton_similarity_network.xgmml"

    # 1). Noise filtering (Mode-based denoising)
    print("\n[Step 1/6] Performing mode-based noise reduction...")
    # Eliminates non-specific baseline noise while retaining valid low-abundance peaks
    filter_mgf_by_nonzero_mode(raw_input, denoised_mgf)

    # 2). Low-abundance signal filtering (TIC filtering)
    print("\n[Step 2/6] Filtering signals by Total Ion Current (Threshold >= 2000)...")
    # Removes data with extremely low content to retain valid precursor signals
    filter_mgf_by_total_intensity(denoised_mgf, tic_filtered_mgf, min_total_intensity=2000)

    # 3). Initial target molecular family screening
    print(f"\n[Step 3/6] Constructing initial family using seed m/z {target_mz}...")
    # Uses square root transformation of intensities for cosine similarity assessment
    run_final_restoration_v2(tic_filtered_mgf, family_network_mgf, target_mz, target_rt, threshold=0.7)

    # 4). Precursor formula re-screening
    print("\n[Step 4/6] Filtering for target diterpenoids (C18-C25)...")
    # Retains ions with carbon counts between 18 and 25 based on precursor formulas
    filter_mgf_with_nitrogen_rule(family_network_mgf, precursor_filtered_mgf)

    # 5). Skeleton-associated fragment retention
    print("\n[Step 5/6] Extracting skeletal feature fragments (C and H only)...")
    # Selectively retains fragment peaks composed exclusively of C and H atoms
    run_skeleton_fragment_filter(precursor_filtered_mgf, skeleton_filtered_mgf)

    # 6). Skeleton similarity assessment and network generation
    print("\n[Step 6/6] Assessing skeleton similarity and generating XGMML file...")
    # Evaluates topological similarity via cosine similarity without square root transformation
    run_similarity_network_pipeline(skeleton_filtered_mgf, final_xgmml)

    # --- Pipeline Completion ---
    end_time = time.time()
    elapsed = end_time - start_time
    print("\n" + "=" * 50)
    print(f"Pipeline executed successfully! Total time: {elapsed:.2f} seconds")
    print(f"Final network file generated: {final_xgmml}")
    print("=" * 50)


if __name__ == "__main__":
    if os.path.exists("spectral_denoising/Text.mgf"):
        run_fiamn_pipeline()
    else:
        print("Error: Raw input file 'Text.mgf' not found.")
