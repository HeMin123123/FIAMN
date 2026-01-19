import numpy as np
from matchms.importing import load_from_mgf
from matchms.similarity import CosineGreedy
from matchms.filtering import normalize_intensities
from matchms import Spectrum


# --- Internal Helper Functions ---
def merge_peaks_internal(spectrum: Spectrum, tolerance=0.02) -> Spectrum:
    """
    Merges peaks within a certain m/z tolerance, keeping the highest intensity peak.
    """
    if spectrum is None: return None
    mz, intensities = spectrum.peaks.mz, spectrum.peaks.intensities
    # Sort indices by intensity in descending order
    idx_sorted = np.argsort(intensities)[::-1]
    keep_idx, seen_mz = [], []
    for i in idx_sorted:
        # Check if the current peak is within tolerance of already kept peaks
        if not any(abs(mz[i] - ref) <= tolerance for ref in seen_mz):
            keep_idx.append(i)
            seen_mz.append(mz[i])
    keep_idx.sort()
    return Spectrum(mz=mz[keep_idx], intensities=intensities[keep_idx], metadata=spectrum.metadata.copy())


def save_as_pepmass_rt_mgf(spectrums, output_file):
    """
    Saves the list of spectrums back to an MGF file including PEPMASS and RT metadata.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for spec in spectrums:
            f.write("BEGIN IONS\n")
            meta = spec.metadata
            f.write(f"TITLE={meta.get('title', '')}\n")
            # Try both possible metadata keys for Retention Time
            rt = meta.get('retention_time') or meta.get('rtinseconds')
            if rt is not None: f.write(f"RTINSECONDS={rt}\n")
            # Try both possible metadata keys for Precursor M/Z
            pepmass = meta.get('precursor_mz') or meta.get('pepmass')
            if pepmass is not None: f.write(f"PEPMASS={pepmass}\n")
            for m, i in zip(spec.peaks.mz, spec.peaks.intensities):
                f.write(f"{m} {i}\n")
            f.write("END IONS\n\n")


# --- Improved Search Logic Function ---
def run_final_restoration_v2(input_path, output_path, target_mz, target_rt, threshold=0.7):
    """
    Finds a seed spectrum by MZ/RT and performs a network-based similarity search
    to extract related spectra.
    """
    print(f"--- Loading MGF file ---")
    # Enable metadata_harmonization to ensure standard access to mz/rt
    raw_spectrums = list(load_from_mgf(input_path, metadata_harmonization=True))

    # 1. Search for Seed Index using MZ and RT
    seed_idx = None
    print(f"🔍 Searching for target: m/z={target_mz}, RT={target_rt}...")

    for i, s in enumerate(raw_spectrums):
        curr_mz = s.get("precursor_mz")
        curr_rt = s.get("retention_time")
        if curr_mz and curr_rt:
            # Set a small error tolerance for matching
            if abs(curr_mz - target_mz) < 0.005 and abs(curr_rt - target_rt) < 2.0:
                seed_idx = i
                print(f"✅ Seed spectrum found! Title: {s.get('title')[:50]}...")
                break

    if seed_idx is None:
        print("❌ Could not find spectrum via m/z and RT. Please check your parameters.")
        return

    # 2. Pre-processing and Clustering
    processed_spectrums = [normalize_intensities(merge_peaks_internal(s)) for s in raw_spectrums]
    cosine_sim = CosineGreedy(tolerance=0.02)
    selected_indices, to_explore = {seed_idx}, [seed_idx]

    print(f"🚀 Searching for similar spectra (threshold: {threshold})...")
    # Breadth-first search for connected similarity network
    while to_explore:
        curr = to_explore.pop(0)
        for i, other in enumerate(processed_spectrums):
            if i in selected_indices: continue
            score = cosine_sim.pair(processed_spectrums[curr], other)['score']
            if score >= threshold:
                selected_indices.add(i)
                to_explore.append(i)

    # 3. Save Results
    final_specs = [processed_spectrums[i] for i in selected_indices]
    save_as_pepmass_rt_mgf(final_specs, output_path)
    print(f"✨ Done! Found {len(final_specs)} related spectra.")


if __name__ == "__main__":
    # --- User Configuration ---
    # Based on your provided info:
    # RT for this TITLE should be 1201.59
    # PEPMASS (m/z) should be obtained from the original file record

    TARGET_MZ = 289.253430034396  # <--- Replace with the PEPMASS value from your file
    TARGET_RT = 1201.59           # Retention Time

    run_final_restoration_v2("filtered_TIC.mgf", "family_network.mgf", TARGET_MZ, TARGET_RT)