import re
import math
import numpy as np
from typing import List, Dict, Tuple
from itertools import combinations

# ===================== Core Configuration =====================
SIM_THRESHOLD = 0.9  # Cosine similarity threshold
MZ_TOLERANCE_DA = 0.01  # Mass tolerance (Da)
ENTROPY_RATIO_LIMIT = 0.75  # Spectral entropy ratio threshold; values below this set line to dashed


def calculate_spectral_entropy(peaks_dict: Dict[float, float]) -> float:
    """Calculates the Spectral Entropy of a single spectrum."""
    intensities = np.array(list(peaks_dict.values()))
    # Remove anomalies with intensity <= 0
    p = intensities[intensities > 0]
    if len(p) <= 1:
        return 0.0
    # Normalize intensity distribution
    prob = p / np.sum(p)
    # Shannon entropy formula
    return -np.sum(prob * np.log(prob + 1e-12))


def parse_mgf(file_path: str) -> List[Dict]:
    """Parses an MGF file and extracts spectral data."""
    mgf_entries = []
    current_entry, peaks = None, []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('BEGIN IONS'):
                current_entry, peaks = {}, []
            elif line.startswith('END IONS'):
                if current_entry is not None and peaks:
                    current_entry['peaks'] = peaks
                    mgf_entries.append(current_entry)
            elif '=' in line and current_entry is not None:
                key, value = line.split('=', 1)
                current_entry[key.strip()] = value.strip()
            elif current_entry is not None and re.match(r'^\d+\.?\d*\s+\d+\.?\d*$', line):
                try:
                    mz, intensity = map(float, line.split())
                    peaks.append((mz, intensity))
                except ValueError:
                    continue
    return mgf_entries


def generate_spectrum_id(entry: Dict) -> str:
    """Generates a unique ID for a spectrum based on Pepmass and RT."""
    pepmass = entry.get('PEPMASS', '0.0').split()[0]
    rt = entry.get('RTINSECONDS', '0')
    return f"{round(float(pepmass), 2):.2f}-{int(round(float(rt)))}"


def cosine_similarity(peaks1: Dict[float, float], peaks2: Dict[float, float]) -> float:
    """Calculates the basic cosine similarity between two peak lists."""
    keys1, keys2 = list(peaks1.keys()), list(peaks2.keys())
    dot_product, matched = 0.0, set()
    for mz1, int1 in peaks1.items():
        for idx2, mz2 in enumerate(keys2):
            if idx2 not in matched and abs(mz1 - mz2) <= MZ_TOLERANCE_DA:
                dot_product += int1 * peaks2[mz2]
                matched.add(idx2)
                break
    n1 = math.sqrt(sum(i ** 2 for i in peaks1.values()))
    n2 = math.sqrt(sum(i ** 2 for i in peaks2.values()))
    return dot_product / (n1 * n2) if n1 and n2 else 0.0


def save_cytoscape_xgmml(file_path: str, all_nodes: List[Dict], edges: List[Dict]):
    """Generates an XGMML file containing node/edge logic for Cytoscape visualization."""
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<graph xmlns="http://www.cs.rpi.edu/XGMML" directed="0" label="Diterpene Entropy Network">\n')

        # Write nodes
        for node in all_nodes:
            label = f"{node['pepmass']:.2f}&#10;{node['retention_time']}"
            f.write(f'  <node id="{node["id"]}" label="{label}">\n')
            f.write(f'    <att name="pepmass" type="real" value="{node["pepmass"]}"/>\n')
            f.write(f'    <att name="rt" type="integer" value="{node["retention_time"]}"/>\n')
            f.write(f'    <att name="spectral_entropy" type="real" value="{node["entropy"]}"/>\n')
            f.write('  </node>\n')

        # Write edges (including dashed line logic)
        for i, edge in enumerate(edges):
            f.write(f'  <edge id="{i}" source="{edge["source"]}" target="{edge["target"]}">\n')
            f.write(f'    <att name="cosine_similarity" type="real" value="{edge["similarity"]}"/>\n')
            f.write(f'    <att name="entropy_ratio" type="real" value="{edge["entropy_ratio"]}"/>\n')
            # Set line style based on entropy ratio (Cytoscape attribute)
            line_style = "dash" if edge["entropy_ratio"] < ENTROPY_RATIO_LIMIT else "solid"
            f.write(f'    <att name="line_style" type="string" value="{line_style}"/>\n')
            f.write('  </edge>\n')

        f.write('</graph>\n')
    print(f"✅ XGMML visualization file saved: {file_path} (Edges: {len(edges)})")


def run_similarity_network_pipeline(input_mgf: str, output_xgmml: str):
    """
    Main pipeline: Parses spectra, calculates entropy, computes similarity network,
    and exports to XGMML.
    """
    print("🔍 Parsing spectra and calculating spectral entropy...")
    entries = parse_mgf(input_mgf)
    processed_data = []

    for entry in entries:
        sid = generate_spectrum_id(entry)
        peaks = {mz: i for mz, i in entry.get('peaks', [])}
        entropy = calculate_spectral_entropy(peaks)

        # Extract metadata for node display
        pepmass = float(sid.split('-')[0])
        rt = int(sid.split('-')[1])

        processed_data.append({
            'id': sid, 'peaks': peaks, 'entropy': entropy,
            'pepmass': pepmass, 'retention_time': rt
        })

    edges = []
    print(f"📈 Calculating similarity and entropy ratios (Threshold: {SIM_THRESHOLD})...")
    for (i, d1), (j, d2) in combinations(enumerate(processed_data), 2):
        cos_sim = cosine_similarity(d1['peaks'], d2['peaks'])

        if cos_sim >= SIM_THRESHOLD:
            # Calculate entropy ratio: min(entropy) / max(entropy)
            e1, e2 = d1['entropy'], d2['entropy']
            if e1 > 0 and e2 > 0:
                ratio = min(e1, e2) / max(e1, e2)
            else:
                ratio = 1.0  # Avoid division by zero, default to solid line

            edges.append({
                'source': d1['id'], 'target': d2['id'],
                'similarity': round(cos_sim, 4),
                'entropy_ratio': round(ratio, 4)
            })

    save_cytoscape_xgmml(output_xgmml, processed_data, edges)


if __name__ == "__main__":
    # File configuration
    INPUT_FILE = "skeleton_fragment_filter.mgf"
    OUTPUT_FILE = "skeleton_similarity_network.xgmml"

    run_similarity_network_pipeline(INPUT_FILE, OUTPUT_FILE)