"""
Fetch NHANES data for IPD-QMA demonstration.

Uses pyreadstat to download SAS XPT transport files from CDC.
Fetches demographics + blood pressure + body measures + medications.

This gives us real IPD with:
  - Continuous outcomes: systolic BP, BMI
  - Binary treatment: antihypertensive medication use
  - Multiple "studies" = NHANES survey cycles
"""

import os
import tempfile
import urllib.request
import numpy as np
import pandas as pd

try:
    import pyreadstat
except ImportError:
    raise ImportError("Install pyreadstat: pip install pyreadstat")


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# NHANES cycles and their file suffixes
NHANES_CYCLES = [
    {"cycle": "2017-2018", "yr": "2017_2018", "suf": "J"},
    {"cycle": "2015-2016", "yr": "2015_2016", "suf": "I"},
    {"cycle": "2013-2014", "yr": "2013_2014", "suf": "H"},
    {"cycle": "2011-2012", "yr": "2011_2012", "suf": "G"},
]

# Files to fetch per cycle
NHANES_FILES = {
    "demo": "DEMO",        # Demographics
    "bpx": "BPX",          # Blood pressure
    "bmx": "BMX",          # Body measures
    "bpq": "BPQ",          # Blood pressure questionnaire (medication)
    "rxq_rx": "RXQ_RX",    # Prescription medications
}


def download_xpt(url: str) -> pd.DataFrame:
    """Download a SAS XPT file from URL and return as DataFrame."""
    fd, tf_path = tempfile.mkstemp(suffix=".xpt")
    os.close(fd)
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
        })
        resp = urllib.request.urlopen(req, timeout=60)
        with open(tf_path, 'wb') as f:
            f.write(resp.read())
        df = pd.read_sas(tf_path, format='xport')
        df.columns = [c.lower() for c in df.columns]
        return df
    finally:
        try:
            os.unlink(tf_path)
        except OSError:
            pass


def fetch_cycle(cycle_info: dict) -> pd.DataFrame:
    """Fetch all files for one NHANES cycle and merge."""
    suf = cycle_info["suf"]
    cycle_name = cycle_info["cycle"]
    yr = cycle_info["yr"]

    # Correct CDC URL: /Nchs/Data/Nhanes/Public/{start_year}/DataFiles/
    start_year = yr.split("_")[0]
    base_url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{start_year}/DataFiles"

    dfs = {}
    for key, prefix in NHANES_FILES.items():
        url = f"{base_url}/{prefix}_{suf}.XPT"

        print(f"  Fetching {prefix}_{suf}.XPT ...", end=" ")
        try:
            df = download_xpt(url)
            dfs[key] = df
            print(f"OK ({len(df)} rows)")
        except Exception as e:
            print(f"FAILED ({e})")

    if "demo" not in dfs:
        print(f"  ERROR: No demographics for {cycle_name}")
        return pd.DataFrame()

    # Start with demographics
    merged = dfs["demo"]

    # Merge blood pressure
    if "bpx" in dfs:
        bp_cols = ["seqn"]
        # Systolic BP readings
        for col in dfs["bpx"].columns:
            if col.startswith("bpxs") or col.startswith("bpxd") or col == "seqn":
                if col not in bp_cols:
                    bp_cols.append(col)
        merged = merged.merge(dfs["bpx"][bp_cols], on="seqn", how="left")

    # Merge body measures
    if "bmx" in dfs:
        bmx_cols = ["seqn", "bmxbmi", "bmxht", "bmxwt"]
        bmx_cols = [c for c in bmx_cols if c in dfs["bmx"].columns]
        merged = merged.merge(dfs["bmx"][bmx_cols], on="seqn", how="left")

    # Merge BP questionnaire (medication use)
    if "bpq" in dfs:
        bpq_cols = ["seqn"]
        for col in dfs["bpq"].columns:
            if col.startswith("bpq") or col == "seqn":
                if col not in bpq_cols:
                    bpq_cols.append(col)
        merged = merged.merge(dfs["bpq"][bpq_cols], on="seqn", how="left")

    merged["cycle"] = cycle_name
    merged["study_id"] = cycle_name
    merged["patient_id"] = merged["seqn"].astype(int)

    return merged


def compute_mean_sbp(row):
    """Compute mean systolic BP from readings 2-4 (skip first per NHANES protocol)."""
    readings = []
    for i in range(2, 5):  # Skip bpxsy1 (white-coat effect)
        col = f"bpxsy{i}"
        if col in row.index and pd.notna(row[col]) and row[col] > 0:
            readings.append(row[col])
    if len(readings) >= 1:
        return np.mean(readings)
    return np.nan


def compute_mean_dbp(row):
    """Compute mean diastolic BP from readings 2-4 (skip first per NHANES protocol)."""
    readings = []
    for i in range(2, 5):  # Skip bpxdi1 (white-coat effect)
        col = f"bpxdi{i}"
        if col in row.index and pd.notna(row[col]) and row[col] >= 0:
            readings.append(row[col])
    if len(readings) >= 1:
        return np.mean(readings)
    return np.nan


def fetch_all_nhanes(n_cycles: int = 4) -> pd.DataFrame:
    """
    Fetch and combine NHANES data across cycles.

    Returns
    -------
    pd.DataFrame with columns:
        patient_id, study_id, cycle, age, sex, mean_sbp, mean_dbp,
        bmi, on_bp_meds (binary treatment indicator)
    """
    os.makedirs(DATA_DIR, exist_ok=True)
    cache_path = os.path.join(DATA_DIR, "nhanes_combined.csv")

    if os.path.exists(cache_path):
        print(f"Loading cached data from {cache_path}")
        return pd.read_csv(cache_path)

    print(f"Fetching NHANES data ({n_cycles} cycles)...\n")
    all_dfs = []

    for i, cycle_info in enumerate(NHANES_CYCLES[:n_cycles]):
        print(f"[{i+1}/{n_cycles}] Cycle: {cycle_info['cycle']}")
        df = fetch_cycle(cycle_info)
        if len(df) > 0:
            all_dfs.append(df)
        print()

    if not all_dfs:
        raise RuntimeError("No NHANES data fetched successfully")

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"Combined: {len(combined)} rows, {len(combined.columns)} columns")

    # --- Derive clean analysis variables ---
    # Age
    if "ridageyr" in combined.columns:
        combined["age"] = combined["ridageyr"]
    elif "ridagemn" in combined.columns:
        combined["age"] = combined["ridagemn"] / 12.0

    # Sex (1=Male, 2=Female in NHANES)
    if "riagendr" in combined.columns:
        combined["sex"] = combined["riagendr"].map({1: "Male", 2: "Female"})

    # Mean systolic/diastolic BP
    combined["mean_sbp"] = combined.apply(compute_mean_sbp, axis=1)
    combined["mean_dbp"] = combined.apply(compute_mean_dbp, axis=1)

    # BMI
    if "bmxbmi" in combined.columns:
        combined["bmi"] = combined["bmxbmi"]

    # On BP medications (BPQ050A: "Now taking prescribed medicine for HBP?")
    # 1 = Yes, 2 = No; only asked when BPQ020 (told have HBP) = 1
    # Preserve NaN for people never asked (no HBP history) to avoid
    # conflating normotensives with unmedicated hypertensives.
    if "bpq050a" in combined.columns:
        combined["on_bp_meds"] = np.where(
            combined["bpq050a"] == 1, 1,
            np.where(combined["bpq050a"] == 2, 0, np.nan)
        )
    else:
        combined["on_bp_meds"] = np.nan

    # Keep analysis columns
    analysis_cols = [
        "patient_id", "study_id", "cycle", "age", "sex",
        "mean_sbp", "mean_dbp", "bmi", "on_bp_meds"
    ]
    analysis_cols = [c for c in analysis_cols if c in combined.columns]
    result = combined[analysis_cols].copy()

    # Save cache
    result.to_csv(cache_path, index=False)
    print(f"Saved to {cache_path}")

    return result


def prepare_ipd_qma_data(
    df: pd.DataFrame,
    outcome_col: str = "mean_sbp",
    treatment_col: str = "on_bp_meds",
    study_col: str = "study_id",
    min_per_arm: int = 30,
):
    """
    Prepare NHANES data for IPD-QMA analysis.

    Returns
    -------
    studies_data : list of (control_array, treatment_array) tuples
    labels : list of str
    summary : dict
    """
    # NOTE: Survey weights, PSUs, and stratification are ignored for simplicity.
    # Results are not population-representative. Design effects would increase SEs.
    # This is an observational comparison (confounding by indication expected).
    # Filter to complete cases
    subset = df.dropna(subset=[outcome_col, treatment_col, study_col])
    subset = subset[subset[outcome_col] > 0]

    studies_data = []
    labels = []
    summary = {"studies": []}

    for study_id, group in subset.groupby(study_col):
        control = group[group[treatment_col] == 0][outcome_col].values
        treatment = group[group[treatment_col] == 1][outcome_col].values

        if len(control) >= min_per_arm and len(treatment) >= min_per_arm:
            studies_data.append((control, treatment))
            labels.append(str(study_id))
            summary["studies"].append({
                "study_id": study_id,
                "n_control": len(control),
                "n_treatment": len(treatment),
                "mean_control": np.mean(control),
                "mean_treatment": np.mean(treatment),
                "sd_control": np.std(control, ddof=1),
                "sd_treatment": np.std(treatment, ddof=1),
            })

    summary["n_studies"] = len(studies_data)
    summary["total_patients"] = sum(s["n_control"] + s["n_treatment"] for s in summary["studies"])
    summary["outcome"] = outcome_col
    summary["treatment"] = treatment_col

    return studies_data, labels, summary


if __name__ == "__main__":
    df = fetch_all_nhanes(n_cycles=4)
    print(f"\nDataset shape: {df.shape}")
    print(f"\nColumn types:\n{df.dtypes}")
    print(f"\nMissing values:\n{df.isnull().sum()}")
    print(f"\nSample:\n{df.head()}")
