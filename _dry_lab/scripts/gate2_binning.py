                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           binning.py
#!/usr/bin/env python3

import os
import glob
import gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def filter_reads(
    summary_file,
    min_length=0,
    max_length=None,
    min_qscore=0,
    max_qscore=None
):
    """
    Filter sequencing reads by read length and quality score.

    Returns
    -------
    filtered_reads : pandas.DataFrame
        Reads satisfying filtering criteria.
    length_stats : pandas.DataFrame
        Read counts, read share, and base share across length bins.
    total_reads : int
        Total reads in the summary file before filtering.
    """

    # Load sequencing summary
    df = pd.read_csv(summary_file, sep="\t")
    total_reads = len(df)

    # Filter by read length
    mask = df["sequence_length_template"] >= min_length

    if max_length is not None:
        mask &= df["sequence_length_template"] <= max_length

    # Filter by quality score
    mask &= df["mean_qscore_template"] >= min_qscore

    if max_qscore is not None:
        mask &= df["mean_qscore_template"] <= max_qscore

    filtered_reads = df[mask].copy()

    # Define length bins
    bins = [0, 200, 500, 651, 1001, 2001, float("inf")]

    labels = [
        "Under 200",
        "200–500",
        "500–650",
        "650–1000",
        "1000–2000",
        "Over 2000"
    ]

    filtered_reads["length_bin"] = pd.cut(
        filtered_reads["sequence_length_template"],
        bins=bins,
        labels=labels,
        right=False
    )

    # Reads and bases per bin
    grouped = (
        filtered_reads
        .groupby("length_bin", observed=False)["sequence_length_template"]
        .agg(["count", "sum"])
    )

    total_passed_reads = len(filtered_reads)
    total_bases = filtered_reads["sequence_length_template"].sum()

    length_stats = (
        grouped
        .rename(columns={"count": "Read count", "sum": "Bases"})
        .rename_axis("Length bin (bp)")
        .reset_index()
    )

    length_stats["Length bin (bp)"] = length_stats["Length bin (bp)"].astype(str)

    # Share of passed reads and passed bases (%)
    if total_passed_reads > 0:
        length_stats["Read share (%)"] = 100 * length_stats["Read count"] / total_passed_reads
        length_stats["Base share (%)"] = 100 * length_stats["Bases"] / total_bases
    else:
        length_stats["Read share (%)"] = 0.0
        length_stats["Base share (%)"] = 0.0

    return filtered_reads, length_stats, total_reads


# ==========================================================
# FILE PATHS
# ==========================================================

file1 = "/users/escott60/batspleen/T23_66_Diyo_RNA_spleen/20261102_1734_MD-101579_FBC09298_cb4c8893/sequencing_summary_FBC09298_cb4c8893_9b6d8c39.txt"
file2 = "/users/escott60/batspleen/T23_66_Diyo_RNA_spleen_2ndRun/BINF6350_Fall2026_RNA01/T23_66_Diyo_RNA_spleen_2ndrun/20261103_1357_MD-101579_FBC09298_013e8950/sequencing_summary_FBC09298_013e8950_973790c0.txt"


# ==========================================================
# FILTERING PARAMETERS
# ==========================================================

min_length = 0
max_length = None

min_qscore = 9
max_qscore = None


# ==========================================================
# PROCESSING RUN 1
# ==========================================================

filtered_reads1, stats1, total_reads1 = filter_reads(
    file1,
    min_length=min_length,
    max_length=max_length,
    min_qscore=min_qscore,
    max_qscore=max_qscore
)


# ==========================================================
# PROCESSING RUN 2
# ==========================================================

filtered_reads2, stats2, total_reads2 = filter_reads(
    file2,
    min_length=min_length,
    max_length=max_length,
    min_qscore=min_qscore,
    max_qscore=max_qscore
)


# ==========================================================
# RUN-LEVEL SUMMARY
# ==========================================================

passed1 = len(filtered_reads1)
passed2 = len(filtered_reads2)

total_bases1 = int(filtered_reads1["sequence_length_template"].sum())
total_bases2 = int(filtered_reads2["sequence_length_template"].sum())

run_summary = pd.DataFrame({
    "Metric": ["Total reads", "Passed reads", "Passed reads (%)", "Passed bases"],
    "Run 1": [f"{total_reads1:,}", f"{passed1:,}", f"{100 * passed1 / total_reads1:.2f}", f"{total_bases1:,}"],
    "Run 2": [f"{total_reads2:,}", f"{passed2:,}", f"{100 * passed2 / total_reads2:.2f}", f"{total_bases2:,}"],
})


# ==========================================================
# COMPARISON TABLE
# ==========================================================

comparison = pd.DataFrame({
    "Length bin (bp)": stats1["Length bin (bp)"],
    "Run 1 reads": stats1["Read count"],
    "Run 2 reads": stats2["Read count"],
    "Run 1 read share (%)": stats1["Read share (%)"],
    "Run 2 read share (%)": stats2["Read share (%)"],
    "Run 1 base share (%)": stats1["Base share (%)"],
    "Run 2 base share (%)": stats2["Base share (%)"],
})

# Add totals row
total_row = pd.DataFrame([{
    "Length bin (bp)": "Total",
    "Run 1 reads": passed1,
    "Run 2 reads": passed2,
    "Run 1 read share (%)": stats1["Read share (%)"].sum(),
    "Run 2 read share (%)": stats2["Read share (%)"].sum(),
    "Run 1 base share (%)": stats1["Base share (%)"].sum(),
    "Run 2 base share (%)": stats2["Base share (%)"].sum(),
}])

comparison = pd.concat([comparison, total_row], ignore_index=True)

share_cols = [c for c in comparison.columns if "share (%)" in c]
comparison[share_cols] = comparison[share_cols].round(2)


# ==========================================================
# OUTPUT
# ==========================================================

print("\n==========================================================")
print("              SEQUENCING RUN COMPARISON")
print("==========================================================")

print(f"\nFiltering criteria:")
print(f"Minimum length: {min_length} bp")
print(f"Maximum length: {max_length if max_length is not None else 'None'}")
print(f"Minimum Q-score: {min_qscore}")
print(f"Maximum Q-score: {max_qscore if max_qscore is not None else 'None'}")

print("\nRun-level summary:")
print(run_summary.to_string(index=False))

print("\nRead counts and shares of passed reads/bases per length bin:")
print(comparison.to_string(index=False))

print("\n==========================================================")

# Optional: Export results
comparison.to_csv("sequencing_run_comparison.csv", index=False)
run_summary.to_csv("sequencing_run_summary.csv", index=False)


# ==========================================================
# PLOT: READ LENGTH DISTRIBUTIONS WITH BINS MARKED
# ==========================================================

# Keep in sync with the bins in filter_reads()
BIN_EDGES = [0, 200, 500, 651, 1001, 2001]
BIN_LABELS = ["Under 200", "200–500", "500–650", "650–1000", "1000–2000", "Over 2000"]

def plot_length_distributions(runs, out_png):
    """runs: list of (run_name, filtered_reads_df)."""
    all_lengths = np.concatenate(
        [df["sequence_length_template"].clip(lower=1).to_numpy() for _, df in runs]
    )
    xmin, xmax = all_lengths.min(), all_lengths.max()
    hist_bins = np.logspace(np.log10(xmin), np.log10(xmax), 150)

    fig, axes = plt.subplots(len(runs), 1, figsize=(11, 4.2 * len(runs)), sharex=True)
    axes = np.atleast_1d(axes)

    # Bin spans: [lower, upper) for each label, clipped to the data range
    spans = list(zip(BIN_EDGES, BIN_EDGES[1:] + [xmax * 1.05]))

    for ax, (name, df) in zip(axes, runs):
        lengths = df["sequence_length_template"].clip(lower=1)
        ax.hist(lengths, bins=hist_bins, color="#4C72B0", edgecolor="none")
        ax.set_xscale("log")
        ax.set_xlim(xmin, xmax * 1.05)

        counts = df["length_bin"].value_counts(sort=False)

        # Add headroom above the tallest bar so labels sit clear of the data
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
        ymax = ax.get_ylim()[1]
        levels = [0.97, 0.88, 0.79]   # stagger label heights so neighbours don't collide

        for i, ((lo, hi), label) in enumerate(zip(spans, BIN_LABELS)):
            lo_c, hi_c = max(lo, xmin), min(hi, xmax * 1.05)
            if lo_c >= hi_c:
                continue
            if i % 2 == 0:
                ax.axvspan(lo_c, hi_c, color="grey", alpha=0.08, lw=0)
            if lo > xmin:
                ax.axvline(lo, color="black", ls="--", lw=0.8)
            x_text = np.sqrt(lo_c * hi_c)  # geometric centre suits a log axis
            ax.text(
                x_text, ymax * levels[i % 3],
                f"{label}\n{int(counts.get(label, 0)):,} reads",
                ha="center", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
            )

        ax.set_ylabel("Passed reads")
        ax.set_title(f"{name}: read length distribution ({len(df):,} passed reads)")

    axes[-1].set_xlabel("Read length (bp, log scale)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved plot: {out_png}")

plot_length_distributions(
    [("Run 1", filtered_reads1), ("Run 2", filtered_reads2)],
    "read_length_distributions.png",
)


# ==========================================================
# EXTRACT FASTQ RECORDS INTO ONE FILE PER BIN
# ==========================================================

# Where to look for FASTQ files (searched recursively, .fastq and .fastq.gz).
# Default: the folder holding each summary file. If your FASTQs live elsewhere
# (e.g. a Dorado output folder), point these at it.
FASTQ_DIR1 = os.path.dirname(file1)
FASTQ_DIR2 = os.path.dirname(file2)
FASTQ_OUT_DIR = "fastq_by_bin"

def safe_name(label):
    return label.replace("–", "-").replace(" ", "_")

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def extract_fastq_by_bin(filtered_reads, search_dir, out_dir, run_tag):
    os.makedirs(out_dir, exist_ok=True)

    # read_id -> bin label, for passed reads only
    bin_of = dict(
        zip(filtered_reads["read_id"], filtered_reads["length_bin"].astype(str))
    )
    labels = list(filtered_reads["length_bin"].cat.categories)

    handles = {
        lab: open(os.path.join(out_dir, f"{run_tag}_{safe_name(lab)}.fastq"), "w")
        for lab in labels
    }
    written = dict.fromkeys(labels, 0)
    seen = set()

    fastq_files = sorted(
        p for p in glob.glob(os.path.join(search_dir, "**", "*"), recursive=True)
        if p.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
    )
    print(f"{run_tag}: scanning {len(fastq_files)} FASTQ file(s) under {search_dir}")

    try:
        for path in fastq_files:
            with open_text(path) as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline()
                    plus = f.readline()
                    qual = f.readline()

                    read_id = header[1:].split()[0].strip()
                    lab = bin_of.get(read_id)
                    if lab is None or read_id in seen:
                        continue
                    seen.add(read_id)

                    # rstrip handles any CRLF line endings from Windows-written files
                    handles[lab].write(
                        header.rstrip("\r\n") + "\n"
                        + seq.rstrip("\r\n") + "\n"
                        + plus.rstrip("\r\n") + "\n"
                        + qual.rstrip("\r\n") + "\n"
                    )
                    written[lab] += 1
    finally:
        for h in handles.values():
            h.close()

    # Check against the summary-based counts
    expected = filtered_reads["length_bin"].value_counts(sort=False)
    print(f"\n{run_tag}: FASTQ records written vs passed reads in summary")
    for lab in labels:
        exp = int(expected[lab])
        flag = "" if written[lab] == exp else "   <-- MISMATCH"
        print(f"  {lab:<10} {written[lab]:>10,} / {exp:>10,}{flag}")
    print(f"  {'Total':<10} {sum(written.values()):>10,} / {len(filtered_reads):>10,}\n")

extract_fastq_by_bin(filtered_reads1, FASTQ_DIR1, FASTQ_OUT_DIR, "run1")
extract_fastq_by_bin(filtered_reads2, FASTQ_DIR2, FASTQ_OUT_DIR, "run2")










