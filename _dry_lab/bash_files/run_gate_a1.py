#!/usr/bin/env python3
"""
Stage A1 gate numbers.
Reports passed reads, passed bases (Gb, 3 s.f.), and % unblock_mux_change
for each run. Also demonstrates and fixes the CRLF bug.
"""
import argparse
import glob
import os
import sys

import pandas as pd


def find_summary(run_dir):
    hits = sorted(glob.glob(os.path.join(run_dir, "**", "sequencing_summary_*.txt"), recursive=True))
    if len(hits) != 1:
        sys.exit(f"Expected 1 sequencing_summary in {run_dir}, found {len(hits)}: {hits}")
    return hits[0]


def crlf_check(path):
    with open(path, "rb") as fh:
        return fh.readline().endswith(b"\r\n")


def load_buggy(path):
    return pd.read_csv(path, sep="\t", lineterminator="\n", low_memory=False)


def load_clean(path):
    df = pd.read_csv(path, sep="\t", lineterminator="\n", low_memory=False)
    df.columns = df.columns.str.strip()
    for col in df.columns:
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.strip()
    if pd.api.types.is_string_dtype(df["passes_filtering"]):
        df["passes_filtering"] = df["passes_filtering"].str.upper().map({"TRUE": True, "FALSE": False})
    return df


def gate_numbers(run_dir):
    path = find_summary(run_dir)
    run = os.path.basename(os.path.normpath(run_dir))

    buggy = load_buggy(path)
    last_col = buggy.columns[-1]
    buggy_hits = (buggy[last_col] == "unblock_mux_change").sum()

    df = load_clean(path)
    n_total = len(df)
    passed = df[df["passes_filtering"]]
    n_pass = len(passed)
    pass_bases = passed["sequence_length_template"].sum()
    n_unblock = (df["end_reason"] == "unblock_mux_change").sum()
    n_unblock_pass = (passed["end_reason"] == "unblock_mux_change").sum()

    res = {
        "run": run,
        "summary_file": path,
        "n_columns": df.shape[1],
        "header_has_CRLF": crlf_check(path),
        "last_column_raw_name": repr(last_col),
        "unblock_hits_BUGGY_read": int(buggy_hits),
        "total_reads": n_total,
        "duplicate_read_ids": int(df["read_id"].duplicated().sum()),
        "passed_reads": n_pass,
        "passed_bases": int(pass_bases),
        "passed_bases_Gb_3sf": float(f"{pass_bases / 1e9:.3g}"),
        "unblock_mux_change_reads": int(n_unblock),
        "unblock_pct_of_all_reads": round(100 * n_unblock / n_total, 2),
        "unblock_pct_of_passed_reads": round(100 * n_unblock_pass / n_pass, 2),
    }
    return res, df["end_reason"].value_counts()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--out", default="gate_a1_results.tsv")
    args = ap.parse_args()

    print(f"pandas {pd.__version__}  python {sys.version.split()[0]}\n")
    rows = []
    for rd in args.run_dirs:
        res, end_reasons = gate_numbers(rd)
        rows.append(res)
        print(f"===== {res['run']} =====")
        for k, v in res.items():
            print(f"  {k:30s} {v}")
        print("  end_reason counts:")
        for reason, n in end_reasons.items():
            print(f"    {reason:30s} {n}")
        print()

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()