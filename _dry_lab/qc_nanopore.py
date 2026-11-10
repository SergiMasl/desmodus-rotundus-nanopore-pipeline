import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # write files, do not open windows
import matplotlib.pyplot as plt

PASS_Q = 9.0                    # MinKNOW's pass/fail threshold for RNA004

# ---------------------------------------------------------------- 1. LOAD
# NOTE: MinKNOW runs on Windows, so these files use CRLF line endings.
# If you read them without handling that, the LAST column silently carries a
# stray '\r' and string comparisons on it fail. This is a real and common bug.
# pandas handles it; awk and cut do not. Always check your last column.

def load(path):
    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
        usecols=["read_id", "parent_read_id", "channel", "mux", "start_time",
                 "duration", "passes_filtering", "sequence_length_template",
                 "mean_qscore_template", "end_reason"],
    )
    # strip whitespace/CR from any object column
    for c in df.select_dtypes("object"):
        df[c] = df[c].str.strip()
    df["passes_filtering"] = df["passes_filtering"].astype(str).str.upper().eq("TRUE")
    df["hour"] = df["start_time"] / 3600.0
    return df


# ------------------------------------------------------- 2. SUMMARY NUMBERS
def summarise(df, label, out):
    p = df[df.passes_filtering]
    L = p.sequence_length_template.values
    srt = np.sort(L)[::-1]
    n50 = srt[np.searchsorted(np.cumsum(srt), srt.sum() / 2)]

    rows = [
        ("total reads",            f"{len(df):,}"),
        ("passed reads",           f"{len(p):,} ({100*len(p)/len(df):.1f}%)"),
        ("failed reads",           f"{len(df)-len(p):,}"),
        ("passed bases",           f"{L.sum()/1e9:.3f} Gb"),
        ("mean read length",       f"{L.mean():.0f} bp"),
        ("median read length",     f"{np.median(L):.0f} bp"),
        ("read length N50",        f"{n50:,} bp"),
        ("longest read",           f"{L.max():,} bp"),
        ("median qscore (passed)", f"{p.mean_qscore_template.median():.1f}"),
        ("run duration",           f"{df.hour.max():.1f} h"),
        ("channels producing reads", f"{df.channel.nunique():,}"),
        ("split reads",            f"{(df.read_id != df.parent_read_id).sum():,}"),
    ]
    print(f"\n{'='*58}\n  {label}\n{'='*58}")
    for k, v in rows:
        print(f"  {k:<26} {v}")
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(
        f"{out}_summary.csv", index=False)
    return rows


# ------------------------------------------------------------- 3. FIGURES
def figures(df, label, out):
    p = df[df.passes_filtering]
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(f"Nanopore direct RNA QC  |  {label}", fontsize=13)

    # (a) read length, log x. Transcript lengths span orders of magnitude,
    #     so a linear axis hides everything above ~2 kb.
    a = ax[0, 0]
    bins = np.logspace(1.5, 4.7, 70)
    a.hist(p.sequence_length_template.clip(lower=1), bins=bins, color="#4a7ba7")
    a.set_xscale("log")
    a.set_xlabel("read length (bp, log)"); a.set_ylabel("passed reads")
    a.set_title("(a) Read length distribution")

    # (b) quality. The dip at Q9 is the pass/fail filter, not biology.
    a = ax[0, 1]
    a.hist(df.mean_qscore_template, bins=80, color="#8a9a5b")
    a.axvline(PASS_Q, color="crimson", ls="--", lw=1.2)
    a.annotate("pass threshold Q9", (PASS_Q, a.get_ylim()[1]*0.9),
               xytext=(6, 0), textcoords="offset points",
               color="crimson", fontsize=8)
    a.set_xlabel("mean qscore"); a.set_ylabel("all reads")
    a.set_title("(b) Quality distribution")

    # (c) length vs quality. Short+low-quality reads are usually adapter-only
    #     or aborted translocations.
    a = ax[0, 2]
    # zero-length reads exist (adapter-only / aborted) and break a log axis
    nz = df[df.sequence_length_template > 0]
    s = nz.sample(min(60000, len(nz)), random_state=0)
    a.hexbin(s.sequence_length_template, s.mean_qscore_template,
             xscale="log", gridsize=45, mincnt=1, cmap="viridis")
    a.axhline(PASS_Q, color="crimson", ls="--", lw=1.0)
    a.set_xlabel("read length (bp, log)"); a.set_ylabel("mean qscore")
    a.set_title("(c) Length vs quality")

    # (d) THE KEY PLOT: yield over time. Direct RNA declines because pores
    #     block irreversibly. A smooth decay is normal. A cliff is a failure.
    a = ax[1, 0]
    hr = np.floor(df.hour).astype(int)
    a.plot(hr.value_counts().sort_index(), color="#333", lw=1.6, label="all")
    hp = np.floor(p.hour).astype(int)
    a.plot(hp.value_counts().sort_index(), color="#4a7ba7", lw=1.6, label="passed")
    a.set_xlabel("hour of run"); a.set_ylabel("reads started")
    a.set_title("(d) Throughput over time"); a.legend(fontsize=8)

    # (e) cumulative yield. Shows how much of your data arrived early.
    a = ax[1, 1]
    o = p.sort_values("hour")
    a.plot(o.hour, np.cumsum(o.sequence_length_template)/1e9, color="#4a7ba7")
    a.set_xlabel("hour of run"); a.set_ylabel("cumulative Gb (passed)")
    a.set_title("(e) Cumulative yield")

    # (f) end_reason: WHY each read stopped.
    #   signal_positive     = clean end, molecule finished
    #   unblock_mux_change  = pore was blocked and actively ejected
    #   device_data_error   = electronics/acquisition problem
    a = ax[1, 2]
    er = df.end_reason.value_counts()
    a.barh(range(len(er)), er.values, color="#b5714f")
    a.set_yticks(range(len(er))); a.set_yticklabels(er.index, fontsize=8)
    a.set_xscale("log"); a.invert_yaxis()
    a.set_xlabel("reads (log)"); a.set_title("(f) Why reads ended")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{out}_qc.png", dpi=130)
    plt.close(fig)
    print(f"\n  wrote {out}_qc.png")

    print("\n  end_reason breakdown:")
    for k, v in er.items():
        print(f"    {k:<26} {v:>9,}  {100*v/len(df):5.2f}%")


# ---------------------------------------------------------------- 4. MAIN
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(path).split(".")[0]
    label = os.path.basename(path)[:60]
    df = load(path)
    summarise(df, label, out)
    figures(df, label, out)
    print("\n  Done.\n")