#!/usr/bin/env python3
"""
BINF 6350 Lab Quest - Stage B1: locate a gene family in the reference annotation.

Input: a tab-separated gene list with a header:
    symbol    description_regex    note
  symbol            official gene symbol to look for (e.g. HBB)
  description_regex regex matched against the GTF 'description' attribute,
                    to catch genes annotated only as LOC numbers
                    (e.g. "hemoglobin subunit beta"); may be empty
  note              free text (why it is on the list); may be empty

For each list entry the script searches the GTF gene records by:
  1. symbol   : gene_id / gene attribute equals the symbol (case-insensitive)
  2. description : description matches the regex (only reported if 1 found nothing,
                   or as extra LOC paralogs)
and reports coordinates, transcript IDs and transcript lengths (sum of exons).

Outputs:
  <prefix>_gene_table.tsv     one row per located annotation gene
  <prefix>_transcripts.tsv    one row per transcript (ID, length, exons)
  <prefix>_not_found.tsv      list entries with no match
  and a printed gate summary.

Usage:
  python b1_find_genes.py <genomic.gtf> <gene_list.tsv> <out_prefix>
"""
import csv
import os
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def attrs(field):
    return {k: v for k, v in ATTR_RE.findall(field)}


def parse_gtf(path):
    genes = {}                                  # gene_id -> info
    tx = defaultdict(lambda: {"exons": 0, "length": 0, "gene_id": None,
                              "biotype": ""})
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                continue
            seqid, _, feat, start, end, _, strand, _, a = f
            if feat == "gene":
                d = attrs(a)
                gid = d.get("gene_id")
                genes[gid] = {
                    "gene_id": gid,
                    "gene": d.get("gene", gid),
                    "description": d.get("description", ""),
                    "biotype": d.get("gene_biotype", ""),
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "strand": strand,
                }
            elif feat == "exon":
                d = attrs(a)
                tid = d.get("transcript_id")
                if not tid:
                    continue
                t = tx[tid]
                t["exons"] += 1
                t["length"] += int(end) - int(start) + 1
                t["gene_id"] = d.get("gene_id")
                t["biotype"] = d.get("transcript_biotype", t["biotype"])
    by_gene = defaultdict(list)
    for tid, t in tx.items():
        by_gene[t["gene_id"]].append((tid, t["length"], t["exons"]))
    return genes, by_gene


def read_list(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")),
                                delimiter="\t"):
            r = {k.strip(): (v or "").strip() for k, v in r.items() if k}
            if r.get("symbol"):
                rows.append(r)
    return rows


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    gtf, glist, out = sys.argv[1:]
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    genes, by_gene = parse_gtf(gtf)
    n_loc = sum(1 for g in genes.values() if g["gene_id"].startswith("LOC"))
    print(f"GTF: {len(genes):,} genes ({n_loc:,} named only LOC...)")

    entries = read_list(glist)
    print(f"gene list: {len(entries)} entries\n")

    by_symbol = defaultdict(list)
    for g in genes.values():
        by_symbol[g["gene_id"].upper()].append(g)
        if g["gene"].upper() != g["gene_id"].upper():
            by_symbol[g["gene"].upper()].append(g)

    gene_rows, tx_rows, missing = [], [], []
    n_sym = n_desc = 0
    for e in entries:
        sym = e["symbol"]
        rx = e.get("description_regex", "")
        hits = [(g, "symbol") for g in by_symbol.get(sym.upper(), [])]
        seen = {g["gene_id"] for g, _ in hits}
        if rx:
            pat = re.compile(rx, re.I)
            for g in genes.values():
                if g["gene_id"] not in seen and pat.search(g["description"]):
                    hits.append((g, "description"))
                    seen.add(g["gene_id"])
        if not hits:
            missing.append(e)
            print(f"  NOT FOUND  {sym}")
            continue
        how = "symbol" if any(m == "symbol" for _, m in hits) else "description"
        if how == "symbol":
            n_sym += 1
        else:
            n_desc += 1
        print(f"  found      {sym:<12} by {how:<11} "
              f"-> {', '.join(g['gene_id'] for g, _ in hits)}")
        for g, m in hits:
            txs = sorted(by_gene.get(g["gene_id"], []))
            lens = [l for _, l, _ in txs]
            gene_rows.append({
                "list_symbol": sym, "matched_by": m,
                "gene_id": g["gene_id"], "description": g["description"],
                "biotype": g["biotype"], "seqid": g["seqid"],
                "start": g["start"], "end": g["end"], "strand": g["strand"],
                "n_transcripts": len(txs),
                "transcript_ids": ";".join(t for t, _, _ in txs),
                "tx_len_min": min(lens) if lens else "",
                "tx_len_max": max(lens) if lens else "",
                "note": e.get("note", ""),
            })
            for t, l, x in txs:
                tx_rows.append({"list_symbol": sym, "gene_id": g["gene_id"],
                                "transcript_id": t, "length_nt": l, "exons": x})

    def write(path, rows, cols):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write(f"{out}_gene_table.tsv", gene_rows, list(gene_rows[0]) if gene_rows else ["list_symbol"])
    write(f"{out}_transcripts.tsv", tx_rows,
          ["list_symbol", "gene_id", "transcript_id", "length_nt", "exons"])
    write(f"{out}_not_found.tsv", missing,
          list(missing[0]) if missing else ["symbol", "description_regex", "note"])

    print("\n===== GATE B1 =====")
    print(f"  genes on final list            {len(entries)}")
    print(f"  located in reference annotation {n_sym + n_desc}"
          f"   (by symbol {n_sym}, only by description/LOC {n_desc})")
    print(f"  not found                      {len(missing)}")
    print(f"\nWrote {out}_gene_table.tsv, _transcripts.tsv, _not_found.tsv")


if __name__ == "__main__":
    main()
