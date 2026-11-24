import argparse
import csv
import re
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    raise SystemExit(
        "Biopython not found. Install with:\n"
        "  pip install biopython --break-system-packages\n"
        "or in your conda env:\n"
        "  conda install -n phylo biopython"
    )

STOP_CODONS = {"TAA", "TAG", "TGA"}
COVERAGE_RE = re.compile(r"^(\d+)-(\d+)/(\d+)$")


def load_contigs(fasta_path):
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(fasta_path), "fasta")}


def find_fasta_for_file(fasta_dir, file_field):
    candidate = Path(fasta_dir) / file_field
    if candidate.exists():
        return candidate
    stem = Path(file_field).stem
    matches = list(Path(fasta_dir).glob(f"{stem}*"))
    return matches[0] if matches else None


def parse_coverage(coverage_str):
    m = COVERAGE_RE.match(coverage_str.strip())
    if not m:
        return None
    ref_start, ref_end, ref_length = (int(x) for x in m.groups())
    return ref_start, ref_end, ref_length


def translate_check(nt_seq):
    n_codons = len(nt_seq) // 3
    protein = []
    stop_positions = []
    for i in range(n_codons):
        codon = nt_seq[i * 3:i * 3 + 3]
        if codon in STOP_CODONS:
            if i < n_codons - 1:
                stop_positions.append(i + 1)  # premature stop, keep scanning
                protein.append("*")
                continue
            else:
                protein.append("*")
                continue  # stop at expected end, not premature
        aa = str(Seq(codon).translate())
        protein.append(aa)
    return "".join(protein), stop_positions


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abricate", required=True, help="ABRicate .tab file for one lineage")
    ap.add_argument("--fasta_dir", required=True, help="Folder containing that lineage's assembly fastas")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--gene", default="inlA", help="Gene name to screen (default: inlA)")
    args = ap.parse_args()

    all_rows = []
    with open(args.abricate) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("GENE") != args.gene:
                continue
            all_rows.append(row)

    if not all_rows:
        raise SystemExit(f"No rows with GENE == '{args.gene}' found in {args.abricate}")

    from collections import defaultdict
    grouped = defaultdict(list)
    for row in all_rows:
        grouped[row["#FILE"]].append(row)

    rows = []
    n_isolates_with_multiple = 0
    for file_field, group in grouped.items():
        if len(group) > 1:
            n_isolates_with_multiple += 1
        best = max(group, key=lambda r: float(r["%COVERAGE"]))
        rows.append(best)

    if n_isolates_with_multiple:
        print(f"NOTE: {n_isolates_with_multiple} isolate(s) had multiple {args.gene} "
              f"hit rows; kept only the highest-%COVERAGE row for each.")

    results = []
    fasta_cache = {}
    missing_fasta = 0

    for row in rows:
        file_field = row["#FILE"]
        contig_id = row["SEQUENCE"]
        start, end = int(row["START"]), int(row["END"])
        strand = "+" if start <= end else "-"
        lo, hi = min(start, end), max(start, end)

        cov = parse_coverage(row.get("COVERAGE", ""))
        ref_start, ref_end, ref_length = cov if cov else (None, None, None)

        entry = {
            "isolate": file_field,
            "contig": contig_id,
            "start": start,
            "end": end,
            "strand": strand,
            "ref_start": ref_start,
            "ref_end": ref_end,
            "ref_length": ref_length,
            "abricate_pct_coverage": row.get("%COVERAGE", ""),
            "abricate_pct_identity": row.get("%IDENTITY", ""),
        }

        fasta_path = find_fasta_for_file(args.fasta_dir, file_field)
        if fasta_path is None:
            missing_fasta += 1
            entry.update({"frame_offset": "", "extracted_nt_len": "",
                          "frameshift_risk": "", "pmsc_detected": "",
                          "all_stop_positions": "", "early_stop_positions_leq40": "",
                          "later_stop_positions_gt40": "", "likely_real_pmsc": "",
                          "translated_protein_length": "",
                          "note": "FASTA_NOT_FOUND"})
            results.append(entry)
            continue

        if fasta_path not in fasta_cache:
            fasta_cache[fasta_path] = load_contigs(fasta_path)
        contigs = fasta_cache[fasta_path]

        if contig_id not in contigs:
            entry.update({"frame_offset": "", "extracted_nt_len": "",
                          "frameshift_risk": "", "pmsc_detected": "",
                          "all_stop_positions": "", "early_stop_positions_leq40": "",
                          "later_stop_positions_gt40": "", "likely_real_pmsc": "",
                          "translated_protein_length": "",
                          "note": f"CONTIG_{contig_id}_NOT_FOUND"})
            results.append(entry)
            continue

        seq_fwd = contigs[contig_id][lo - 1:hi]

        seq_rev = str(Seq(seq_fwd).reverse_complement())

        def apply_frame_offset(s):
            if ref_start is not None and ref_start != 1:
                off = (ref_start - 1) % 3
                return s[off:], off
            return s, 0

        nominal_seq = seq_fwd if strand == "+" else seq_rev
        flipped_seq = seq_rev if strand == "+" else seq_fwd

        nominal_seq, frame_offset = apply_frame_offset(nominal_seq)
        flipped_seq, flipped_offset = apply_frame_offset(flipped_seq)

        note = ""
        if ref_start is not None and ref_start != 1 and frame_offset != 0:
            note = f"frame_corrected_offset_{frame_offset}_ref_start_{ref_start}"
        elif ref_start is not None and ref_start != 1:
            note = f"ref_start_{ref_start}_but_in_frame"

        nominal_protein, nominal_stops = translate_check(nominal_seq)
        flipped_protein, flipped_stops = translate_check(flipped_seq)

        STRAND_FLIP_THRESHOLD = 5  # nominal must have at least this many more stops to trigger a flip
        if len(nominal_stops) - len(flipped_stops) >= STRAND_FLIP_THRESHOLD:
            seq = flipped_seq
            protein, stop_positions = flipped_protein, flipped_stops
            frame_offset = flipped_offset
            note = (note + ";STRAND_AUTO_FLIPPED_orig_had_" + str(len(nominal_stops)) + "_stops").lstrip(";")
        else:
            seq = nominal_seq
            protein, stop_positions = nominal_protein, nominal_stops

        extracted_len = len(seq)
        frameshift_risk = "YES" if extracted_len % 3 != 0 else "NO"

        early_stops = [p for p in stop_positions if p <= 40]
        later_stops = [p for p in stop_positions if p > 40]

        entry.update({
            "frame_offset": frame_offset,
            "extracted_nt_len": extracted_len,
            "frameshift_risk": frameshift_risk,
            "pmsc_detected": "YES" if stop_positions else "NO",
            "all_stop_positions": ";".join(str(p) for p in stop_positions),
            "early_stop_positions_leq40": ";".join(str(p) for p in early_stops),
            "later_stop_positions_gt40": ";".join(str(p) for p in later_stops),
            "likely_real_pmsc": "YES" if later_stops else ("SUSPECT_ARTIFACT_ONLY" if early_stops else "NO"),
            "translated_protein_length": len(protein.replace("*", "")),
            "note": note if not (frameshift_risk == "YES") else (note + ";FRAMESHIFT_RISK_VERIFY_MANUALLY").lstrip(";"),
        })
        results.append(entry)

    fieldnames = ["isolate", "contig", "start", "end", "strand",
                  "ref_start", "ref_end", "ref_length", "frame_offset",
                  "extracted_nt_len", "frameshift_risk", "pmsc_detected",
                  "all_stop_positions", "early_stop_positions_leq40",
                  "later_stop_positions_gt40", "likely_real_pmsc",
                  "translated_protein_length",
                  "abricate_pct_coverage", "abricate_pct_identity", "note"]

    with open(args.out, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    n_pmsc = sum(1 for r in results if r.get("pmsc_detected") == "YES")
    n_real = sum(1 for r in results if r.get("likely_real_pmsc") == "YES")
    n_artifact_only = sum(1 for r in results if r.get("likely_real_pmsc") == "SUSPECT_ARTIFACT_ONLY")
    n_frameshift = sum(1 for r in results if r.get("frameshift_risk") == "YES")
    print(f"Processed {len(results)} isolates.")
    print(f"  Any stop detected: {n_pmsc}")
    print(f"  Likely REAL PMSC (stop after codon 40): {n_real}")
    print(f"  Early-only (codon <=40, suspect artifact): {n_artifact_only}")
    print(f"  Frameshift-risk hits (verify manually): {n_frameshift}")
    if missing_fasta:
        print(f"  WARNING: {missing_fasta} isolates had no matching FASTA found in {args.fasta_dir}")
    print(f"Output written to {args.out}")

if __name__ == "__main__":
    main()
