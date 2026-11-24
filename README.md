# inlA_pmsc_check.py
inlA pmsc checker for Listeria monocytogenes

USAGE (run separately per lineage):
    python inlA_pmsc_check.py \
        --abricate virulence_results.tab \
        --fasta_dir/path/to/fastas/ \
        --out inlA_PMSC.tsv

REQUIREMENTS:
Requires Python ≥ 3.8 and Biopython (v1.85 used)
OUTPUT COLUMNS (tab-separated, one row per isolate; for isolates with several inlA hits, the hit with the highest ABRicate %COVERAGE is kept):
isolate, contig, start, end, strand: ABRicate file name, contig and alignment coordinates; strand is inferred from start > end.
ref_start, ref_end, ref_length: the aligned span of the inlA reference, parsed from the ABRicate COVERAGE field.
frame_offset: the number of bases trimmed from the 5′ end to restore the reading frame when the alignment does not start at reference position 1.
extracted_nt_len: the length of the extracted region after trimming.
frameshift_risk: YES if extracted_nt_len is not a multiple of three.
pmsc_detected: YES if any in-frame stop codon occurs before the final codon.
all_stop_positions: codon positions of all premature stops, 1-based from the start of the extracted region.
early_stop_positions_leq40: premature stops at codon ≤ 40, treated as possible alignment-edge artefacts.
later_stop_positions_gt40: premature stops at codon > 40.
likely_real_pmsc: YES if any stop is after codon 40; SUSPECT_ARTIFACT_ONLY if all stops are at codon ≤ 40; otherwise NO.
translated_protein_length: the number of translated codons excluding stops.
abricate_pct_coverage, abricate_pct_identity: copied from ABRicate.
note: flags for frame correction, automatic strand reversal (applied when the nominal strand has at least five more stops than the reverse), frameshift risk, or a missing FASTA or contig.
