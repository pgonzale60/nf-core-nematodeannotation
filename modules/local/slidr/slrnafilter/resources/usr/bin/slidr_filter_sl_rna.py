#!/usr/bin/env python3
"""
slidr_filter_sl_rna.py — Identify genomic SL RNA candidate genes by scanning splice donor and Sm-binding motifs.

Modernized component of SLIDR (Wenzel et al., BMC Bioinformatics 2021).
Scans genomic regions flanking tail centroid BLASTN hits for:
  1. Canonical 5' splice donor site (GT)
  2. Sm ribonucleoprotein binding site motif (.{20,60}AT{4,6}G)
Outputs candidate SL RNA gene loci and sequences for secondary structure evaluation (RNAfold).
"""

import argparse
import gzip
import os
import re
import sys
from typing import Dict

REV_COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq: str) -> str:
    return seq.translate(REV_COMP_TABLE)[::-1]


class FastaIndexedReader:
    """Fast indexed FASTA reader."""
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.fai_path = fasta_path + ".fai"
        self._index = {}
        self._file = None

        if os.path.exists(self.fai_path):
            with open(self.fai_path, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        chrom = parts[0]
                        length = int(parts[1])
                        offset = int(parts[2])
                        line_bases = int(parts[3])
                        line_bytes = int(parts[4])
                        self._index[chrom] = (length, offset, line_bases, line_bytes)
            self._file = open(fasta_path, "rb")

        if not self._index:
            # No usable .fai beside the FASTA. Build one instead of proceeding
            # with an empty index: fetch() would then return "" for every request
            # and the caller would silently analyse nothing. That failure mode
            # cost a production run its entire candidate set -- 30,341,155
            # alignments in, zero candidates out, and no error reported anywhere --
            # because the task directory held only a symlink to the FASTA, while
            # the .fai next to the real file was never staged.
            self._build_index()
        if not self._index or self._file is None:
            raise RuntimeError(
                f"FastaIndexedReader could not index {self.fasta_path}: no readable "
                f"{self.fai_path} and the FASTA itself could not be scanned"
            )

    def _build_index(self) -> None:
        """
        Construct an fai-style index by scanning the FASTA once, and cache it
        beside the FASTA when that directory is writable (harmless when not).

        Line widths are taken from each record's first sequence line, which is
        what an fai assumes.
        """
        index = {}
        try:
            with open(self.fasta_path, "rb") as fh:
                chrom = None
                length = 0
                offset = 0
                line_bases = 0
                line_bytes = 0
                first_of_record = True
                while True:
                    line = fh.readline()
                    if not line:
                        break
                    if line.startswith(b">"):
                        if chrom is not None:
                            index[chrom] = (length, offset, line_bases, line_bytes)
                        chrom = line[1:].split()[0].decode("ascii", "replace")
                        length = 0
                        offset = fh.tell()
                        line_bases = 0
                        line_bytes = 0
                        first_of_record = True
                        continue
                    if chrom is None:
                        continue
                    seq = line.rstrip(b"\n").rstrip(b"\r")
                    if first_of_record and seq:
                        line_bases = len(seq)
                        line_bytes = len(line)
                        first_of_record = False
                    length += len(seq)
                if chrom is not None:
                    index[chrom] = (length, offset, line_bases, line_bytes)
            self._file = open(self.fasta_path, "rb")
        except OSError:
            return

        self._index = index
        try:
            with open(self.fai_path, "w") as out:
                for name, (ln, off, lb, nlb) in index.items():
                    out.write(f"{name}\t{ln}\t{off}\t{lb}\t{nlb}\n")
        except OSError:
            pass

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if start < 0:
            start = 0
        if chrom not in self._index or self._file is None:
            return ""

        chrom_len, offset, line_bases, line_bytes = self._index[chrom]
        if end > chrom_len:
            end = chrom_len
        if start >= end:
            return ""

        nl_len = line_bytes - line_bases
        start_line = start // line_bases
        start_col = start % line_bases
        start_byte = offset + start_line * line_bytes + start_col

        total_bases = end - start
        bytes_to_read = total_bases + ((start_col + total_bases) // line_bases + 1) * nl_len
        self._file.seek(start_byte)
        chunk = self._file.read(bytes_to_read).decode("ascii", errors="replace")
        seq = "".join(chunk.split())[:total_bases]
        return seq.upper()

    def close(self):
        if self._file is not None:
            self._file.close()


def read_fasta_dict(fasta_path: str) -> Dict[str, str]:
    """Read FASTA into a dict {seq_id: seq}."""
    res = {}
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        curr_id = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_id:
                    res[curr_id] = "".join(seq_parts)
                curr_id = line[1:].split()[0].split(";")[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if curr_id:
            res[curr_id] = "".join(seq_parts)
    return res


def main():
    parser = argparse.ArgumentParser(
        description="Filter tail centroids for SL RNA genomic loci (donor GT + Sm motif)."
    )
    parser.add_argument("--blast", required=True, help="BLASTN tabular output file (outfmt 6)")
    parser.add_argument("--centroids", required=True, help="Tail centroids FASTA file")
    parser.add_argument("--fasta", required=True, help="Reference assembly FASTA file")
    parser.add_argument("--prefix", default="slidr", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--donor", default="GT", help="Splice donor regex (default: GT)")
    parser.add_argument("--sm", default=r".{20,60}AT{4,6}G", help="Sm-binding site motif regex (default: .{20,60}AT{4,6}G)")
    parser.add_argument("--max-rna-len", type=int, default=80, help="Max SL RNA length excluding SL (default: 80 bp)")
    parser.add_argument("--max-overlap", type=int, default=10, help="Max SL overlap with donor site (default: 10 bp)")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fa_reader = FastaIndexedReader(args.fasta)
    centroids = read_fasta_dict(args.centroids)

    pattern = re.compile(rf"({args.donor})({args.sm})", re.IGNORECASE)

    out_fa_path = os.path.join(args.outdir, f"{args.prefix}.candidate_sl_rnas.fa")
    out_tsv_path = os.path.join(args.outdir, f"{args.prefix}.sl_rna_candidates.tsv")
    out_gff_path = os.path.join(args.outdir, f"{args.prefix}.sl_rna_candidates.gff3")

    candidates_count = 0
    blast_records_count = 0

    opener = gzip.open if args.blast.endswith(".gz") else open

    with opener(args.blast, "rt") as bf, \
         open(out_fa_path, "w") as fa_out, \
         open(out_tsv_path, "w") as tsv_out, \
         open(out_gff_path, "w") as gff_out:

        tsv_out.write("centroid_id\tchrom\tsl_start\tsl_end\tstrand\tsl_seq\tdonor_motif\tsm_motif\tsl_rna_seq\n")
        gff_out.write("##gff-version 3\n")

        for line in bf:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 14:
                continue

            blast_records_count += 1

            # Fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qseq, qlen
            qseqid = parts[0].split(";")[0]
            chrom = parts[1]
            qstart = int(parts[6])
            qend = int(parts[7])
            sstart = int(parts[8])
            send = int(parts[9])
            qseq = parts[12]
            qlen = int(parts[13])

            # In SLIDR, we require that the 3' end of the centroid tail aligns: qend == qlen
            if qend != qlen:
                continue

            is_plus_strand = (sstart < send)
            strand = "+" if is_plus_strand else "-"

            sl_seq = centroids.get(qseqid, qseq)

            # Query flanking genomic window downstream of alignment
            window_len = args.max_overlap + args.max_rna_len
            if is_plus_strand:
                # sstart and send are 1-based coordinates
                # Search downstream: from send onwards
                window_start = send  # 0-based is send
                window_end = send + window_len
                flank_seq = fa_reader.fetch(chrom, window_start, window_end)
            else:
                # On negative strand, downstream of 3' end of transcript is upstream on reference (+) strand
                window_start = max(0, send - 1 - window_len)
                window_end = send - 1
                raw_flank = fa_reader.fetch(chrom, window_start, window_end)
                flank_seq = revcomp(raw_flank)

            match = pattern.search(flank_seq)
            if match:
                match_start = match.start()
                # Check overlap constraint: donor site must occur within max_overlap bases
                if match_start <= args.max_overlap:
                    donor_seq = match.group(1).upper()
                    sm_seq = match.group(2).upper()
                    full_rna_trailer = flank_seq[:match.end()].upper()
                    full_sl_rna = (sl_seq + full_rna_trailer).upper()

                    candidates_count += 1
                    cid = f"{qseqid}_loc{candidates_count}"

                    # Compute genomic gene coordinates (1-based closed)
                    if is_plus_strand:
                        gene_start = sstart
                        gene_end = send + match.end()
                    else:
                        gene_start = send - match.end()
                        gene_end = sstart

                    fa_out.write(f">{cid} {chrom}:{gene_start}-{gene_end}({strand}) sl={qseqid}\n{full_sl_rna}\n")
                    tsv_out.write(f"{qseqid}\t{chrom}\t{gene_start}\t{gene_end}\t{strand}\t{sl_seq}\t{donor_seq}\t{sm_seq}\t{full_sl_rna}\n")

                    gff_out.write(
                        f"{chrom}\tSLIDR\tSL_RNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
                        f"ID={cid};Parent_SL={qseqid};Donor={donor_seq};Sm_Motif={sm_seq}\n"
                    )

    fa_reader.close()
    print(f"[slidr_filter_sl_rna] Processed {blast_records_count:,} alignments -> Identified {candidates_count:,} candidate SL RNA loci.", file=sys.stderr)


if __name__ == "__main__":
    main()
