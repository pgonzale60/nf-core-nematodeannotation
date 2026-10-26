#!/usr/bin/env python3
"""
slidr_extract_tails.py — Extract 5' soft-clipped tails and trans-splice acceptor sites from RNA-seq alignments.

Modernized component of SLIDR (Wenzel et al., BMC Bioinformatics 2021).
Extracts soft-clipped read tails from coordinate-sorted spliced BAM alignments
(compatible with STAR and HISAT2), maps genomic junction coordinates, and validates
canonical trans-splice acceptor motifs (AG / CT on reverse strand).
"""

import argparse
import gzip
import json
import os
import re
import sys
from typing import Iterator, Optional

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


REV_COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    return seq.translate(REV_COMP_TABLE)[::-1]


class FastaIndexedReader:
    """Fast indexed FASTA reader with fallback for systems without pysam."""
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.fai_path = fasta_path + ".fai"
        self._pysam_fa = None
        self._index = {}
        self._file = None

        if HAS_PYSAM:
            try:
                self._pysam_fa = pysam.FastaFile(fasta_path)
                return
            except Exception:
                self._pysam_fa = None

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

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch 0-based [start, end) subsequence from chromosome."""
        if start < 0:
            start = 0
        if start >= end:
            return ""

        if self._pysam_fa is not None:
            try:
                return self._pysam_fa.fetch(chrom, start, end).upper()
            except Exception:
                return ""

        if chrom not in self._index or self._file is None:
            return ""

        chrom_len, offset, line_bases, line_bytes = self._index[chrom]
        if end > chrom_len:
            end = chrom_len
        if start >= end:
            return ""

        # Calculate byte offsets in FASTA file
        nl_len = line_bytes - line_bases
        start_line = start // line_bases
        start_col = start % line_bases
        start_byte = offset + start_line * line_bytes + start_col

        total_bases = end - start
        # Read slightly more bytes to account for newlines
        bytes_to_read = total_bases + ((start_col + total_bases) // line_bases + 1) * nl_len
        self._file.seek(start_byte)
        chunk = self._file.read(bytes_to_read).decode("ascii", errors="replace")
        seq = "".join(chunk.split())[:total_bases]
        return seq.upper()

    def close(self):
        if self._pysam_fa is not None:
            self._pysam_fa.close()
        if self._file is not None:
            self._file.close()


def parse_sam_line(line: str) -> Optional[dict]:
    """Parse a single SAM alignment line into a record dict."""
    if line.startswith("@"):
        return None
    fields = line.strip().split("\t")
    if len(fields) < 11:
        return None

    flag = int(fields[1])
    if flag & 0x4:  # Unmapped
        return None

    return {
        "qname": fields[0],
        "flag": flag,
        "rname": fields[2],
        "pos": int(fields[3]) - 1,  # 0-based
        "mapq": int(fields[4]),
        "cigar": fields[5],
        "seq": fields[9],
        "qual": fields[10],
    }


CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHPX=])")

def parse_cigar_tuples(cigar_str: str):
    """Parse CIGAR string into list of (op_char, length)."""
    return [(op, int(length)) for length, op in CIGAR_REGEX.findall(cigar_str)]


def extract_alignment_tails(
    read_name: str,
    flag: int,
    chrom: str,
    pos: int,
    cigar_str: str,
    seq: str,
    min_tail_len: int,
    fa_reader: Optional[FastaIndexedReader],
    strandedness: int = 0
) -> Optional[dict]:
    """
    Evaluate if an alignment has a 5' soft clip representing an SL tail.
    Returns detail dict if valid, else None.
    """
    is_reverse = bool(flag & 0x10)
    is_read1 = bool(flag & 0x40)
    is_read2 = bool(flag & 0x80)

    # In stranded RNA-seq:
    # strandedness == 2 (reverse/dUTP): R2 is sense (5' end is true mRNA 5'), R1 is antisense
    # strandedness == 1 (forward): R1 is sense, R2 is antisense
    # strandedness == 0 (unstranded): treat both as sense
    is_sense = True
    if strandedness == 2:
        is_sense = is_read2 or (not is_read1 and not is_read2)  # single-end or R2
    elif strandedness == 1:
        is_sense = is_read1 or (not is_read1 and not is_read2)

    cigar_tuples = parse_cigar_tuples(cigar_str)
    if not cigar_tuples:
        return None

    clip_len = 0
    tail_seq = ""
    junction_pos = -1
    strand = "+"
    tail_type = ""

    if not is_reverse:
        # Mapped to (+) strand
        # If sense: 5' clip is at start of alignment.
        # Aligners (e.g. HISAT2) absorb the leader's 3'-terminal 'AG' into the genomic acceptor AG,
        # placing the alignment start directly at the AG. Therefore, the 0-based junction position
        # is pos + 2 (the true TSS / first exonic base), and the preceding acceptor dinucleotide is [pos, pos+2).
        if is_sense and cigar_tuples[0][0] == "S":
            clip_len = cigar_tuples[0][1]
            if clip_len >= min_tail_len:
                tail_seq = seq[:clip_len]
                junction_pos = pos + 2
                strand = "+"
                tail_type = "+5p"
        # If antisense (short insert): 3' clip revcomp can capture opposite SL
        elif (not is_sense) and cigar_tuples[-1][0] == "S":
            clip_len = cigar_tuples[-1][1]
            if clip_len >= min_tail_len:
                tail_seq = revcomp(seq[-clip_len:])
                ref_len = sum(l for op, l in cigar_tuples if op in ("M", "D", "N", "=", "X"))
                junction_pos = pos + ref_len
                strand = "+"
                tail_type = "+3p"
    else:
        # Mapped to (-) strand
        # Note: In SAM/BAM, seq of reverse read is already reverse-complemented.
        # The 5' end of the transcript on the negative strand corresponds to the 3' end (rightmost mapped)
        # of the reverse read. The aligner aligns up to the reverse-complemented AG (CT on reference + strand).
        # Therefore, the 0-based junction position is pos + ref_len - 2, and the preceding acceptor
        # dinucleotide on the reference (+) strand is [junction_pos, junction_pos + 2), whose reverse complement is AG.
        if is_sense and cigar_tuples[-1][0] == "S":
            clip_len = cigar_tuples[-1][1]
            if clip_len >= min_tail_len:
                tail_seq = revcomp(seq[-clip_len:])
                ref_len = sum(l for op, l in cigar_tuples if op in ("M", "D", "N", "=", "X"))
                junction_pos = pos + ref_len - 2
                strand = "-"
                tail_type = "-5p"
        elif (not is_sense) and cigar_tuples[0][0] == "S":
            clip_len = cigar_tuples[0][1]
            if clip_len >= min_tail_len:
                tail_seq = seq[:clip_len]
                junction_pos = pos
                strand = "-"
                tail_type = "-3p"

    if not tail_seq or junction_pos < 0:
        return None

    # Check genomic acceptor dinucleotide
    # Canonical acceptor is AG on transcript sense
    genomic_dinuc = "NN"
    is_canonical_acceptor = False
    if fa_reader is not None:
        if strand == "+":
            # Preceding dinucleotide is [junction_pos - 2, junction_pos)
            genomic_dinuc = fa_reader.fetch(chrom, junction_pos - 2, junction_pos)
            is_canonical_acceptor = (genomic_dinuc == "AG")
        else:
            # On negative strand, preceding dinucleotide in transcript sense is [junction_pos, junction_pos + 2) on reference (+) strand
            raw_dinuc = fa_reader.fetch(chrom, junction_pos, junction_pos + 2)
            genomic_dinuc = revcomp(raw_dinuc)
            is_canonical_acceptor = (genomic_dinuc == "AG")

    return {
        "read_id": read_name,
        "chrom": chrom,
        "junction_pos": junction_pos,  # 0-based
        "strand": strand,
        "tail_type": tail_type,
        "clip_len": clip_len,
        "tail_seq": tail_seq,
        "genomic_dinuc": genomic_dinuc,
        "is_canonical": is_canonical_acceptor,
    }


def stream_bam_or_sam(bam_path: str) -> Iterator[dict]:
    """Stream alignments from BAM (via pysam or samtools view) or SAM."""
    if HAS_PYSAM:
        try:
            with pysam.AlignmentFile(bam_path, "rb") as samfile:
                for read in samfile.fetch(until_eof=True):
                    if read.is_unmapped:
                        continue
                    cigar_str = read.cigarstring or ""
                    yield {
                        "qname": read.query_name,
                        "flag": read.flag,
                        "rname": read.reference_name,
                        "pos": read.reference_start,
                        "mapq": read.mapping_quality,
                        "cigar": cigar_str,
                        "seq": read.query_sequence or "",
                        "qual": read.qual or "",
                    }
            return
        except Exception:
            pass

    # Fallback: invoke samtools view -h
    import subprocess
    cmd = ["samtools", "view", "-h", bam_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        rec = parse_sam_line(line)
        if rec:
            yield rec
    proc.wait()


def main():
    parser = argparse.ArgumentParser(
        description="Extract soft-clipped 5' tails and trans-splice acceptor sites from spliced RNA-seq alignments."
    )
    parser.add_argument("--bam", required=True, help="Input coordinate-sorted BAM file")
    parser.add_argument("--fasta", required=True, help="Reference assembly FASTA file")
    parser.add_argument("--prefix", default="slidr", help="Output filename prefix")
    parser.add_argument("--min-tail-len", type=int, default=8, help="Minimum tail length (default: 8 bp)")
    parser.add_argument("--strandedness", type=str, default="0", choices=["0", "1", "2", "auto"],
                        help="Strandedness (0=unstranded, 1=forward, 2=reverse, auto=infer; default: 0)")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fa_reader = FastaIndexedReader(args.fasta)

    # Output file paths
    tails_fa_path = os.path.join(args.outdir, f"{args.prefix}.tails.fa.gz")
    sites_tsv_path = os.path.join(args.outdir, f"{args.prefix}.acceptor_sites.tsv.gz")
    stats_json_path = os.path.join(args.outdir, f"{args.prefix}.tail_stats.json")

    stranded = int(args.strandedness) if args.strandedness != "auto" else 0

    total_reads = 0
    total_clipped = 0
    canonical_ag_count = 0
    unique_junctions = set()

    with gzip.open(tails_fa_path, "wt") as fa_out, gzip.open(sites_tsv_path, "wt") as tsv_out:
        # TSV Header
        tsv_out.write("read_id\tchrom\tjunction_pos_0based\tjunction_pos_1based\tstrand\ttail_type\tclip_len\ttail_seq\tgenomic_dinuc\tis_canonical_ag\n")

        for rec in stream_bam_or_sam(args.bam):
            total_reads += 1
            if total_reads % 500000 == 0:
                print(f"[slidr_extract_tails] Processed {total_reads:,} reads | {total_clipped:,} SL tails extracted ({canonical_ag_count:,} at AG)...", file=sys.stderr)

            res = extract_alignment_tails(
                read_name=rec["qname"],
                flag=rec["flag"],
                chrom=rec["rname"],
                pos=rec["pos"],
                cigar_str=rec["cigar"],
                seq=rec["seq"],
                min_tail_len=args.min_tail_len,
                fa_reader=fa_reader,
                strandedness=stranded
            )

            if res is not None:
                total_clipped += 1
                if res["is_canonical"]:
                    canonical_ag_count += 1

                junc_id = (res["chrom"], res["junction_pos"], res["strand"])
                unique_junctions.add(junc_id)

                # Write FASTA record: header includes read_id, junction coords, clip length
                fa_header = f">{res['read_id']}|{res['chrom']}:{res['junction_pos']}:{res['strand']}|len={res['clip_len']}|dinuc={res['genomic_dinuc']}"
                fa_out.write(f"{fa_header}\n{res['tail_seq']}\n")

                # Write TSV record
                tsv_out.write(
                    f"{res['read_id']}\t{res['chrom']}\t{res['junction_pos']}\t{res['junction_pos'] + 1}\t"
                    f"{res['strand']}\t{res['tail_type']}\t{res['clip_len']}\t{res['tail_seq']}\t"
                    f"{res['genomic_dinuc']}\t{res['is_canonical']}\n"
                )

    fa_reader.close()

    metrics = {
        "total_processed_alignments": total_reads,
        "extracted_sl_tails": total_clipped,
        "canonical_ag_tails": canonical_ag_count,
        "canonical_ag_percentage": round(100.0 * canonical_ag_count / max(1, total_clipped), 2),
        "unique_junction_sites": len(unique_junctions),
        "min_tail_length": args.min_tail_len,
        "strandedness": args.strandedness,
    }

    with open(stats_json_path, "w") as jf:
        json.dump(metrics, jf, indent=2)

    print(f"[slidr_extract_tails] Complete! Processed {total_reads:,} reads. Emitted {total_clipped:,} tails across {len(unique_junctions):,} unique sites ({canonical_ag_count:,} canonical AG, {metrics['canonical_ag_percentage']}%).", file=sys.stderr)


if __name__ == "__main__":
    main()
