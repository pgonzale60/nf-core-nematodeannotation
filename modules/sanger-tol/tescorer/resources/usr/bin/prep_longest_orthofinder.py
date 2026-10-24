import os
import sys
import re

out_dir = sys.argv[1]
ref_fasta = sys.argv[2]
braker_aa = sys.argv[3]
braker_gtf = sys.argv[4]

os.makedirs(out_dir, exist_ok=True)

# -------------------------------------------------------------
# 1. BRAKER3: Filter to Longest Isoform per Gene
# -------------------------------------------------------------
print("Mapping BRAKER3 transcripts to genes...")
tx2gene = {}
with open(braker_gtf, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] in ('CDS', 'exon', 'transcript', 'mRNA'):
            attr = parts[8]
            m_tx = re.search(r'transcript_id "([^"]+)"', attr)
            m_gn = re.search(r'gene_id "([^"]+)"', attr)
            if m_tx and m_gn:
                tx2gene[m_tx.group(1)] = m_gn.group(1)

print(f"Mapped {len(tx2gene)} BRAKER3 transcripts to {len(set(tx2gene.values()))} genes.")

# Load BRAKER3 proteins
print("Selecting longest protein isoform per BRAKER3 gene...")
gene_longest = {} # gene -> (tx_id, seq)
current_tx = None
current_seq = []

with open(braker_aa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if current_tx:
                seq_str = "".join(current_seq).replace('*', '')
                gn = tx2gene.get(current_tx, current_tx.split('.')[0])
                if gn not in gene_longest or len(seq_str) > len(gene_longest[gn][1]):
                    gene_longest[gn] = (current_tx, seq_str)
            current_tx = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_tx:
        seq_str = "".join(current_seq).replace('*', '')
        gn = tx2gene.get(current_tx, current_tx.split('.')[0])
        if gn not in gene_longest or len(seq_str) > len(gene_longest[gn][1]):
            gene_longest[gn] = (current_tx, seq_str)

out_braker = os.path.join(out_dir, 'nxAuaRhod1_APS4.fa')
with open(out_braker, 'w') as f:
    for gn, (tx_id, seq_str) in gene_longest.items():
        f.write(f">{tx_id}\n{seq_str}\n")
print(f"Wrote {len(gene_longest)} representative BRAKER3 proteins to {out_braker}.")

# -------------------------------------------------------------
# 2. Reference Proteomes: Longest Isoform per Gene
# -------------------------------------------------------------
species_names = ['Caenorhabditis elegans', 'Pristionchus pacificus', 'Brugia malayi']
species_files = {
    'Caenorhabditis elegans': os.path.join(out_dir, 'c_elegans.fa'),
    'Pristionchus pacificus': os.path.join(out_dir, 'p_pacificus.fa'),
    'Brugia malayi': os.path.join(out_dir, 'b_malayi.fa')
}

ref_longest = {sp: {} for sp in species_names} # sp -> {gene_id -> (seq_id, seq)}

print("Parsing reference nematode proteomes for longest gene isoform...")
current_sp = None
current_id = None
current_gn = None
current_seq = []

with open(ref_fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if current_sp and current_id:
                seq_str = "".join(current_seq).replace('*', '')
                if current_gn not in ref_longest[current_sp] or len(seq_str) > len(ref_longest[current_sp][current_gn][1]):
                    ref_longest[current_sp][current_gn] = (current_id, seq_str)
            
            current_sp = None
            current_id = None
            current_gn = None
            current_seq = []
            
            for sp in species_names:
                if f"OS={sp}" in line:
                    current_sp = sp
                    current_id = line[1:].strip().split()[0].replace('|', '_')
                    m_gn = re.search(r'GN=([^\s]+)', line)
                    if m_gn:
                        current_gn = m_gn.group(1)
                    else:
                        current_gn = current_id
                    break
        else:
            if current_sp:
                current_seq.append(line.strip())

    if current_sp and current_id:
        seq_str = "".join(current_seq).replace('*', '')
        if current_gn not in ref_longest[current_sp] or len(seq_str) > len(ref_longest[current_sp][current_gn][1]):
            ref_longest[current_sp][current_gn] = (current_id, seq_str)

for sp, out_path in species_files.items():
    with open(out_path, 'w') as f:
        for gn, (seq_id, seq_str) in ref_longest[sp].items():
            f.write(f">{seq_id}\n{seq_str}\n")
    print(f"Wrote {len(ref_longest[sp])} representative proteins for {sp} to {out_path}.")

print("Input preparation complete!")
