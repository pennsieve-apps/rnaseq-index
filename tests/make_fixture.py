"""
A tiny synthetic genome: two chromosomes, six genes (one with two isoforms, one on the minus
strand), a GENCODE-style GTF and transcript FASTA, gzipped like the real downloads.

  python -m tests.make_fixture data/fixture            # linked layout: transcriptome/ genome/ annotation/
  python -m tests.make_fixture data/fixture --flat     # everything side by side
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")
CHROMS = {"chr1": 6000, "chr2": 4000}
# gene id, name, type, chrom, strand, transcripts: {tx id: [(start, end), ...]} (1-based, inclusive)
GENES = [
    ("ENSMUSG00000000001.1", "Alpha", "protein_coding", "chr1", "+",
     {"ENSMUST00000000001.1": [(101, 400), (601, 900)], "ENSMUST00000000002.1": [(101, 400), (1001, 1300)]}),
    ("ENSMUSG00000000002.3", "Beta", "protein_coding", "chr1", "-", {"ENSMUST00000000003.2": [(2001, 2600)]}),
    ("ENSMUSG00000000003.1", "Gamma", "lncRNA", "chr1", "+", {"ENSMUST00000000004.1": [(3001, 3200), (3401, 3800)]}),
    ("ENSMUSG00000000004.2", "Delta", "protein_coding", "chr2", "+", {"ENSMUST00000000005.1": [(201, 900)]}),
    ("ENSMUSG00000000005.1", "Epsilon", "protein_coding", "chr2", "-", {"ENSMUST00000000006.1": [(1501, 2100)]}),
    ("ENSMUSG00000000006.1", "Zeta", "snRNA", "chr2", "+", {"ENSMUST00000000007.1": [(3001, 3150)]}),
]


def genome(seed: int = 7) -> dict[str, str]:
    rng = random.Random(seed)
    return {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS.items()}


def transcripts(seqs: dict[str, str]) -> dict[str, tuple[str, str, str]]:
    """{tx id: (sequence, gene id, GENCODE header)}"""
    out = {}
    for gid, name, gtype, chrom, strand, txs in GENES:
        for i, (tid, exons) in enumerate(txs.items()):
            s = "".join(seqs[chrom][a - 1:b] for a, b in exons)
            if strand == "-":
                s = s.translate(COMP)[::-1]
            header = f"{tid}|{gid}|-|-|{name}-{201 + i}|{name}|{len(s)}|{gtype}|"
            out[tid] = (s, gid, header)
    return out


def gtf_lines() -> list[str]:
    lines = ["##description: evidence-based annotation of the mouse genome (GRCm39), version M39 (Ensembl 115)",
             "##provider: GENCODE"]
    for gid, name, gtype, chrom, strand, txs in GENES:
        start = min(a for ex in txs.values() for a, _ in ex)
        end = max(b for ex in txs.values() for _, b in ex)
        ga = f'gene_id "{gid}"; gene_type "{gtype}"; gene_name "{name}";'
        lines.append(f"{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\t{ga}")
        for tid, exons in txs.items():
            ta =f'gene_id "{gid}"; transcript_id "{tid}"; gene_type "{gtype}"; gene_name "{name}";'
            lines.append(f"{chrom}\tHAVANA\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{ta}")
            for n, (a, b) in enumerate(exons, 1):
                lines.append(f"{chrom}\tHAVANA\texon\t{a}\t{b}\t.\t{strand}\t.\t{ta} exon_number {n};")
    return lines


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def make_fixture(root: Path, flat: bool = False) -> Path:
    seqs = genome()
    txs = transcripts(seqs)
    d = {k: root if flat else root / k for k in ("transcriptome", "genome", "annotation")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    with gzip.open(d["transcriptome"] / "gencode.vM39.transcripts.fa.gz", "wt") as fh:
        for s, _, header in txs.values():
            fh.write(f">{header}\n{wrap(s)}\n")
    with gzip.open(d["genome"] / "GRCm39.primary_assembly.genome.fa.gz", "wt") as fh:
        for c, s in seqs.items():
            fh.write(f">{c} {c}\n{wrap(s)}\n")
    with gzip.open(d["annotation"] / "gencode.vM39.annotation.gtf.gz", "wt") as fh:
        fh.write("\n".join(gtf_lines()) + "\n")
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--flat", action="store_true")
    a = ap.parse_args()
    make_fixture(a.root, a.flat)
