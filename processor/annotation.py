"""
Transcript-to-gene table and gene table from a GTF (GENCODE or Ensembl), keyed by the
transcript ids salmon indexes (the FASTA's ids, versions kept as there).
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')
GENCODE_DESC_RE = re.compile(r"annotation of the (\w+) genome \((\S+)\), version (\w+)")
ENSEMBL_NAME_RE = re.compile(r"^([A-Z][a-z]+_[a-z]+)\.([^.]+(?:\.[^.]+)?)\.(\d+)\.(?:chr\.)?gtf")
GENCODE_ORGANISMS = {"human": "Homo sapiens", "mouse": "Mus musculus"}


def open_text(path: Path):
    """A FASTA/GTF, gzipped or plain, as text (decided by the gzip magic bytes, not the name)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="ascii", errors="replace") if gz else open(path, encoding="ascii", errors="replace")


def fasta_ids(path: Path, gencode: bool) -> list[str]:
    """The record ids of a FASTA in file order: the first word of each header, cut at '|' for GENCODE."""
    ids = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0] if len(line) > 2 else ""
                ids.append(name.split("|", 1)[0] if gencode else name)
    return ids


def is_gencode_fasta(path: Path) -> bool:
    """GENCODE transcript FASTA headers are '|'-separated (ENST...|ENSG...|...|name|length|biotype|)."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                return "|" in line.split(maxsplit=1)[0]
    return False


def _versioned(attrs: dict[str, str], key: str) -> str:
    """An id with its version: GENCODE writes it into the id, Ensembl in a separate <kind>_version attribute."""
    ident = attrs.get(f"{key}_id", "")
    version = attrs.get(f"{key}_version")
    return f"{ident}.{version}" if version and "." not in ident else ident


def _unversioned(ident: str) -> str:
    head, dot, tail = ident.rpartition(".")
    return head if dot and tail.isdigit() else ident


@dataclass
class Gene:
    name: str
    type: str
    chrom: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    def length(self) -> int:
        """Exonic length: the union of all exons of all transcripts (overlaps counted once)."""
        total, end = 0, 0
        for s, e in sorted(self.exons):
            if e > end:
                total += e - max(s - 1, end)
                end = e
        return total


@dataclass
class Annotation:
    tx_gene: dict[str, str]
    genes: dict[str, Gene]
    header: list[str]


def read_gtf(path: Path) -> Annotation:
    """Transcripts and genes of a GTF; gene length from the exon lines."""
    tx_gene: dict[str, str] = {}
    genes: dict[str, Gene] = {}
    header: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if len(header) < 20:
                    header.append(line.strip())
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in ("transcript", "exon"):
                continue
            attrs = dict(ATTR_RE.findall(cols[8]))
            gene_id = _versioned(attrs, "gene")
            if not gene_id:
                continue
            g = genes.get(gene_id)
            if g is None:
                g = genes[gene_id] = Gene(attrs.get("gene_name", gene_id),
                                          attrs.get("gene_type") or attrs.get("gene_biotype", ""), cols[0])
            if cols[2] == "transcript":
                tx_gene[_versioned(attrs, "transcript")] = gene_id
            else:
                g.exons.append((int(cols[3]), int(cols[4])))
    if not tx_gene:
        raise ValueError(f"{path.name}: no transcript lines with transcript_id and gene_id (is it a GTF?)")
    return Annotation(tx_gene, genes, header)


def tx2gene_rows(ids: list[str], ann: Annotation) -> tuple[list[dict], list[str]]:
    """
    One row per FASTA transcript: matched to the GTF by id (with version), else by id without
    version. Transcripts the GTF does not have become their own gene (as salmon and tximport
    treat them), and are returned as the second value.
    """
    loose = {}
    for tx, gene in ann.tx_gene.items():
        loose.setdefault(_unversioned(tx), gene)
    rows, missing = [], []
    for tx in dict.fromkeys(ids):
        gene_id = ann.tx_gene.get(tx) or loose.get(_unversioned(tx))
        if gene_id is None:
            missing.append(tx)
            rows.append({"transcript_id": tx, "gene_id": tx, "gene_name": tx, "gene_type": ""})
            continue
        g = ann.genes[gene_id]
        rows.append({"transcript_id": tx, "gene_id": gene_id, "gene_name": g.name, "gene_type": g.type})
    return rows, missing


def gene_rows(tx_rows: list[dict], ann: Annotation) -> list[dict]:
    """genes.csv: every gene with at least one indexed transcript."""
    out = {}
    for r in tx_rows:
        gid = r["gene_id"]
        if gid in out:
            continue
        g = ann.genes.get(gid)
        out[gid] = ({"gene_id": gid, "gene_name": g.name, "gene_type": g.type, "chrom": g.chrom, "length": g.length()}
                    if g else {"gene_id": gid, "gene_name": gid, "gene_type": "", "chrom": "", "length": ""})
    return list(out.values())


def describe(ann: Annotation, gtf_name: str) -> tuple[str, str, str]:
    """(organism, source, assembly) from a GENCODE GTF header or an Ensembl GTF file name; ('', gtf_name, '') otherwise."""
    for line in ann.header:
        m = GENCODE_DESC_RE.search(line)
        if m:
            return GENCODE_ORGANISMS.get(m.group(1), m.group(1)), f"GENCODE {m.group(3)}", m.group(2)
    m = ENSEMBL_NAME_RE.match(gtf_name)
    if m:
        return m.group(1).replace("_", " "), f"Ensembl {m.group(3)}", m.group(2)
    return "", gtf_name, ""
