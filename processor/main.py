"""
rnaseq-index: a salmon index plus gene tables for one genome, built once and kept as a
workspace reference for rnaseq-quant.

  OUTPUT_DIR/salmon/          the salmon index (genome as decoys when a genome FASTA is given)
  OUTPUT_DIR/tx2gene.csv      transcript_id, gene_id, gene_name, gene_type
  OUTPUT_DIR/genes.csv        gene_id, gene_name, gene_type, chrom, length
  OUTPUT_DIR/index.json       organism, source, inputs (file, sha256), salmon version, decoys
  OUTPUT_DIR/provenance.json
"""
from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from processor import annotation as A
from processor import resources

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("rnaseq-index")

APP_NAME = "rnaseq-index"
REPO_ROOT = Path(__file__).resolve().parent.parent
SALMON = os.environ.get("SALMON_BIN", "salmon")
FASTA_RE = re.compile(r"\.(fa|fasta|fna)(\.gz)?$", re.IGNORECASE)
GTF_RE = re.compile(r"\.gtf(\.gz)?$", re.IGNORECASE)
GENOME_NAME_RE = re.compile(r"genome|primary_assembly|toplevel|dna_sm|dna_rm|\.dna\.", re.IGNORECASE)
PORTS = ("transcriptome", "genome", "annotation")
MAX_MISSING_FRACTION = 0.5
CITATIONS = {
    "salmon": {"doi": "10.1038/nmeth.4197", "pmid": "28263959", "text": "Patro et al. Nature Methods 2017"},
    "decoys": {"doi": "10.1186/s13059-020-02151-8", "pmid": "32894187", "text": "Srivastava et al. Genome Biology 2020"},
    "gencode": {"doi": "10.1093/nar/gkae1078", "pmid": "39565199", "text": "Mudge et al. Nucleic Acids Research 2025"},
}
DEFAULTS = {"KMER_LEN": "31", "ORGANISM": "", "SOURCE": "", "THREADS": "0"}


class InputError(Exception):
    """No usable input: exit 1."""


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) or v).strip() for k, v in DEFAULTS.items()}
    try:
        k, threads = int(raw["KMER_LEN"]), int(raw["THREADS"])
    except ValueError as exc:
        raise ValueError(f"KMER_LEN and THREADS must be integers: {exc}") from exc
    if not (15 <= k <= 31 and k % 2 == 1):
        raise ValueError(f"KMER_LEN must be odd and between 15 and 31, got {k}")
    if threads < 0:
        raise ValueError("THREADS must be >= 0")
    for name in ("ORGANISM", "SOURCE"):
        if len(raw[name]) > 200 or any(c in raw[name] for c in "\n\r\t"):
            raise ValueError(f"{name} must be one line of at most 200 characters")
    return {"KMER_LEN": k, "ORGANISM": raw["ORGANISM"], "SOURCE": raw["SOURCE"], "THREADS": threads}


def walk_files(root: Path) -> list[Path]:
    """Every file under root (symlinks followed, loops and hidden entries skipped), sorted."""
    out, seen = [], set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        out += [Path(dirpath, f) for f in sorted(filenames) if not f.startswith(".") and Path(dirpath, f).is_file()]
    return out


def find_inputs(root: Path) -> dict[str, Path | None]:
    """
    transcriptome, genome (optional) and annotation: from INPUT_DIR/<port>/ when linked, else by
    name anywhere under INPUT_DIR (one GTF; one FASTA = the transcriptome; with two, the one named
    like a genome (genome, primary_assembly, toplevel, .dna.) is the genome).
    """
    files = walk_files(root)
    by_port: dict[str, list[Path]] = {p: [] for p in PORTS}
    loose_fasta, loose_gtf = [], []
    for f in files:
        top = f.relative_to(root).parts[0] if len(f.relative_to(root).parts) > 1 else ""
        if top in ("transcriptome", "genome") and FASTA_RE.search(f.name):
            by_port[top].append(f)
        elif top == "annotation" and GTF_RE.search(f.name):
            by_port["annotation"].append(f)
        elif FASTA_RE.search(f.name):
            loose_fasta.append(f)
        elif GTF_RE.search(f.name):
            loose_gtf.append(f)
    if not by_port["annotation"]:
        by_port["annotation"] = loose_gtf
    if not by_port["transcriptome"] and not by_port["genome"]:
        genomes = [f for f in loose_fasta if GENOME_NAME_RE.search(f.name)]
        others = [f for f in loose_fasta if f not in genomes]
        if len(loose_fasta) == 1:
            by_port["transcriptome"] = loose_fasta
        elif len(others) == 1 and len(genomes) == 1:
            by_port["transcriptome"], by_port["genome"] = others, genomes
        else:
            by_port["transcriptome"] = loose_fasta
    elif not by_port["transcriptome"]:
        by_port["transcriptome"] = loose_fasta
    out: dict[str, Path | None] = {}
    for port in PORTS:
        found = by_port[port]
        if len(found) > 1:
            raise InputError(f"{port}: expected one file, found {len(found)} "
                             f"({', '.join(str(f.relative_to(root)) for f in found[:5])}); link only one")
        out[port] = found[0] if found else None
    if out["transcriptome"] is None:
        raise InputError("no transcriptome FASTA (.fa/.fasta/.fna, optionally .gz) under INPUT_DIR")
    if out["annotation"] is None:
        raise InputError("no annotation GTF (.gtf or .gtf.gz) under INPUT_DIR")
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def salmon_version() -> str:
    out = subprocess.run([SALMON, "--version"], capture_output=True, text=True, check=True).stdout
    return out.split()[-1]


def write_csv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def build(input_dir: Path, output_dir: Path, params: dict, timings: dict) -> dict:
    t0 = time.monotonic()
    found = find_inputs(input_dir)
    rel = {k: (str(v.relative_to(input_dir)) if v else None) for k, v in found.items()}
    log.info("transcriptome %s, genome %s, annotation %s", rel["transcriptome"], rel["genome"] or "none (no decoys)",
             rel["annotation"])

    tx, genome, gtf = found["transcriptome"], found["genome"], found["annotation"]
    gencode = A.is_gencode_fasta(tx)
    ids = A.fasta_ids(tx, gencode)
    if not ids:
        raise InputError(f"{rel['transcriptome']} has no FASTA records")
    n_dup = len(ids) - len(set(ids))
    if n_dup:
        raise InputError(f"{rel['transcriptome']}: {n_dup} transcript ids occur more than once")
    ann = A.read_gtf(gtf)
    tx_rows, missing = A.tx2gene_rows(ids, ann)
    timings["annotation"] = round(time.monotonic() - t0, 1)
    if len(missing) > MAX_MISSING_FRACTION * len(ids):
        raise InputError(f"{len(missing)} of {len(ids)} transcripts (e.g. {missing[0]}) are not in {gtf.name}: "
                         "the annotation does not belong to this transcriptome (same source and release?)")
    if missing:
        log.warning("%d of %d transcripts are not in the annotation; each is its own gene (e.g. %s)",
                    len(missing), len(ids), ", ".join(missing[:3]))
    genes = A.gene_rows(tx_rows, ann)
    log.info("%d transcripts, %d genes, GENCODE headers: %s (%.1fs)", len(tx_rows), len(genes), gencode,
             timings["annotation"])

    build_dir = Path("/tmp") / "rnaseq-index-build"
    shutil.rmtree(build_dir, ignore_errors=True)
    build_dir.mkdir(parents=True)
    n_threads = resources.threads(params["THREADS"])
    cmd = [SALMON, "index", "-t", str(tx), "-i", str(output_dir / "salmon"), "-k", str(params["KMER_LEN"]),
           "-p", str(n_threads), "--tmpdir", str(build_dir)]
    if gencode:
        cmd.append("--gencode")
    decoys: list[str] = []
    if genome is not None:
        t1 = time.monotonic()
        decoys = A.fasta_ids(genome, False)
        clash = set(decoys) & set(ids)
        if not decoys or clash:
            raise InputError(f"{rel['genome']}: " + ("no FASTA records" if not decoys else
                             f"{len(clash)} sequence names are also transcript ids (e.g. {sorted(clash)[0]})"))
        (build_dir / "decoys.txt").write_text("\n".join(decoys) + "\n")
        cmd[3:4] = [str(tx), str(genome)]          # decoy records must come last
        cmd += ["-d", str(build_dir / "decoys.txt")]
        timings["decoy_names"] = round(time.monotonic() - t1, 1)
        log.info("genome: %d sequences used as decoys", len(decoys))
    limit = resources.memory_limit()
    if limit:
        # sshash's external sort: half the container's memory, the rest is for the other build stages
        cmd += ["--ramLimit", str(max(1, limit // 2**31))]
    log.info("running: %s", " ".join(cmd))
    t2 = time.monotonic()
    proc = subprocess.run(cmd, stdout=sys.stdout, stderr=subprocess.STDOUT, check=False,
                          env={**os.environ, "NO_COLOR": "1"})
    timings["salmon_index"] = round(time.monotonic() - t2, 1)
    shutil.rmtree(build_dir, ignore_errors=True)
    if proc.returncode != 0:
        raise InputError(f"salmon index failed with exit code {proc.returncode} (see the log above)")
    log.info("salmon index built in %.1fs", timings["salmon_index"])

    write_csv(output_dir / "tx2gene.csv", ["transcript_id", "gene_id", "gene_name", "gene_type"], tx_rows)
    write_csv(output_dir / "genes.csv", ["gene_id", "gene_name", "gene_type", "chrom", "length"], genes)
    organism, source, assembly = A.describe(ann, gtf.name)
    t3 = time.monotonic()
    inputs = {k: {"file": rel[k], "sha256": sha256(v)} for k, v in found.items() if v is not None}
    timings["checksums"] = round(time.monotonic() - t3, 1)
    index_doc = {
        "organism": params["ORGANISM"] or organism or "unknown",
        "source": params["SOURCE"] or source,
        "assembly": assembly,
        "salmon_version": salmon_version(),
        "kmer_len": params["KMER_LEN"],
        "decoys": genome is not None,
        "n_decoy_sequences": len(decoys),
        "gencode_headers": gencode,
        "n_transcripts": len(tx_rows),
        "n_genes": len(genes),
        "n_transcripts_not_in_annotation": len(missing),
        "inputs": inputs,
        "files": {"salmon": "salmon/", "tx2gene": "tx2gene.csv", "genes": "genes.csv"},
    }
    write_json(output_dir / "index.json", index_doc)
    return index_doc


def citations(index_doc: dict | None) -> list[dict]:
    """The works this run used: salmon; the decoy method with a genome; GENCODE for a GENCODE annotation."""
    out = [CITATIONS["salmon"]]
    if index_doc and index_doc["decoys"]:
        out.append(CITATIONS["decoys"])
    if index_doc and index_doc["source"].startswith("GENCODE"):
        out.append(CITATIONS["gencode"])
    return out


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir or not Path(input_dir).is_dir():
        log.error("INPUT_DIR (an existing folder) and OUTPUT_DIR must be set")
        return 2
    try:
        params = read_parameters()
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}
    status, error, index_doc = "ok", None, None
    try:
        index_doc = build(input_dir, output_dir, params, timings)
    except (InputError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        status, error = "failed", str(exc)
        log.error("%s", exc)
    peak = resources.peak_rss_mb()
    provenance = {
        "app": APP_NAME,
        "version": (REPO_ROOT / "VERSION").read_text().strip(),
        "tools": {"salmon": index_doc["salmon_version"] if index_doc else None, "python": sys.version.split()[0]},
        "parameters": params,
        "inputs": index_doc["inputs"] if index_doc else None,
        "started": started.isoformat(timespec="seconds"),
        "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "samples": [],
        "status": status,
        "error": error,
        "seconds": round(time.monotonic() - t_start, 1),
        "timings": timings,
        "threads": resources.threads(params["THREADS"]),
        "peak_rss_mb": max(peak.values()),
        "peak_rss_mb_by_process": peak,
        "cgroup_peak_mb": resources.cgroup_peak_mb(),
        "citations": citations(index_doc),
    }
    write_json(output_dir / "provenance.json", provenance)
    log.info("%s in %.1fs; peak RSS %s MB (salmon %s MB)", status, provenance["seconds"], peak["self"],
             peak["children"])
    return 0 if status == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
