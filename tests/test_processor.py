"""The synthetic genome through the processor (both layouts; needs salmon), and the annotation pieces."""
from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
from pathlib import Path

import pytest

from processor import annotation as A
from processor import main as app
from tests.make_fixture import make_fixture

needs_salmon = pytest.mark.skipif(not shutil.which(os.environ.get("SALMON_BIN", "salmon")),
                                  reason="salmon not on PATH (set SALMON_BIN)")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.DEFAULTS:
        monkeypatch.delenv(name, raising=False)
    for k, v in params.items():
        monkeypatch.setenv(k, v)
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@needs_salmon
@pytest.mark.parametrize("flat", [False, True])
def test_build_with_decoys(tmp_path, monkeypatch, flat):
    fx = make_fixture(tmp_path / "in", flat=flat)
    out = tmp_path / "out"
    assert run(fx, out, monkeypatch) == 0
    assert {p.name for p in out.iterdir()} == {"salmon", "tx2gene.csv", "genes.csv", "index.json", "provenance.json"}
    info = json.loads((out / "index.json").read_text())
    assert info["organism"] == "Mus musculus" and info["source"] == "GENCODE M39" and info["assembly"] == "GRCm39"
    assert info["decoys"] is True and info["n_decoy_sequences"] == 2 and info["gencode_headers"] is True
    assert info["n_transcripts"] == 7 and info["n_genes"] == 6 and info["n_transcripts_not_in_annotation"] == 0
    assert set(info["inputs"]) == {"transcriptome", "genome", "annotation"}
    assert all(len(v["sha256"]) == 64 for v in info["inputs"].values())
    tx = rows(out / "tx2gene.csv")
    assert tx[0] == {"transcript_id": "ENSMUST00000000001.1", "gene_id": "ENSMUSG00000000001.1",
                     "gene_name": "Alpha", "gene_type": "protein_coding"}
    genes = {r["gene_id"]: r for r in rows(out / "genes.csv")}
    assert genes["ENSMUSG00000000001.1"]["length"] == "900"          # 300 shared + 300 + 300
    assert genes["ENSMUSG00000000006.1"] == {"gene_id": "ENSMUSG00000000006.1", "gene_name": "Zeta",
                                             "gene_type": "snRNA", "chrom": "chr2", "length": "150"}
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["app"] == "rnaseq-index" and prov["version"] == "0.1.0" and prov["status"] == "ok"
    assert [c["doi"] for c in prov["citations"]] == ["10.1038/nmeth.4197", "10.1186/s13059-020-02151-8",
                                                     "10.1093/nar/gkae1078"]
    assert all(c["pmid"] and c["text"] for c in prov["citations"]) and prov["peak_rss_mb"] > 0


@needs_salmon
def test_build_without_genome(tmp_path, monkeypatch):
    fx = make_fixture(tmp_path / "in")
    shutil.rmtree(fx / "genome")
    assert run(fx, tmp_path / "out", monkeypatch, ORGANISM="mouse (test)", SOURCE="custom") == 0
    info = json.loads((tmp_path / "out" / "index.json").read_text())
    assert info["decoys"] is False and info["organism"] == "mouse (test)" and info["source"] == "custom"
    assert info["detected_source"] == "GENCODE M39" and info["index_version"] and info["sshash_format_version"]
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    # no decoys: no decoy paper; the GENCODE citation follows what was detected, not SOURCE
    assert [c["doi"] for c in prov["citations"]] == ["10.1038/nmeth.4197", "10.1093/nar/gkae1078"]


def test_wrong_annotation_fails(tmp_path, monkeypatch):
    fx = make_fixture(tmp_path / "in")
    with gzip.open(fx / "annotation" / "gencode.vM39.annotation.gtf.gz", "wt") as fh:
        fh.write('chr1\tX\ttranscript\t1\t100\t.\t+\t.\tgene_id "OTHER1"; transcript_id "OTHERT1";\n')
    assert run(fx, tmp_path / "out", monkeypatch) == 1
    assert json.loads((tmp_path / "out" / "provenance.json").read_text())["status"] == "failed"


def test_missing_inputs_and_bad_parameters(tmp_path, monkeypatch):
    (tmp_path / "empty").mkdir()
    assert run(tmp_path / "empty", tmp_path / "out", monkeypatch) == 1
    for params in ({"KMER_LEN": "32"}, {"KMER_LEN": "x"}, {"THREADS": "-1"}, {"SOURCE": "a\nb"}):
        assert run(tmp_path / "empty", tmp_path / "out", monkeypatch, **params) == 2


def test_root_layout_tells_genome_from_transcriptome(tmp_path):
    fx = make_fixture(tmp_path, flat=True)
    found = app.find_inputs(fx)
    assert found["transcriptome"].name == "gencode.vM39.transcripts.fa.gz"
    assert found["genome"].name == "GRCm39.primary_assembly.genome.fa.gz"
    (fx / "other.fa").write_text(">x\nACGT\n")
    with pytest.raises(app.InputError):
        app.find_inputs(fx)


def test_ensembl_gtf_versions_and_names(tmp_path):
    gtf = tmp_path / "Rattus_norvegicus.mRatBN7.2.115.gtf"
    gtf.write_text(
        "#!genome-build mRatBN7.2\n"
        '1\tensembl\ttranscript\t10\t60\t.\t+\t.\tgene_id "ENSRNOG1"; gene_version "3"; transcript_id "ENSRNOT1"; '
        'transcript_version "2"; gene_name "Abc1"; gene_biotype "protein_coding";\n'
        '1\tensembl\texon\t10\t30\t.\t+\t.\tgene_id "ENSRNOG1"; gene_version "3"; transcript_id "ENSRNOT1"; '
        'transcript_version "2";\n'
        '1\tensembl\texon\t21\t60\t.\t+\t.\tgene_id "ENSRNOG1"; gene_version "3"; transcript_id "ENSRNOT1"; '
        'transcript_version "2";\n')
    ann = A.read_gtf(gtf)
    tx, missing = A.tx2gene_rows(["ENSRNOT1.2", "ENSRNOT9.1"], ann)
    assert tx[0] == {"transcript_id": "ENSRNOT1.2", "gene_id": "ENSRNOG1.3", "gene_name": "Abc1",
                     "gene_type": "protein_coding"}
    assert missing == ["ENSRNOT9.1"] and tx[1]["gene_id"] == "ENSRNOT9.1"
    assert A.gene_rows(tx, ann)[0]["length"] == 51                      # overlapping exons counted once
    assert A.describe(ann, gtf.name) == ("Rattus norvegicus", "Ensembl 115", "mRatBN7.2")
    # unversioned FASTA ids still match
    assert A.tx2gene_rows(["ENSRNOT1"], ann)[0][0]["gene_id"] == "ENSRNOG1.3"


def test_truncated_gzip_fails_with_provenance_and_no_partial_outputs(tmp_path, monkeypatch):
    fx = make_fixture(tmp_path / "in")
    gtf = fx / "annotation" / "gencode.vM39.annotation.gtf.gz"
    data = gtf.read_bytes()
    gtf.write_bytes(data[: len(data) // 2])
    out = tmp_path / "out"
    assert run(fx, out, monkeypatch) == 1
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["status"] == "failed" and "EOFError" in prov["error"]
    assert [p.name for p in out.iterdir()] == ["provenance.json"]


def test_symlinks_out_of_the_input_are_not_read(tmp_path, caplog):
    fx = make_fixture(tmp_path / "in")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "secret.gtf").write_text("x")
    # links the platform places at the root or in a port's folder are trusted; deeper ones must stay inside
    (fx / "annotation" / "deeper").mkdir()
    (fx / "annotation" / "deeper" / "escape.gtf").symlink_to(outside / "secret.gtf")
    (fx / "annotation" / "deeper" / "escape-dir").symlink_to(outside, target_is_directory=True)
    found = app.find_inputs(fx)
    assert found["annotation"].name == "gencode.vM39.annotation.gtf.gz"
    assert "annotation/deeper/escape.gtf: resolves outside the input" in caplog.text
    assert "annotation/deeper/escape-dir: resolves outside the input" in caplog.text
    # port folders that are themselves links (as the platform links inputs) are followed
    linked = tmp_path / "linked"
    linked.mkdir()
    for port in ("transcriptome", "annotation"):
        (linked / port).symlink_to(fx / port, target_is_directory=True)
    found = app.find_inputs(linked)
    assert found["transcriptome"] is not None and found["annotation"] is not None


def test_threads_without_a_quota_use_the_available_cpus(monkeypatch):
    from processor import safe_inputs

    monkeypatch.setattr(safe_inputs, "_cgroup_cpus", lambda: None)
    assert safe_inputs.cpu_threads(3) == 3
    assert safe_inputs.cpu_threads(0) == (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
                                          else os.cpu_count())
