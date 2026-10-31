# rnaseq-index

Builds the reference that bulk RNA-seq quantification needs, **once per genome and annotation
release**, so it can be kept as a workspace reference and linked into every
`rnaseq-quant` run (workflow A of the bulk RNA-seq workflows:
`rnaseq-quant (+ index from rnaseq-index) → bulk-de → bulk-report`).

**For a biologist.** RNA-seq reads are counted per gene by matching them to the known
transcript sequences of the organism. This app prepares that lookup once: it indexes the
transcript sequences with salmon and records which transcript belongs to which gene, what the
gene is called and how long it is. Given the genome as well, it adds the whole genome as
"decoys", so reads from introns or unannotated regions are recognised as such instead of being
forced onto a similar transcript (this lowers false expression of e.g. pseudogenes and
paralogues, at the cost of a much bigger build).

**For a computational analyst.** `salmon index -k 31` (salmon 2.8.0, selective alignment;
`--gencode` for GENCODE headers) over the transcript FASTA, with the genome appended and its
sequence names passed as `-d decoys.txt` when a genome is given (Srivastava et al. 2020).
`tx2gene.csv` is keyed by the FASTA's transcript ids (versions kept as in the headers;
Ensembl GTFs' `transcript_version` / `gene_version` attributes are appended to match), so
rnaseq-quant can aggregate without `ignoreTxVersion`. Gene length in `genes.csv` is the
exonic length (union of all exons of the gene, overlaps counted once).

## Ports

| Port | Direction | What |
|---|---|---|
| `transcriptome` | input (subfolder), required | One transcript FASTA (`.fa`, `.fasta`, `.fna`, optionally `.gz`): GENCODE `gencode.v*.transcripts.fa.gz` or Ensembl `*.cdna.all.fa.gz` (+ `ncrna` in the same file for non-coding genes). |
| `genome` | input (subfolder), optional | The genome FASTA of the same assembly (`*.primary_assembly.genome.fa.gz`), used as decoys. |
| `annotation` | input (subfolder), required | The GTF of the same release (`gencode.v*.annotation.gtf.gz`, `*.gtf.gz`). |
| `index` | output | `salmon/`, `tx2gene.csv`, `genes.csv`, `index.json`, `provenance.json` (the whole `OUTPUT_DIR`). |

Without links, everything may sit side by side in `INPUT_DIR`: the one `.gtf(.gz)` is the
annotation; one FASTA is the transcriptome; with two, the one named like a genome (`genome`,
`primary_assembly`, `toplevel`, `.dna.`) is the genome. Anything ambiguous stops the run
(exit 1) with a message.

## Parameters

| Name | Default | Effect |
|---|---|---|
| `KMER_LEN` | `31` | salmon's k (odd, 15–31). 31 for reads ≥ 75 bp; 23–25 for reads under ~50 bp. |
| `ORGANISM` | empty | Written to `index.json`. Empty: from the GENCODE GTF header (human, mouse) or an Ensembl GTF file name (`Genus_species.Assembly.Release.gtf.gz`). |
| `SOURCE` | empty | Written to `index.json` (e.g. `GENCODE M39`). Empty: from the GTF header / file name, else the GTF file name. |
| `THREADS` | `0` | Threads for salmon; 0 = the container's CPU quota. |

Bad values exit 2 before anything is read.

## Outputs

- `salmon/`: the salmon 2.x index. rnaseq-quant refuses an index from another salmon major
  version (salmon 2 cannot read C++ salmon ≤ 1.x indices and vice versa).
- `tx2gene.csv`: `transcript_id, gene_id, gene_name, gene_type`, one row per FASTA transcript.
  Transcripts the GTF lacks are their own gene (counted in `index.json`; more than half
  missing stops the run: annotation and transcriptome are from different releases).
- `genes.csv`: `gene_id, gene_name, gene_type, chrom, length` (exonic length in bp).
- `index.json`: `organism`, `source`, `assembly`, `salmon_version`, `kmer_len`, `decoys`
  (true/false), `n_decoy_sequences`, `gencode_headers`, `n_transcripts`, `n_genes`,
  `n_transcripts_not_in_annotation`, `inputs` (`file`, `sha256` of each input), `files`.
- `provenance.json`: app version, salmon version, parameters, inputs, timings, peak memory,
  `citations` (salmon; the decoy paper when a genome was given; GENCODE for a GENCODE GTF).

## Runtime

`computeTypes: [standard]`, 8 vCPU / 30 GiB, 4 h timeout. salmon's external k-mer sort gets
half of the container's memory (`--ramLimit`); the build's scratch files go to `/tmp`
(ephemeral disk: ~8 GB for mouse transcripts + genome).

Measured with `docker run --platform linux/amd64 --cpus 4 --memory 6g` (amd64 emulation on an
M-series Mac; native amd64 is faster), GENCODE mouse M39 (481,956 transcripts, 78,374 genes),
timings from `provenance.json`:

| Build | wall time | peak memory | index size |
|---|---|---|---|
| Transcripts only (`--ramLimit 3`) | **174 s** (GTF + FASTA tables 41 s, salmon 130 s) | salmon 3.8 GB RSS, Python 0.8 GB; cgroup 5.3 GB incl. page cache | 1.3 GB (`salmon/`), 179 M distinct k-mers, 1.98 M unitigs |
| Transcripts + GRCm39 genome as decoys (61 sequences) | **not measured to completion**: killed by the 7 GB limit after 15 min, in sshash's final dictionary step | > 6.3 GB salmon RSS at `--ramLimit 3` | 2.21 billion distinct k-mers (12× the transcriptome's) |

The decoy build could not be measured within the 8 GB local Docker. salmon's documentation
says only that whole-genome decoys "require considerably more memory to build the index"
(salmon docs, mapping-based mode) and gives no number for salmon 2. The 30 GiB / 8 vCPU
default is therefore a conservative estimate (12× the k-mers; failed above 6.3 GB), not a
measurement: run the
decoy build once on `standard` with that size and check `peak_rss_mb` in `provenance.json`.
The transcripts-only index fits in 8 GiB.

Image: 213 MB (salmon 22 MB + python:3.12-slim).

## Limitations

- An index is only as good as its annotation: genes or isoforms missing from the GTF / FASTA
  cannot be quantified, and their reads end up on similar transcripts (or as unmapped).
- Transcriptome and GTF must be the same release; genome and transcriptome the same assembly.
  The app checks ids (transcriptome vs GTF) and name clashes (genome vs transcripts), not
  sequences.
- salmon collapses transcripts with identical sequences (listed in
  `salmon/duplicate_clusters.tsv`); their reads go to the retained copy, which can sit in a
  different gene for a few readthrough / paralogue cases.
- Organism and source are inferred only from GENCODE headers and Ensembl file names; set
  `ORGANISM` / `SOURCE` for anything else.

## Licences

Our code: Apache-2.0 (`LICENSE`). salmon 2.x is BSD-3-Clause (C++ salmon ≤ 1.12 was
GPL-3.0); nothing in the image is GPL. See `THIRD_PARTY.md`.

## Credits

- **salmon** (Rob Patro and the COMBINE-lab; Patro, Duggal, Love, Irizarry, Kingsford) builds
  the index; salmon 2.x is the Rust rewrite.
- **Selective alignment and decoys** (Srivastava, Malik, Sarkar, Zakeri, Almodaresi, Soneson,
  Love, Kingsford, Patro) are the method behind the decoy-aware index.
- **GENCODE** (the GENCODE consortium; Mudge, Frankish et al.) provides the demo reference,
  mouse release M39 (GRCm39).

## Citations

Checked against Crossref (title, authors and their order, journal, year, volume, issue,
pages) and PubMed (E-utilities `esearch` on the DOI, then `esummary`: title, journal, date,
volume, pages, author count, first and last author). Full author lists are in
`CITATION.cff`.

- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware
  quantification of transcript expression. *Nat Methods* 14(4):417-419 (2017).
  doi:10.1038/nmeth.4197. PMID 28263959.
- Srivastava A, Malik L, Sarkar H, Zakeri M, Almodaresi F, Soneson C, Love MI, Kingsford C,
  Patro R. Alignment and mapping methodology influence transcript abundance estimation.
  *Genome Biol* 21(1):239 (2020). doi:10.1186/s13059-020-02151-8. PMID 32894187. (decoys)
- Mudge JM, Carbonell-Sala S, Diekhans M, et al. (56 authors). GENCODE 2025: reference gene
  annotation for human and mouse. *Nucleic Acids Res* 53(D1):D966-D975 (2025).
  doi:10.1093/nar/gkae1078. PMID 39565199. (GENCODE's current reference publication, per
  gencodegenes.org)
