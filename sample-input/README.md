# sample-input

rnaseq-index takes a transcriptome FASTA, a GTF of the same release and, optionally, the genome
FASTA of the same assembly. The demo reference of the bulk RNA-seq workflows is GENCODE mouse
release M39, downloaded by `scripts/fetch-sample-data.sh` in `workflow-bulk-rnaseq`:

```sh
cd ../workflow-bulk-rnaseq && scripts/fetch-sample-data.sh
ls data/reference/gencode-mouse-M39/
# gencode.vM39.transcripts.fa.gz  gencode.vM39.annotation.gtf.gz  GRCm39.primary_assembly.genome.fa.gz
```

Lay them out as the platform does when the three inputs are linked (each in its own folder),
or put them side by side (then the genome is told apart by its name: `genome`,
`primary_assembly`, `toplevel` or `.dna.`):

```sh
R=../workflow-bulk-rnaseq/data/reference/gencode-mouse-M39
mkdir -p in/transcriptome in/annotation in/genome
ln -s $PWD/$R/gencode.vM39.transcripts.fa.gz in/transcriptome/
ln -s $PWD/$R/gencode.vM39.annotation.gtf.gz in/annotation/
ln -s $PWD/$R/GRCm39.primary_assembly.genome.fa.gz in/genome/     # optional: decoys
mkdir -p tmp && docker run --rm --platform linux/amd64 --memory 6g --network none --read-only -v $PWD/tmp:/tmp \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v $PWD/in:/data/input:ro -v $PWD/out:/data/output rnaseq-index:dev
```

(With symlinks, mount the folder they point into as well, or copy the files.) The workflow
builders keep the transcripts-only index under
`workflow-bulk-rnaseq/data/work/gencode-mouse-M39/index/`; rnaseq-quant reads it as its `index`
input.

Ensembl works the same way: `Mus_musculus.GRCm39.cdna.all.fa.gz` (plus `ncrna` for non-coding
genes, concatenated into one file), `Mus_musculus.GRCm39.115.gtf.gz` and
`Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`; transcript versions are taken from the GTF's
`transcript_version` attribute.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python    # data/fixture/{transcriptome,genome,annotation}/
make run PYTHON=.venv/bin/python        # Docker over data/fixture/ -> data/output/
```
