# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components; their licences apply to the image. Nothing in it is GPL.

| Component | What | Licence | Source |
|---|---|---|---|
| salmon 2.8.0 | `salmon index`: the selective-alignment index (piscem-rs / sshash-rs based). The release's prebuilt `salmon-cli-x86_64-unknown-linux-gnu.tar.xz`, checked against the sha256 published with the release (`10a7f799…1efa6`). Licence file at `/usr/local/share/licenses/salmon/LICENSE`. salmon 2.x is a Rust rewrite; the binary statically includes Rust crates under MIT / Apache-2.0 / BSD licences (listed in salmon's `Cargo.lock` at tag `v2.8.0`). | BSD-3-Clause (salmon ≥ 2.0; C++ salmon ≤ 1.12 was GPL-3.0) | https://github.com/COMBINE-lab/salmon/releases/tag/v2.8.0 |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12), pinned by digest. The app uses only the Python standard library. | PSF-2.0 (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

No reference data is bundled: the transcriptome, genome and annotation are inputs.

## Data used in testing (not in the image)

| Data | Licence | Source |
|---|---|---|
| GENCODE mouse release M39 (GRCm39): `gencode.vM39.transcripts.fa.gz`, `gencode.vM39.annotation.gtf.gz`, `GRCm39.primary_assembly.genome.fa.gz` | Public data; EMBL-EBI terms of use (https://www.ebi.ac.uk/about/terms-of-use); cite GENCODE 2025 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/ |
