# rnaseq-index: salmon (official release binary) + the Python standard library.
# linux/amd64; no network at run time; runs as uid 1000.
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS fetch

# salmon 2.8.0 (BSD-3-Clause), the release's x86_64 Linux build, checked against the sha256 the
# release publishes (salmon-cli-x86_64-unknown-linux-gnu.tar.xz.sha256).
ARG SALMON_VERSION=2.8.0
ARG SALMON_SHA256=10a7f7994df7ef0592ac191cc79a1abd41316463000b307647a9dcab21e1efa6
ADD --checksum=sha256:${SALMON_SHA256} \
    https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-cli-x86_64-unknown-linux-gnu.tar.xz /tmp/salmon.tar.xz
RUN apt-get update && apt-get install -y --no-install-recommends xz-utils && rm -rf /var/lib/apt/lists/* \
    && tar -xJf /tmp/salmon.tar.xz -C /tmp \
    && install -m 0755 /tmp/salmon-cli-x86_64-unknown-linux-gnu/salmon /usr/local/bin/salmon \
    && install -D -m 0644 /tmp/salmon-cli-x86_64-unknown-linux-gnu/LICENSE /usr/local/share/licenses/salmon/LICENSE

FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

COPY --from=fetch /usr/local/bin/salmon /usr/local/bin/salmon
COPY --from=fetch /usr/local/share/licenses/salmon /usr/local/share/licenses/salmon
RUN salmon --version

WORKDIR /processor
COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH=/processor \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    HOME=/tmp

# The platform runs the task as uid 1000; nothing is written outside OUTPUT_DIR and /tmp.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.main"]
