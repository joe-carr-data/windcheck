# Lockfile-frozen build recipe for the windcheck transaction and census
# tools.
#
#   docker build -t windcheck .
#   docker run --rm -v /data:/data windcheck \
#       transaction /data/candidate.tifxyz --out /data/final.tifxyz
#
# The transform stage runs the frozen operator scripts under bench/,
# which every published certificate cites by path; this image carries
# the whole source tree so the operator is exactly the published one.
# Python dependencies install with --frozen from the committed uv.lock.

FROM python:3.13-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ git \
    && rm -rf /var/lib/apt/lists/*
COPY --from=ghcr.io/astral-sh/uv:0.6.11 /uv /usr/local/bin/uv

WORKDIR /opt/windcheck
COPY . .
RUN g++ -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp \
    && uv sync --frozen --quiet

ENTRYPOINT ["uv", "run", "windcheck"]
CMD ["--help"]
