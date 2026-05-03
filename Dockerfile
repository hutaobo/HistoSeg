FROM python:3.12-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

FROM base AS core

RUN python -m pip install -U pip \
    && python -m pip install -e ".[threed]"

ENTRYPOINT ["histoseg-3d"]
CMD ["--help"]

FROM core AS viz

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libegl1 \
        libgl1 \
        libglib2.0-0 \
        libglx-mesa0 \
        xvfb \
    && rm -rf /var/lib/apt/lists/* \
    && python -m pip install -e ".[viz,docs]"

ENV PYVISTA_OFF_SCREEN=true

ENTRYPOINT ["histoseg-3d"]
CMD ["--help"]
