# hash:sha256:8de4e5315b71b114cf1c9c761f2a43edbcb6cf57253126239d2dea1ffcc7d2f8
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/miniconda3:25.3.0-3-python3.12.10-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
