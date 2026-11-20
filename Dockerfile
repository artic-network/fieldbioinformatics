FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Clair3 >= 2.0.5 bundles models into bin/models (artic's default model dir); drop them,
# plus the hardlinked pkgs/ copy, so the default image has none and INCLUDE_MODELS gets artic's manifest only.
RUN sed -i 's/name: artic/name: base/' /tmp/env.yml && \
    micromamba install --yes --file /tmp/env.yml && \
    micromamba clean --all --yes && \
    rm -rf "$MAMBA_ROOT_PREFIX/bin/models" "$MAMBA_ROOT_PREFIX"/pkgs/clair3-*/bin/models && \
    rm /tmp/env.yml

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY . ./fieldbioinformatics/

USER root

RUN python3 -m pip install --no-cache-dir ./fieldbioinformatics && \
    rm -rf ./fieldbioinformatics

ARG INCLUDE_MODELS=false
RUN if [ "$INCLUDE_MODELS" = "true" ]; then artic_get_models; fi

USER $MAMBA_USER

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

CMD ["/bin/bash"]
