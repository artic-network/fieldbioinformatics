FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

RUN sed -i 's/name: artic/name: base/' /tmp/env.yml && \
    micromamba install --yes --file /tmp/env.yml && \
    micromamba clean --all --yes && \
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
