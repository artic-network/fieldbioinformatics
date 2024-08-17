FROM mambaorg/micromamba:1.5.8

COPY . ./fieldbioinformatics/

USER root

RUN apt-get update && apt-get install -y --no-install-recommends build-essential wget

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

RUN sed -i 's/name: artic/name: base/' /tmp/env.yml

RUN micromamba install --yes --file /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN python3 -m pip install ./fieldbioinformatics

USER $MAMBA_USER

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

CMD ["/bin/bash"]