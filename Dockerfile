# start with an image with conda installed
FROM condaforge/mambaforge AS compile-image

WORKDIR /data

COPY . ./fieldbioinformatics/

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt-get clean && apt-get autoclean

# copy in artic
RUN cd /data/fieldbioinformatics && \
  mamba install conda -n base -c conda-forge -c defaults -c bioconda && \
  mamba env create -f /data/fieldbioinformatics/environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "artic", "/bin/bash", "-c"]
RUN mamba install -c conda-forge -n artic python=3.9 conda-pack && \
  cd /data/fieldbioinformatics && \
  pip install .

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda list && \
  conda-pack -n artic -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar && /venv/bin/conda-unpack

SHELL ["/bin/bash", "-c"]

RUN conda clean --all &&\
  conda remove --name artic --all

# build artic
WORKDIR /data/fieldbioinformatics
RUN source /venv/bin/activate && pip install --user --no-cache-dir . \ 
  && pip uninstall -y tensorflow keras pyabpoa \
  && conda install -y -c conda-forge tensorflow~=2.11 keras~=2.11

# build image
FROM debian:bullseye-slim AS runtime-image

COPY --from=compile-image /root/.local /root/.local
ENV PATH=/root/.local/bin:$PATH

# Copy /venv from the previous stage:
COPY --from=compile-image /venv /venv

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt install -y procps && \
  apt-get clean && apt-get autoclean

WORKDIR /data

SHELL ["/bin/bash", "-c"]

# to allow streamed log output
ENV PYTHONUNBUFFERED=1
ENV PATH=/venv/bin:$PATH

CMD ["/bin/bash"]
