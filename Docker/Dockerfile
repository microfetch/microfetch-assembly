FROM continuumio/miniconda3:4.9.2
LABEL authors="Anthony Underwood" \
      description="Docker image for de novo assembly"

RUN apt update; apt install -y  gcc procps jq
COPY conda_environments/kat.yml /
RUN conda env create -f /kat.yml && conda clean -a

COPY conda_environments/assembly.yml /
RUN conda env create -f /assembly.yml && conda clean -a

ENV PATH  /opt/conda/envs/assembly/bin:/opt/conda/envs/kat/bin:$PATH
