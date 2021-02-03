FROM continuumio/miniconda3:4.9.2
LABEL authors="Anthony Underwood" \
      description="Docker image for de novo assembly"

RUN apt update; apt install -y  gcc procps jq

COPY conda_environments/kat.yml /
RUN conda env create -f /kat.yml && conda clean -a

COPY conda_environments/confindr.yml /
RUN conda env create -f /confindr.yml && conda clean -a

COPY conda_environments/quast.yml /
RUN conda env create -f /quast.yml && conda clean -a

COPY conda_environments/assembly.yml /
RUN conda env create -f /assembly.yml && conda clean -a

COPY requirements.txt /
RUN pip install -r /requirements.txt

COPY confindr_database.tar.gz /confindr_database.tar.gz
RUN tar xfz /confindr_database.tar.gz && rm /confindr_database.tar.gz && chown -R root:root /confindr_database && chmod -R o+w /confindr_database

COPY scripts/filter_contigs.py /usr/local/bin/

ENV PATH /opt/conda/bin:/opt/conda/condabin:/opt/conda/envs/assembly/bin:/opt/conda/envs/quast/bin:/opt/conda/envs/confindr/bin:/opt/conda/envs/kat/bin:$PATH
