FROM continuumio/miniconda3:4.10.3
LABEL authors="Matt Jaquiery" \
      description="Docker image for development"

RUN mkdir /usr/share/man/man1/  # for Java -- https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=863199#23
RUN apt-get update; apt-get install -y  gcc procps jq default-jre

COPY Docker/conda_environments/kat.yml /
RUN conda env create -f /kat.yml && conda clean -a

COPY Docker/conda_environments/confindr.yml /
RUN conda env create -f /confindr.yml && conda clean -a

COPY Docker/conda_environments/quast.yml /
RUN conda env create -f /quast.yml && conda clean -a

COPY Docker/conda_environments/assembly.yml /
RUN conda env create -f /assembly.yml && conda clean -a

COPY Docker/conda_environments/trimmomatic.yml /
RUN conda env create -f /trimmomatic.yml && conda clean -a

COPY Docker/requirements.txt /
RUN pip install -r /requirements.txt

COPY Docker/confindr_database.tar.gz /confindr_database.tar.gz
RUN tar xfz /confindr_database.tar.gz && rm /confindr_database.tar.gz && chown -R root:root /confindr_database && chmod -R o+w /confindr_database

COPY Docker/scripts/filter_contigs.py /usr/local/bin/

RUN wget -qO- get.nextflow.io | bash
RUN mv nextflow /usr/local/bin/

ENV PATH /opt/conda/bin:/opt/conda/condabin:/opt/conda/envs/trimmomatic/bin:/opt/conda/envs/assembly/bin:/opt/conda/envs/quast/bin:/opt/conda/envs/confindr/bin:/opt/conda/envs/kat/bin:$PATH
