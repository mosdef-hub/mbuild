ARG PY_VERSION=3.7
FROM continuumio/miniconda3

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
      maintainer.url="https://mosdef.org"

SHELL ["/bin/bash", "-c"]

ADD . /mbuild

WORKDIR /mbuild

RUN conda update conda -yq

RUN conda config --set always_yes yes --set changeps1 no

RUN conda config --add channels omnia

RUN conda config --add channels janschulz

RUN conda config --add channels conda-forge

RUN conda config --add channels mosdef

RUN conda create -n mbuild-docker python=$PY_VERSION --file requirements-dev.txt

RUN source activate mbuild-docker && python setup.py install

RUN echo "source activate mbuild-docker" >> ~/.bashrc

ENV PATH "/opt/conda/envs/mbuild-docker/bin:$PATH"

WORKDIR /scratch