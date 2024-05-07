ARG PY_VERSION=3.10
FROM continuumio/miniconda3:4.10.3-alpine AS builder

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
  maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /mbuild

WORKDIR /mbuild

# Create a group and user
RUN addgroup -S anaconda && adduser -S anaconda -G anaconda

RUN apk update && apk add libarchive &&\
  conda update conda -yq && \
  conda config --set always_yes yes --set changeps1 no && \
  . /opt/conda/etc/profile.d/conda.sh && \
  sed -i -E "s/python.*$/python="$(PY_VERSION)"/" environment-dev.yml && \
  conda install -c conda-forge mamba && \
  mamba env create --file environment-dev.yml && \
  conda activate mbuild-dev && \
  mamba install -c conda-forge jupyter python="$PY_VERSION" && \
  python setup.py install && \
  echo "source activate mbuild-dev" >> \
  /home/anaconda/.profile && \
  conda clean -afy && \
  mkdir -p /home/anaconda/data && \
  chown -R anaconda:anaconda /mbuild && \
  chown -R anaconda:anaconda /opt && \
  chown -R anaconda:anaconda /home/anaconda

WORKDIR /home/anaconda

COPY devtools/docker-entrypoint.sh /entrypoint.sh

RUN chmod a+x /entrypoint.sh

USER anaconda

ENTRYPOINT ["/entrypoint.sh"]
CMD ["jupyter"]
