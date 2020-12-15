ARG PY_VERSION=3.7
FROM continuumio/miniconda3:4.8.2-alpine AS builder

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
      maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /mbuild

WORKDIR /mbuild

RUN conda update conda -yq && \
	conda config --set always_yes yes --set changeps1 no && \
	conda config --add channels omnia && \
	conda config --add channels conda-forge && \
	. /opt/conda/etc/profile.d/conda.sh && \
	conda create -n mbuild-docker python=$PY_VERSION nomkl --file requirements-dev.txt && \
	conda activate mbuild-docker && \
        python setup.py install && \
	echo "source activate mbuild-docker" >> \
	/home/anaconda/.profile && \
	conda clean -afy && \
	mkdir /home/anaconda/data && \
	chown -R anaconda:anaconda /mbuild && \
	chown -R anaconda:anaconda /opt && \
	chown -R anaconda:anaconda /home/anaconda


WORKDIR /home/anaconda

COPY devtools/docker-entrypoint.sh /entrypoint.sh

RUN chmod a+x /entrypoint.sh

USER anaconda

ENTRYPOINT ["/entrypoint.sh"]
CMD ["jupyter"]
