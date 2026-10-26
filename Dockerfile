# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:1.5.8

ARG ENVIRONMENT_FILE=environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER ${ENVIRONMENT_FILE} /tmp/deepcscn-environment.yml
RUN micromamba create --yes --file /tmp/deepcscn-environment.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/DeepCSCN
RUN micromamba run --name deepcscn R CMD INSTALL --no-multiarch /opt/DeepCSCN

ENV PATH=/opt/conda/envs/deepcscn/bin:${PATH}
ENV RETICULATE_PYTHON=/opt/conda/envs/deepcscn/bin/python
WORKDIR /work
CMD ["R"]
