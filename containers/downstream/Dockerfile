# Build with scripts/build_downstream_container.py and the accepted Linux locks.
# BASE_IMAGE is resolved to an immutable digest by the build helper.
ARG BASE_IMAGE
FROM ${BASE_IMAGE} AS native
USER root
COPY locks/linux-64.explicit.txt /opt/katana-build/locks/linux-64.explicit.txt
COPY conda-cache/ /opt/conda/pkgs/
RUN micromamba --no-rc install --yes --name base --file /opt/katana-build/locks/linux-64.explicit.txt \
    && micromamba clean --all --yes
FROM native AS builder
COPY locks/pip-linux-64.lock /opt/katana-build/locks/pip-linux-64.lock
ENV PATH=/opt/conda/bin:$PATH \
    KATANA_CONTAINER=1 \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    HOME=/tmp
RUN python -m pip install --no-cache-dir --no-build-isolation --require-hashes --no-deps -r /opt/katana-build/locks/pip-linux-64.lock \
    && python -m pip check
COPY wheel/ /opt/katana-build/wheel/
RUN python -m pip install --no-cache-dir --no-deps /opt/katana-build/wheel/*.whl \
    && python -m pip check
FROM ${BASE_IMAGE}
USER root
COPY --from=builder /opt/conda/ /opt/conda/
COPY system-packages/ /opt/katana-build/system-packages/
COPY locks/system-packages-linux-amd64.json /opt/katana/system-packages-linux-amd64.json
# Exact official Debian packages are hash-verified by the context helper.
# No apt update, upgrade, solver or unpinned system dependency is run here.
RUN dpkg --install /opt/katana-build/system-packages/*.deb \
    && ps -p 1 -o pid=,comm= \
    && rm -rf /opt/katana-build/system-packages
ENV PATH=/opt/conda/bin:$PATH \
    KATANA_CONTAINER=1 \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    HOME=/tmp
COPY check_downstream_environment.py /opt/katana/check_downstream_environment.py
COPY build-inputs.json /opt/katana/build-inputs.json
# Nextflow supplies /bin/bash and the task command. No host activation needed.
ENTRYPOINT []
CMD ["katana-flow", "--help"]
