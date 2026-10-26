FROM python:3.11-slim@sha256:b7a9e5f6a93d6d4a31f129bbd4bdee69efccb6ded0fb230d78aafaf37af98377

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR=/tmp/matplotlib
COPY requirements/runtime-linux-py311.lock /opt/runtime.lock
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip --isolated install --require-hashes --no-deps -r /opt/runtime.lock \
    && python -m pip check
WORKDIR /app
COPY src/oncorna/ /app/oncorna/
USER 65532:65532
CMD ["python", "-B", "-m", "oncorna.synthetic_demo"]
