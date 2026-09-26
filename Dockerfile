FROM python:3.12-slim
WORKDIR /app
COPY pyproject.toml README.md LICENSE /app/
COPY src /app/src
RUN python -m pip install --no-cache-dir '.[de]' && useradd --uid 10001 --create-home analyst
USER analyst
WORKDIR /data
ENV MPLCONFIGDIR=/tmp/matplotlib
ENTRYPOINT ["sigtrellis"]

