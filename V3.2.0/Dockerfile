# GRACE-FO PPP-AR POD V3.2.0
# ========================
# Build:  docker build -t gnss-pod:v3.2.0 .
# Run:
#   docker run --rm \
#     -v /path/to/gnss_pod/data:/data:ro \
#     -v /path/to/output:/data/results:rw \
#     -e OREKIT_DATA_PATH=/data/orekit \
#     gnss-pod:v3.2.0 \
#     --config /app/config.yaml
#
# Without Docker (direct run):
#   pip install -r requirements.txt
#   export OREKIT_DATA_PATH=/path/to/data/orekit
#   python eval_5day_orekit.py --config config.yaml

FROM eclipse-temurin:17-jre-jammy AS java-base
FROM python:3.12-slim

COPY --from=java-base /opt/java/openjdk /opt/java/openjdk
ENV JAVA_HOME=/opt/java/openjdk
ENV PATH="${JAVA_HOME}/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/
COPY eval_5day_orekit.py run_sequential_pod.py config.yaml ./
COPY references/ ./references/
COPY VERSION.md ./

RUN mkdir -p /data/results
ENV OREKIT_DATA_PATH=/data/orekit PYTHONUNBUFFERED=1

ENTRYPOINT ["python", "-u", "eval_5day_orekit.py"]
CMD ["--config", "/app/config.yaml"]
