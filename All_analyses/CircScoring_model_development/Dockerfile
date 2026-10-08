FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake libcurl4-openssl-dev libssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/CircScoring_model_development
COPY . .

RUN Rscript install_exact_packages.R

CMD ["Rscript", "CircScoring_model_development.R", \
     "--input", "data/P1N1_training_data.tsv", \
     "--output", "results", "--threads", "1"]

