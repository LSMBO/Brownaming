FROM condaforge/miniforge3:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV LOCAL_DB_PATH=/db

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml /app/environment.yml

# Create conda environment
RUN mamba env create -f environment.yml && \
    mamba clean -afy

# Activate environment by default
SHELL ["conda", "run", "-n", "brownaming", "/bin/bash", "-c"]

# Create necessary directories
RUN mkdir -p /app/runs

# Note: Python files, models, and time_prediction_model/ are mounted as volumes
# This allows modifications without rebuilding the image

# Set the entrypoint to use the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "brownaming"]

# Default command
CMD ["python", "/app/main.py", "--help"]
