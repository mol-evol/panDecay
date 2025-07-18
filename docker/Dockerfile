# panDecay Docker Container
# Production-ready container with PAUP*, MrBayes, and all dependencies
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    cmake \
    wget \
    curl \
    unzip \
    gzip \
    tar \
    git \
    libopenmpi-dev \
    openmpi-bin \
    openmpi-common \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libdb5.3-dev \
    libbz2-dev \
    libexpat1-dev \
    liblzma-dev \
    tk-dev \
    libffi-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN useradd -m -s /bin/bash pandecay && \
    mkdir -p /app /data /output && \
    chown -R pandecay:pandecay /app /data /output

# Install PAUP* (academic version)
RUN mkdir -p /opt/paup && \
    cd /opt/paup && \
    wget -q https://phylosolutions.com/paup-test/paup4a168_centos64.gz && \
    gunzip paup4a168_centos64.gz && \
    chmod +x paup4a168_centos64 && \
    ln -s /opt/paup/paup4a168_centos64 /usr/local/bin/paup && \
    rm -f *.gz

# Install MrBayes from source for better compatibility
RUN cd /tmp && \
    wget -q https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/MrBayes-3.2.7.tar.gz && \
    tar -xzf MrBayes-3.2.7.tar.gz && \
    cd MrBayes-3.2.7 && \
    ./configure --enable-mpi=yes --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd / && \
    rm -rf /tmp/MrBayes-3.2.7*

# Install Python dependencies
COPY ../requirements.txt /app/
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r requirements.txt

# Copy application files
COPY ../panDecay.py /app/
COPY ../src/ /app/src/
COPY ../examples/ /app/examples/

# Create config directory
RUN mkdir -p /app/config

# Create entrypoint script
COPY <<'EOF' /app/docker-entrypoint.sh
#!/bin/bash

set -e

# Function to display help
show_help() {
    echo "panDecay Docker Container"
    echo "========================="
    echo ""
    echo "Usage: docker run [OPTIONS] pandecay [COMMAND] [ARGS...]"
    echo ""
    echo "Commands:"
    echo "  analyze    Run phylogenetic decay analysis"
    echo "  config     Generate configuration templates"
    echo "  validate   Validate configuration files"
    echo "  help       Show this help message"
    echo "  bash       Open interactive shell"
    echo ""
    echo "Examples:"
    echo "  # Basic analysis with ML"
    echo "  docker run -v \$(pwd)/data:/data -v \$(pwd)/output:/output pandecay analyze /data/alignment.fas --analysis ml"
    echo ""
    echo "  # Analysis with YAML config"
    echo "  docker run -v \$(pwd):/workspace pandecay analyze --config /workspace/config.yaml"
    echo ""
    echo "  # Generate YAML config template"
    echo "  docker run -v \$(pwd):/workspace pandecay config --yaml /workspace/config.yaml"
    echo ""
    echo "  # Interactive mode"
    echo "  docker run -it -v \$(pwd):/workspace pandecay bash"
    echo ""
    echo "Environment Variables:"
    echo "  PANDECAY_THREADS    Number of threads (default: auto)"
    echo "  PANDECAY_DEBUG      Enable debug mode (default: false)"
    echo "  PANDECAY_FORMAT     Visualization format (default: both)"
    echo ""
}

# Set default environment variables
export PANDECAY_THREADS=${PANDECAY_THREADS:-auto}
export PANDECAY_DEBUG=${PANDECAY_DEBUG:-false}
export PANDECAY_FORMAT=${PANDECAY_FORMAT:-both}

# Handle different commands
case "${1:-help}" in
    analyze)
        shift
        echo "Running panDecay analysis..."
        python3 /app/panDecay.py "$@" \
            --threads "${PANDECAY_THREADS}" \
            --viz-format "${PANDECAY_FORMAT}" \
            ${PANDECAY_DEBUG:+--debug}
        ;;
    
    config)
        shift
        case "${1:-help}" in
            --yaml|yaml)
                shift
                output_file="${1:-/output/pandecay_config.yaml}"
                echo "📝 Generating YAML configuration template..."
                python3 /app/panDecay.py --generate-yaml-config "${output_file}"
                ;;
            --ini|ini)
                shift
                output_file="${1:-/output/pandecay_config.ini}"
                echo "📝 Generating INI configuration template..."
                python3 /app/panDecay.py --generate-config "${output_file}"
                ;;
            *)
                echo "Usage: config [yaml|ini] [output_file]"
                exit 1
                ;;
        esac
        ;;
    
    validate)
        shift
        config_file="${1:-/app/config/example_config.yaml}"
        echo "Validating configuration file: ${config_file}"
        python3 /app/config_loader.py validate "${config_file}"
        ;;
    
    migrate)
        shift
        ini_file="${1}"
        yaml_file="${2}"
        if [[ -z "${ini_file}" || -z "${yaml_file}" ]]; then
            echo "Usage: migrate <ini_file> <yaml_file>"
            exit 1
        fi
        echo "Migrating configuration from INI to YAML..."
        python3 /app/config_loader.py migrate "${ini_file}" "${yaml_file}"
        ;;
    
    test)
        echo "🧪 Running system tests..."
        echo "Testing PAUP*..."
        paup -v || echo "PAUP* test completed"
        echo "Testing MrBayes..."
        mb -v || echo "MrBayes test completed"
        echo "Testing Python modules..."
        python3 -c "
import panDecay
import config_loader
import dual_visualization
import async_constraint_processor
print('All modules imported successfully')
print(f'Async processing: {panDecay.HAS_ASYNC_PROCESSING}')
print(f'Enhanced config: {panDecay.HAS_ENHANCED_CONFIG}')
print(f'Dual visualization: {panDecay.HAS_DUAL_VISUALIZATION}')
"
        ;;
    
    bash|shell)
        echo "🐚 Opening interactive shell..."
        exec /bin/bash
        ;;
    
    help|--help|-h)
        show_help
        ;;
    
    *)
        # If first argument looks like a file or starts with -, pass to panDecay directly
        if [[ -f "${1}" || "${1}" == -* ]]; then
            echo "Running panDecay with arguments: $*"
            python3 /app/panDecay.py "$@"
        else
            echo "ERROR: Unknown command: ${1}"
            echo ""
            show_help
            exit 1
        fi
        ;;
esac
EOF

# Make entrypoint executable
RUN chmod +x /app/docker-entrypoint.sh

# Switch to non-root user
USER pandecay

# Create volume mount points
VOLUME ["/data", "/output", "/workspace"]

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python3 -c "import panDecay; print('OK')" || exit 1

# Set metadata
LABEL maintainer="panDecay Project" \
      version="1.1.0" \
      description="Phylogenetic decay analysis container with PAUP* and MrBayes" \
      org.opencontainers.image.title="panDecay" \
      org.opencontainers.image.description="Phylogenetic decay analysis with ML, Bayesian, and Parsimony methods" \
      org.opencontainers.image.version="1.1.0" \
      org.opencontainers.image.authors="panDecay Project" \
      org.opencontainers.image.url="https://github.com/your-repo/panDecay" \
      org.opencontainers.image.documentation="https://github.com/your-repo/panDecay/README.md"

# Expose port for potential web interface
EXPOSE 8080

# Set entrypoint
ENTRYPOINT ["/app/docker-entrypoint.sh"]
CMD ["help"]