# panDecay Docker Guide

Complete guide for using panDecay with Docker containers, including installation, configuration, and advanced deployment scenarios.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Advanced Configuration](#advanced-configuration)
- [Docker Compose](#docker-compose)
- [Production Deployment](#production-deployment)
- [Development Workflow](#development-workflow)
- [Troubleshooting](#troubleshooting)

## Quick Start

### 1. Build Container

```bash
# Clone repository
git clone https://github.com/your-repo/panDecay.git
cd panDecay

# Build Docker image
./docker-deploy.sh build
```

### 2. Run Analysis

```bash
# Create data and output directories
mkdir -p data output

# Copy your alignment file to data/
cp alignment.fas data/

# Run analysis
./docker-deploy.sh run alignment.fas --analysis ml --model GTR
```

### 3. View Results

```bash
# Results are in the output directory
ls output/
```

That's it! Your phylogenetic decay analysis is complete.

## Installation

### Prerequisites

- Docker Engine 20.10+
- Docker Compose 2.0+ (optional, for multi-service deployment)
- At least 4GB available RAM
- 10GB available disk space

### Install Docker

#### Linux (Ubuntu/Debian)
```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Install Docker Compose
sudo apt-get install docker-compose-plugin
```

#### macOS
```bash
# Install Docker Desktop from https://docker.com/products/docker-desktop
# Or using Homebrew
brew install --cask docker
```

#### Windows
Download Docker Desktop from https://docker.com/products/docker-desktop

### Build panDecay Container

```bash
# Method 1: Using deployment script (recommended)
./docker-deploy.sh build

# Method 2: Direct Docker build
docker build -t pandecay:latest .

# Method 3: Using Docker Compose
docker-compose build
```

## Basic Usage

### Command Line Interface

The Docker container provides the same command-line interface as the standalone version:

```bash
# Run analysis with Docker
docker run --rm \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  pandecay analyze /data/alignment.fas --analysis ml

# Using deployment script (easier)
./docker-deploy.sh run alignment.fas --analysis ml
```

### Container Commands

```bash
# Show help
docker run pandecay help

# Run analysis
docker run pandecay analyze [OPTIONS] alignment_file

# Generate configuration
docker run pandecay config yaml output_file

# Validate configuration
docker run pandecay validate config_file

# Open interactive shell
docker run -it pandecay bash

# Run system tests
docker run pandecay test
```

### Volume Mounting

Map your local directories to container paths:

```bash
docker run --rm \
  -v /path/to/your/data:/data:ro \           # Input data (read-only)
  -v /path/to/output:/output \               # Output directory
  -v /path/to/config:/workspace/config:ro \  # Configuration files
  pandecay analyze /data/alignment.fas
```

### Environment Variables

Control container behavior with environment variables:

```bash
docker run --rm \
  -e PANDECAY_THREADS=8 \
  -e PANDECAY_DEBUG=true \
  -e PANDECAY_FORMAT=both \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  pandecay analyze /data/alignment.fas
```

## Advanced Configuration

### Environment Variables Reference

| Variable | Default | Description |
|----------|---------|-------------|
| `PANDECAY_THREADS` | `auto` | Number of analysis threads |
| `PANDECAY_DEBUG` | `false` | Enable debug mode |
| `PANDECAY_FORMAT` | `both` | Visualization format (static/interactive/both) |
| `PANDECAY_ASYNC` | `true` | Enable async constraint processing |

### Resource Limits

```bash
# Limit CPU and memory usage
docker run --rm \
  --cpus="4.0" \
  --memory="8g" \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  pandecay analyze /data/alignment.fas
```

### Custom Configuration

```bash
# Create configuration file
cat > config.yaml << EOF
analysis:
  type: ml+bayesian
  model: GTR
  gamma: true

computational:
  threads: 8
  async_constraints: true
  max_async_workers: 4

visualization:
  format: both
  theme: publication
EOF

# Run with custom configuration
docker run --rm \
  -v $(pwd):/workspace \
  pandecay analyze --config /workspace/config.yaml
```

## Docker Compose

Docker Compose provides easier management of multi-service deployments.

### Configuration File

The `docker-compose.yml` includes three service profiles:

- **production**: Standard single-analysis service
- **development**: Interactive development environment
- **batch**: Batch processing for multiple analyses

### Environment Setup

```bash
# Copy and customize environment
cp .env.example .env

# Edit environment variables
nano .env
```

### Service Profiles

#### Production Service

```bash
# Start production service
docker-compose up -d pandecay

# Run analysis
docker-compose exec pandecay analyze /data/alignment.fas --analysis ml

# View logs
docker-compose logs -f pandecay

# Stop service
docker-compose down
```

#### Development Service

```bash
# Start development environment
COMPOSE_PROFILES=development docker-compose up -d pandecay-dev

# Access interactive shell
docker-compose exec pandecay-dev bash

# Stop development environment
docker-compose down
```

#### Batch Processing

```bash
# Prepare batch data
mkdir -p batch_data batch_output batch_config

# Place multiple alignment files in batch_data/
# Place configuration files in batch_config/

# Run batch processing
COMPOSE_PROFILES=batch docker-compose up pandecay-batch
```

### Volume Management

```yaml
# In docker-compose.yml
volumes:
  - ./data:/data:ro                    # Input data (read-only)
  - ./output:/output                   # Output directory
  - ./config:/workspace/config:ro      # Configuration files
  - pandecay-cache:/home/pandecay/.cache  # Persistent cache
```

## Production Deployment

### Resource Planning

#### Minimum Requirements
- **CPU**: 2 cores
- **RAM**: 4GB
- **Storage**: 10GB
- **Network**: Standard connectivity

#### Recommended for Large Datasets
- **CPU**: 8+ cores
- **RAM**: 16GB+
- **Storage**: 100GB+ SSD
- **Network**: High-speed for data transfer

### Security Considerations

```bash
# Run as non-root user (already configured in container)
docker run --user 1000:1000 pandecay

# Read-only root filesystem
docker run --read-only --tmpfs /tmp pandecay

# Limit network access
docker run --network none pandecay

# Drop capabilities
docker run --cap-drop ALL pandecay
```

### Monitoring and Logging

```bash
# Monitor resource usage
docker stats pandecay

# Collect logs with timestamp
docker logs -t -f pandecay > pandecay.log

# Health check
docker inspect --format='{{.State.Health.Status}}' pandecay
```

### Backup and Recovery

```bash
# Backup output data
tar -czf pandecay-output-$(date +%Y%m%d).tar.gz output/

# Backup configuration
tar -czf pandecay-config-$(date +%Y%m%d).tar.gz config/

# Container image backup
docker save pandecay:latest | gzip > pandecay-image.tar.gz
```

## Development Workflow

### Interactive Development

```bash
# Start development container
./docker-deploy.sh dev

# Access container
docker exec -it pandecay-dev bash

# Inside container - modify code and test
cd /workspace
python3 panDecay.py examples/alignment.fas --debug
```

### Code Mounting

```bash
# Mount source code for live editing
docker run -it --rm \
  -v $(pwd):/workspace \
  -v $(pwd)/data:/data \
  -v $(pwd)/output:/output \
  pandecay bash

# Test changes immediately
python3 panDecay.py /data/alignment.fas --debug
```

### Building Custom Images

```bash
# Build with custom tag
docker build -t pandecay:custom .

# Build with specific Python version
docker build --build-arg PYTHON_VERSION=3.11 -t pandecay:py311 .

# Build without cache for clean build
docker build --no-cache -t pandecay:latest .
```

### Testing

```bash
# Run container tests
./docker-deploy.sh test

# Test specific functionality
docker run --rm pandecay test

# Test with your data
docker run --rm \
  -v $(pwd)/test_data:/data:ro \
  -v $(pwd)/test_output:/output \
  pandecay analyze /data/test_alignment.fas --debug
```

## Troubleshooting

### Common Issues

#### 1. Permission Problems

```bash
# Fix ownership of output files
sudo chown -R $USER:$USER output/

# Run container as current user
docker run --user $(id -u):$(id -g) pandecay

# Fix Docker group permissions
sudo usermod -aG docker $USER
newgrp docker
```

#### 2. Memory Issues

```bash
# Increase Docker memory limit (Docker Desktop)
# Settings → Resources → Memory → 8GB+

# Monitor memory usage
docker stats --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}"

# Reduce resource usage
docker run --memory="4g" --cpus="2.0" pandecay
```

#### 3. Volume Mount Issues

```bash
# Check volume mounts
docker inspect pandecay | grep -A 10 '"Mounts"'

# Verify paths exist
ls -la data/ output/ config/

# Use absolute paths
docker run -v /absolute/path/to/data:/data pandecay
```

#### 4. Network Connectivity

```bash
# Test internet connectivity
docker run --rm pandecay ping -c 3 google.com

# Check DNS resolution
docker run --rm pandecay nslookup github.com

# Use host network if needed
docker run --network host pandecay
```

#### 5. Container Won't Start

```bash
# Check container logs
docker logs pandecay

# Run with debug output
docker run -e PANDECAY_DEBUG=true pandecay

# Test with minimal command
docker run --rm pandecay test

# Check image integrity
docker inspect pandecay:latest
```

### Performance Optimization

#### CPU Optimization

```bash
# Match container CPUs to host
docker run --cpus="$(nproc)" pandecay

# Use specific CPU affinity
docker run --cpuset-cpus="0-3" pandecay

# Enable auto thread detection
docker run -e PANDECAY_THREADS=auto pandecay
```

#### Memory Optimization

```bash
# Set memory limits
docker run --memory="8g" --memory-swap="16g" pandecay

# Monitor memory usage
docker exec pandecay free -h

# Use memory-efficient settings
docker run \
  -e PANDECAY_FORMAT=static \
  -e PANDECAY_ASYNC=false \
  pandecay analyze large_alignment.fas
```

#### Storage Optimization

```bash
# Use SSD for Docker data
# Move /var/lib/docker to SSD mount

# Clean up Docker system
docker system prune -f

# Remove unused images
docker image prune -f

# Use multi-stage builds (already implemented)
```

### Getting Help

#### Check Logs

```bash
# Container logs
docker logs pandecay

# Application logs (if debug enabled)
docker exec pandecay cat /output/pandecay.log

# System logs
journalctl -u docker.service
```

#### Debug Mode

```bash
# Enable debug output
docker run -e PANDECAY_DEBUG=true pandecay analyze alignment.fas

# Keep temporary files
docker run -e PANDECAY_DEBUG=true pandecay analyze alignment.fas --keep-files

# Verbose Docker output
docker run --log-level debug pandecay
```

#### Support Resources

- **Container health**: `docker inspect --format='{{.State.Health}}' pandecay`
- **Resource usage**: `docker stats pandecay`
- **Network info**: `docker network ls`
- **Volume info**: `docker volume ls`
- **System info**: `docker system info`

### Container Maintenance

#### Updates

```bash
# Pull latest base image
docker pull ubuntu:22.04

# Rebuild container
./docker-deploy.sh build --no-cache

# Update Docker
sudo apt-get update && sudo apt-get upgrade docker-ce
```

#### Cleanup

```bash
# Stop all containers
docker stop $(docker ps -q)

# Remove all containers
docker container prune -f

# Remove all images
docker image prune -f

# Complete cleanup
./docker-deploy.sh clean
```

This guide covers the complete Docker workflow for panDecay. For additional help, consult the main README.md or open an issue on the project repository.