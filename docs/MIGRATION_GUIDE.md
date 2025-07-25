# panDecay v1.1 Feature Guide

This guide helps you upgrade from panDecay v1.0 to v1.1, which introduces significant new features including YAML configuration, organized output structure, dual visualization, and Docker containerization.

## Overview of Changes

### 🆕 New Features in v1.1

- **YAML/TOML Configuration**: Modern configuration format with validation
- **Organized Output Structure**: Timestamp-based directory organization for better file management
- **Docker Containerization**: Production-ready containers with PAUP* and MrBayes
- **Enhanced Error Handling**: Improved diagnostics and graceful degradation

### Backward Compatibility

**All existing workflows continue to work unchanged.** The new features are additive and optional.

## Configuration Enhancement

### Modern Configuration Formats

panDecay v1.1 supports YAML and TOML configuration formats alongside the existing INI format. These new formats offer better structure, validation, and readability.

#### Generate Template Configs

```bash
# Generate example YAML configuration
python3 panDecay.py --generate-yaml-config example_config.yaml

# Generate example TOML configuration  
python3 panDecay.py --generate-toml-config example_config.toml
```

## Performance Optimization

### Resource Management

```yaml
computational:
  threads: auto              # Auto-detect optimal thread count
  memory_limit: 8G          # Optional memory limit
```

## Visualization Enhancements

### Dual Visualization System

Choose between static plots (publications) and interactive plots (exploration):

```bash
# Static plots only (matplotlib) - good for publications
python3 panDecay.py alignment.fas --viz-format static

# All visualization is now matplotlib-based static plots
# Output is organized in timestamped directories

# Both formats (default) - maximum flexibility
python3 panDecay.py alignment.fas --viz-format both
```

### New Visualization Features

- **Interactive decay plots** with zoom, pan, and hover information
- **Statistical dashboards** with multiple subplot views
- **Publication-ready static plots** with enhanced styling
- **Animated constraint analysis** showing step-by-step results

```yaml
visualization:
  format: both                    # static, interactive, or both
  static_format: [png, pdf]      # Static output formats
  interactive_format: html       # Interactive output format
  dpi: 300                       # High resolution for publications
  theme: publication             # Styling theme
```

## Docker Deployment

### Quick Start with Docker

```bash
# Build the container
./docker-deploy.sh build

# Run analysis with Docker
./docker-deploy.sh run alignment.fas --analysis ml --viz-format both

# Start development environment
./docker-deploy.sh dev

# Test container functionality
./docker-deploy.sh test
```

### Docker Compose Deployment

```bash
# Copy and customize environment
cp .env.example .env

# Start production service
docker-compose up -d pandecay

# Start development environment
COMPOSE_PROFILES=development docker-compose up -d pandecay-dev

# Batch processing
COMPOSE_PROFILES=batch docker-compose up pandecay-batch
```

### Container Usage Examples

```bash
# Basic analysis
docker run -v $(pwd)/data:/data -v $(pwd)/output:/output \
  pandecay analyze /data/alignment.fas --analysis ml

# With YAML configuration
docker run -v $(pwd):/workspace \
  pandecay analyze --config /workspace/config.yaml

# Interactive shell
docker run -it -v $(pwd):/workspace pandecay bash

# Generate configuration template
docker run -v $(pwd):/workspace \
  pandecay config yaml /workspace/config.yaml
```

## Command Line Changes

### New Arguments

```bash
# Enhanced configuration
--config config.yaml            # YAML/TOML configuration file
--generate-yaml-config file     # Generate YAML template
--generate-toml-config file     # Generate TOML template
--generate-ini-config file      # Generate INI template

# Enhanced visualization  
--viz-format both               # static, interactive, both
--static-formats png,pdf        # Static output formats
--interactive-format html       # Interactive format
```

### Existing Arguments (Unchanged)

All existing command-line arguments work exactly as before:

```bash
# These commands work identically in v1.0 and v1.1
python3 panDecay.py alignment.fas --analysis ml --model GTR
python3 panDecay.py alignment.fas --config old_config.ini --threads 8
python3 panDecay.py alignment.fas --analysis bayesian --debug
```

## Troubleshooting Upgrade

### Common Issues

#### 1. Missing Dependencies
```bash
# Install new dependencies
pip3 install pydantic pyyaml toml

# Or from requirements
pip3 install -r requirements.txt
```

#### 2. Configuration Validation Errors
```bash
# Validate YAML configuration
python3 config_loader.py validate config.yaml

# Check for common issues
python3 panDecay.py --config config.yaml --dry-run
```

#### 3. Docker Permission Issues
```bash
# Fix Docker permission issues
sudo chmod +x docker-deploy.sh
sudo chown $USER:$USER -R output/

# Run as current user
docker run --user $(id -u):$(id -g) -v $(pwd):/workspace pandecay
```

### Performance Tuning

#### Memory Usage
```yaml
computational:
  threads: 4                    # Limit thread count if needed
  memory_limit: 4G             # Set memory limit for constrained systems
```

#### Large Dataset Optimization
```yaml
analysis:
  site_specific: false          # Disable for very large alignments
  constraint_strategy: fast     # Use faster constraint generation

visualization:
  format: static               # Skip interactive plots for speed
  generate_trees: false        # Skip tree generation if not needed
```

## Best Practices

### Configuration Management
- Use YAML for new projects (better structure and validation)
- Keep INI configs for existing automated pipelines
- Store configurations in version control
- Use different configs for development vs. production

### Performance Optimization  
- Use Docker for consistent environments
- Monitor resource usage with large datasets
- Consider batch processing for multiple analyses
- Optimize thread usage based on system capabilities

### Visualization Strategy
- Use interactive plots for exploration and data checking
- Use static plots for publications and reports
- Generate both formats for maximum flexibility
- Customize themes for different audiences

### Docker Deployment
- Use docker-compose for production deployments
- Mount data directories as read-only when possible
- Set appropriate resource limits
- Use environment variables for configuration

## Getting Help

### Documentation
- `README.md` - Updated with v1.1 features
- `DOCKER_GUIDE.md` - Comprehensive Docker documentation
- `CONFIGURATION_GUIDE.md` - Detailed configuration reference

### Testing New Features
```bash
# Test dual visualization
python3 panDecay.py examples/alignment.fas --viz-format both --debug

# Test Docker container
./docker-deploy.sh test

# Validate configuration
python3 config_loader.py validate config.yaml
```

### Support
- Check `CHANGELOG.md` for detailed changes
- Review example configurations in `examples/`
- Use `--debug` mode for troubleshooting
- Test with provided example datasets