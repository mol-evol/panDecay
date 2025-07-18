#!/bin/bash

# panDecay Docker Deployment Script
# Provides easy commands for building, running, and managing panDecay containers

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE_NAME="pandecay"
CONTAINER_NAME="pandecay-main"

# Helper functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

show_help() {
    cat << EOF
panDecay Docker Deployment Script
=================================

Usage: $0 [COMMAND] [OPTIONS]

Commands:
  build          Build the panDecay Docker image
  run            Run a single analysis
  dev            Start development environment
  batch          Run batch processing
  test           Test the container
  shell          Open interactive shell
  clean          Clean up containers and images
  logs           Show container logs
  status         Show container status
  stop           Stop running containers
  help           Show this help message

Options:
  --force        Force rebuild/restart
  --debug        Enable debug mode
  --no-cache     Build without cache
  --profile PROF Use specific Docker Compose profile

Examples:
  # Build the container
  $0 build

  # Run analysis with default settings
  $0 run alignment.fas --analysis ml

  # Start development environment
  $0 dev

  # Run with custom config
  $0 run --config myconfig.yaml

  # Test container functionality
  $0 test

  # View logs
  $0 logs

Environment Variables:
  DATA_DIR       Input data directory (default: ./data)
  OUTPUT_DIR     Output directory (default: ./output)
  CONFIG_DIR     Configuration directory (default: ./config)

EOF
}

# Check dependencies
check_dependencies() {
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi

    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        log_error "Docker Compose is not installed or not in PATH"
        exit 1
    fi
}

# Build Docker image
build_image() {
    local force=false
    local no_cache=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            --force)
                force=true
                shift
                ;;
            --no-cache)
                no_cache=true
                shift
                ;;
            *)
                break
                ;;
        esac
    done

    log_info "Building panDecay Docker image..."

    local build_args=""
    if [[ "$no_cache" == "true" ]]; then
        build_args="--no-cache"
    fi

    if [[ "$force" == "true" ]] && docker images -q "$IMAGE_NAME" &> /dev/null; then
        log_info "Removing existing image..."
        docker rmi "$IMAGE_NAME" || true
    fi

    docker build $build_args -t "$IMAGE_NAME:latest" -f "$SCRIPT_DIR/Dockerfile" "$PROJECT_ROOT"
    log_success "Docker image built successfully"
}

# Run single analysis
run_analysis() {
    local data_dir="${DATA_DIR:-$SCRIPT_DIR/data}"
    local output_dir="${OUTPUT_DIR:-$SCRIPT_DIR/output}"
    local config_dir="${CONFIG_DIR:-$SCRIPT_DIR/config}"

    # Create directories if they don't exist
    mkdir -p "$data_dir" "$output_dir" "$config_dir"

    log_info "Running panDecay analysis..."
    log_info "Data directory: $data_dir"
    log_info "Output directory: $output_dir"
    log_info "Config directory: $config_dir"

    docker run --rm \
        -v "$data_dir:/data:ro" \
        -v "$output_dir:/output" \
        -v "$config_dir:/workspace/config:ro" \
        -e PANDECAY_THREADS="${PANDECAY_THREADS:-auto}" \
        -e PANDECAY_DEBUG="${PANDECAY_DEBUG:-false}" \
        -e PANDECAY_FORMAT="${PANDECAY_FORMAT:-both}" \
        "$IMAGE_NAME:latest" analyze "$@"
}

# Start development environment
start_dev() {
    log_info "Starting development environment..."
    COMPOSE_PROFILES=development docker-compose up -d pandecay-dev
    log_success "Development environment started"
    log_info "Access with: docker exec -it pandecay-dev bash"
}

# Run batch processing
run_batch() {
    log_info "Starting batch processing..."
    COMPOSE_PROFILES=batch docker-compose up pandecay-batch
}

# Test container
test_container() {
    log_info "Testing panDecay container..."
    
    # Test basic functionality
    docker run --rm "$IMAGE_NAME:latest" test
    
    # Test with example data if available
    if [[ -f "$SCRIPT_DIR/examples/test_alignment.fas" ]]; then
        log_info "Running test analysis..."
        docker run --rm \
            -v "$SCRIPT_DIR/examples:/data:ro" \
            -v "/tmp/pandecay_test:/output" \
            "$IMAGE_NAME:latest" analyze /data/test_alignment.fas --analysis ml --debug
        log_success "Test analysis completed"
    else
        log_warning "No test data found, skipping analysis test"
    fi
}

# Open interactive shell
open_shell() {
    local data_dir="${DATA_DIR:-$SCRIPT_DIR/data}"
    local output_dir="${OUTPUT_DIR:-$SCRIPT_DIR/output}"

    mkdir -p "$data_dir" "$output_dir"

    log_info "Opening interactive shell..."
    docker run --rm -it \
        -v "$data_dir:/data" \
        -v "$output_dir:/output" \
        -v "$SCRIPT_DIR:/workspace" \
        "$IMAGE_NAME:latest" bash
}

# Clean up containers and images
cleanup() {
    log_info "Cleaning up panDecay containers and images..."
    
    # Stop and remove containers
    docker-compose down -v || true
    docker stop "$CONTAINER_NAME" 2>/dev/null || true
    docker rm "$CONTAINER_NAME" 2>/dev/null || true
    
    # Remove images
    docker rmi "$IMAGE_NAME:latest" 2>/dev/null || true
    docker system prune -f
    
    log_success "Cleanup completed"
}

# Show logs
show_logs() {
    docker-compose logs -f pandecay
}

# Show status
show_status() {
    log_info "Container status:"
    docker-compose ps
    
    echo ""
    log_info "Image information:"
    docker images | grep "$IMAGE_NAME" || log_warning "No panDecay images found"
}

# Stop containers
stop_containers() {
    log_info "Stopping panDecay containers..."
    docker-compose down
    log_success "Containers stopped"
}

# Main command dispatcher
main() {
    check_dependencies

    case "${1:-help}" in
        build)
            shift
            build_image "$@"
            ;;
        run)
            shift
            run_analysis "$@"
            ;;
        dev|development)
            start_dev
            ;;
        batch)
            run_batch
            ;;
        test)
            test_container
            ;;
        shell|bash)
            open_shell
            ;;
        clean|cleanup)
            cleanup
            ;;
        logs)
            show_logs
            ;;
        status)
            show_status
            ;;
        stop)
            stop_containers
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            log_error "Unknown command: $1"
            echo ""
            show_help
            exit 1
            ;;
    esac
}

# Run main function with all arguments
main "$@"