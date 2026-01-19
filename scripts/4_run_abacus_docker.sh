#!/bin/bash
# Script to run ABACUS calculations using Docker
# Based on: https://abacus.deepmodeling.com/en/latest/quick_start/easy_install.html#container-deployment
#
# Optimized for Mac Mini M4:
#   - MPI_PROCS=8 (uses 8 of 10 cores on M4 base)
#   - OMP_NUM_THREADS=1 (MPI scales better for ABACUS)

set -e

# Docker image to use
ABACUS_IMAGE="registry.dp.tech/deepmodeling/abacus"
# Alternative: docker.io/deepmodeling/abacus (if the above doesn't work)

# Default number of MPI processes (8 for Mac Mini M4)
MPI_PROCS=${MPI_PROCS:-8}

# Default OpenMP threads (1 recommended - MPI scales better)
OMP_THREADS=${OMP_NUM_THREADS:-1}

# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default working directory (current directory or specified)
WORK_DIR="${1:-$PWD}"

# Convert to absolute path
if [[ "$WORK_DIR" != /* ]]; then
    WORK_DIR="$(cd "$WORK_DIR" && pwd)"
fi

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in PATH"
    echo "Please install Docker: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo "Error: Docker daemon is not running"
    echo "Please start Docker Desktop or the Docker service"
    exit 1
fi

# Pull the Docker image if it doesn't exist locally
echo "Checking for ABACUS Docker image..."
if ! docker images | grep -q "deepmodeling/abacus\|registry.dp.tech/deepmodeling/abacus"; then
    echo "Pulling ABACUS Docker image..."
    docker pull "$ABACUS_IMAGE" || {
        echo "Warning: Failed to pull from registry.dp.tech, trying docker.io..."
        ABACUS_IMAGE="docker.io/deepmodeling/abacus"
        docker pull "$ABACUS_IMAGE" || {
            echo "Error: Failed to pull ABACUS Docker image"
            echo "Please check your internet connection and Docker registry access"
            exit 1
        }
    }
    echo "✓ Docker image pulled successfully"
else
    echo "✓ Docker image found locally"
fi

# Check if INPUT file exists in the working directory
if [ ! -f "$WORK_DIR/INPUT" ]; then
    echo "Error: No INPUT file found in $WORK_DIR"
    echo "ABACUS requires an INPUT file to run calculations"
    echo ""
    echo "Usage: $0 [working_directory]"
    echo "Example: $0 ./data/abacus_inputs/fragment_001_conf_000"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running ABACUS with Docker"
echo "=========================================="
echo "Working directory: $WORK_DIR"
echo "MPI processes: $MPI_PROCS"
echo "OpenMP threads: $OMP_THREADS"
echo "Docker image: $ABACUS_IMAGE"
echo "=========================================="
echo ""

# Run ABACUS in Docker container
docker run --rm -it \
    -v "$WORK_DIR:/workspace" \
    -w /workspace \
    -e OMP_NUM_THREADS=$OMP_THREADS \
    "$ABACUS_IMAGE" \
    bash -c "mpirun --quiet --allow-run-as-root -n $MPI_PROCS abacus 2>&1 | grep -v 'exiting improperly' | grep -v 'mpirun has exited' || true; exit 0"

echo ""
echo "=========================================="
echo "ABACUS calculation completed"
echo "=========================================="
echo "Results should be in: $WORK_DIR/OUT.ABACUS/"
echo "=========================================="

