#!/bin/bash
# Batch script to run ABACUS calculations on multiple conformations using Docker

set -e

# If not already running under caffeinate, re-execute with caffeinate
if [ -z "$CAFFEINATE_ACTIVE" ]; then
    export CAFFEINATE_ACTIVE=1
    exec caffeinate -i -d "$0" "$@"
fi

# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_SCRIPT="$SCRIPT_DIR/4_run_abacus_docker.sh"

# Parse command line arguments
BASE_DIR="${1:-data/abacus_inputs}"
MAX_CONFORMATIONS="${2:-}"
MPI_ARG="${3:-}"
OMP_ARG="${4:-}"

# Set MPI processes (argument > env var > default)
if [ -n "$MPI_ARG" ]; then
    export MPI_PROCS="$MPI_ARG"
else
    export MPI_PROCS=${MPI_PROCS:-4}
fi

# Set OpenMP threads (argument > env var > default)
if [ -n "$OMP_ARG" ]; then
    export OMP_NUM_THREADS="$OMP_ARG"
    OMP_THREADS="$OMP_ARG"
else
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
    OMP_THREADS=${OMP_NUM_THREADS:-1}
fi

# Limit number of conformations if specified
if [ -n "$MAX_CONFORMATIONS" ]; then
    if ! [[ "$MAX_CONFORMATIONS" =~ ^[0-9]+$ ]]; then
        echo "Error: MAX_CONFORMATIONS must be a number"
        echo "Usage: $0 [base_directory] [max_conformations] [mpi_procs] [omp_threads]"
        echo "Example: $0 data/abacus_inputs 50 8 4"
        exit 1
    fi
fi

# Check if base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory '$BASE_DIR' does not exist"
    echo "Usage: $0 [base_directory] [max_conformations] [mpi_procs] [omp_threads]"
    echo "Example: $0 data/abacus_inputs 50 8 4"
    exit 1
fi

# Check if Docker script exists
if [ ! -f "$DOCKER_SCRIPT" ]; then
    echo "Error: Docker script not found: $DOCKER_SCRIPT"
    exit 1
fi

echo "=========================================="
echo "Batch ABACUS Calculations with Docker"
echo "=========================================="
echo "Base directory: $BASE_DIR"
echo "MPI processes: $MPI_PROCS"
echo "OpenMP threads: $OMP_THREADS"
if [ -n "$MAX_CONFORMATIONS" ]; then
    echo "Max conformations: $MAX_CONFORMATIONS"
fi
echo "=========================================="
echo ""

# Find all conformation directories (directories with INPUT file)
CONF_DIRS=($(find "$BASE_DIR" -type f -name "INPUT" -exec dirname {} \; | sort))

if [ ${#CONF_DIRS[@]} -eq 0 ]; then
    echo "Error: No INPUT files found in $BASE_DIR"
    echo "Expected structure:"
    echo "  $BASE_DIR/"
    echo "    fragment_XXX_conf_000/"
    echo "      INPUT"
    echo "      STRU"
    echo "      KPT"
    echo "      ..."
    exit 1
fi

# Limit if specified
if [ -n "$MAX_CONFORMATIONS" ] && [ ${#CONF_DIRS[@]} -gt $MAX_CONFORMATIONS ]; then
    CONF_DIRS=("${CONF_DIRS[@]:0:$MAX_CONFORMATIONS}")
fi

TOTAL=${#CONF_DIRS[@]}
SUCCESS=0
FAILED=0

echo "Found $TOTAL conformations to process"
echo ""
echo "Note: Running with caffeinate to prevent system sleep"
echo ""

# Run ABACUS on each conformation
for i in "${!CONF_DIRS[@]}"; do
    conf_dir="${CONF_DIRS[$i]}"
    conf_name=$(basename "$conf_dir")
    current=$((i + 1))
    
    echo "[$current/$TOTAL] Processing $conf_name..."
    
    # Check if INPUT file exists
    if [ ! -f "$conf_dir/INPUT" ]; then
        echo "  ✗ Skipping: No INPUT file found"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Run ABACUS using the Docker script
    "$DOCKER_SCRIPT" "$conf_dir" > /tmp/abacus_batch_${current}.log 2>&1
    DOCKER_EXIT_CODE=$?
    
    # Check if output directory was created and has required files
    if [ -d "$conf_dir/OUT.ABACUS" ] && [ -f "$conf_dir/OUT.ABACUS/MD_dump" ]; then
        SUCCESS=$((SUCCESS + 1))
        echo "  ✓ Completed: $conf_name (MD_dump found)"
    elif [ -d "$conf_dir/OUT.ABACUS" ] && ls "$conf_dir/OUT.ABACUS"/running_*.log > /dev/null 2>&1; then
        # Has log but no MD_dump - might be incomplete or failed
        FAILED=$((FAILED + 1))
        echo "  ✗ Failed: $conf_name (no MD_dump file - check log)"
        echo "    Docker exit code: $DOCKER_EXIT_CODE"
        if [ -f /tmp/abacus_batch_${current}.log ]; then
            echo "    Last 5 lines of log:"
            tail -5 /tmp/abacus_batch_${current}.log | sed 's/^/      /'
        fi
    else
        FAILED=$((FAILED + 1))
        echo "  ✗ Failed: $conf_name (no output directory)"
        echo "    Docker exit code: $DOCKER_EXIT_CODE"
        if [ -f /tmp/abacus_batch_${current}.log ]; then
            echo "    Last 10 lines of log:"
            tail -10 /tmp/abacus_batch_${current}.log | sed 's/^/      /'
        fi
    fi
    
    echo ""
done

# Summary
echo "=========================================="
echo "Summary:"
echo "  Total: $TOTAL"
echo "  Success: $SUCCESS"
echo "  Failed: $FAILED"
echo "=========================================="

if [ $FAILED -eq 0 ]; then
    echo ""
    echo "✓ All calculations completed successfully!"
    exit 0
else
    echo ""
    echo "⚠ Some calculations failed. Check logs above for details."
    exit 1
fi

