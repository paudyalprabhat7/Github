#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=public
#SBATCH --gres=gpu:2
#SBATCH --time=1-16:00:00
#SBATCH --output=cube_drop_2_9May.%j.out
#SBATCH --error=cube_drop_2_9May.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppaudya3@asu.edu

# Load necessary environment modules
module purge
module load cuda-11.7.0-gcc-11.2.0 cmake-3.23.1-gcc-11.2.0 gcc-11.2.0-gcc-11.2.0

# Check for CUDA availability
if ! command -v nvcc &> /dev/null; then
    echo "CUDA is not available on this node."
    exit 1
fi

# Display GPU information
nvidia-smi

# Run the executive
./bin/cube_drop_2