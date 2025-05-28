# Exit if any command fails
set -e

# Name conda environment
ENV_NAME="uht-tooling"

# Set conda subdir to osx-64
export CONDA_SUBDIR=osx-64

echo "[1/4] Creating conda environment: $ENV_NAME with Python 3.10..."
conda create --yes --name "$ENV_NAME" python=3.10

echo "[2/4] Activating environment: $ENV_NAME..."
# Use conda's shell hook to activate the environment
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

echo "[3/4] Installing Python requirements via pip..."
pip install -r requirements.txt

echo "[4/4] Installing accessories..."
conda install --yes -c bioconda mafft
conda install --yes -c bioconda minimap2

echo "Setup complete. To activate your environment later, run:"
echo "   conda activate $ENV_NAME"
