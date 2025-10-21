# Exit if any command fails
set -e

# Name conda environment
ENV_NAME="uht-tooling"

# Set conda subdir to osx-64
export CONDA_SUBDIR=osx-64

echo "[1/5] Creating conda environment: $ENV_NAME with Python 3.10..."
conda create --yes --name "$ENV_NAME" python=3.10

echo "[2/5] Activating environment: $ENV_NAME..."
# Use conda's shell hook to activate the environment
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

echo "[3/5] Installing Python requirements via pip..."
pip install -r requirements.txt

echo "[4/5] Installing accessories..."
conda install --yes -c bioconda mafft
conda install --yes -c bioconda minimap2
conda install --yes -c bioconda nanofilt

echo "[5/5] Setup complete!"
echo ""
echo "To activate your environment later, run:"
echo "   conda activate $ENV_NAME"
