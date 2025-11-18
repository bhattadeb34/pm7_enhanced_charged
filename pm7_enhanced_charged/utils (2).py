"""Utility functions for PM7 Enhanced Charged calculator"""
import subprocess
import sys

def install_colab_dependencies():
    """Install dependencies in Google Colab"""
    print("🔧 Installing PM7 Enhanced dependencies for Colab...")
    
    # Install condacolab if not already installed
    try:
        import condacolab
        print("✓ condacolab already available")
    except ImportError:
        print("📦 Installing condacolab...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "condacolab"])
        import condacolab
        condacolab.install()
    
    # Install MOPAC and dependencies via conda
    print("📦 Installing MOPAC and dependencies...")
    subprocess.check_call([
        "conda", "install", "-c", "conda-forge", 
        "mopac", "rdkit", "ase", "pandas", "numpy", 
        "scikit-learn", "matplotlib", "seaborn", "-y"
    ])
    
    print("✅ All dependencies installed successfully!")
    print("⚠️  IMPORTANT: Please restart your runtime now (Runtime > Restart runtime)")

def check_colab_environment():
    """Check if running in Google Colab"""
    try:
        import google.colab
        return True
    except ImportError:
        return False
