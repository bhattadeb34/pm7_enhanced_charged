import tempfile
import uuid
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
import os
import pandas as pd
import re

# Suppress RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

class ColabPM7Calculator:
    """PM7 calculator optimized for Google Colab environment with charge support."""

    def __init__(self, method="PM7"):
        """
        Initialize the calculator for Colab.

        Args:
            method: Semi-empirical method (PM7, PM6, etc.)
        """
        self.method = method
        self.temp_dir = "/tmp"  # Use /tmp in Colab

        # MOPAC keywords optimized for Colab
        self.keywords = "PRECISE GNORM=0.001 SCFCRT=1.D-8"

        # Proton reference energy (gas phase)
        self.proton_hof = 365.7  # kcal/mol

        # Check if MOPAC is available
        self._check_mopac()

    def _check_mopac(self):
        """Check if MOPAC is properly installed."""
        try:
            result = subprocess.run(["mopac"], capture_output=True, text=True)
            print("✓ MOPAC is available")
            return True
        except FileNotFoundError:
            print("✗ MOPAC not found. Please install using the installation cell above.")
            return False

    def smiles_to_3d(self, smiles, charge=0):
        """
        Convert SMILES to 3D coordinates with charge awareness.

        Args:
            smiles: SMILES string
            charge: Molecular charge (default: 0)

        Returns:
            tuple: (atoms_list, coordinates_array, formal_charge) or (None, None, None) if failed
        """
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, None, None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Calculate formal charge from molecule if not specified
            if charge == 'auto':
                charge = Chem.GetFormalCharge(mol)
                print(f"    Auto-detected charge: {charge:+d}")

            # Generate 3D coordinates
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useSmallRingTorsions = True

            status = AllChem.EmbedMolecule(mol, params)
            if status != 0:
                return None, None, None

            # Optimize with force field
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            except:
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
                except:
                    pass  # Use unoptimized structure

            # Extract atoms and coordinates
            atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
            conf = mol.GetConformer()
            coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

            return atoms, coords, charge

        except Exception as e:
            print(f"✗ Failed to generate 3D structure for {smiles}: {e}")
            return None, None, None

    def write_mopac_input(self, atoms, coordinates, label, charge=0, multiplicity=None):
        """
        Write MOPAC input file with charge and multiplicity support.
        
        Args:
            atoms: List of atom symbols
            coordinates: Array of 3D coordinates
            label: Label for the calculation
            charge: Molecular charge (default: 0)
            multiplicity: Spin multiplicity (default: None, auto-determined)
        """
        input_file = os.path.join(self.temp_dir, f"{label}.mop")

        # Build keyword line
        keyword_line = f"{self.method} {self.keywords}"
        
        # Add charge
        if charge != 0:
            keyword_line += f" CHARGE={charge}"
        
        # Add multiplicity/UHF for radicals and charged species
        if multiplicity is not None and multiplicity > 1:
            keyword_line += f" UHF DOUBLET" if multiplicity == 2 else f" UHF TRIPLET" if multiplicity == 3 else f" UHF"
        elif charge != 0:
            # For charged species, use UHF by default for better convergence
            keyword_line += " UHF"

        # Write MOPAC input
        with open(input_file, 'w') as f:
            f.write(f"{keyword_line}\n")
            f.write(f"{self.method} calculation for {label} (charge={charge:+d})\n")
            f.write("\n")

            for atom, coord in zip(atoms, coordinates):
                f.write(f"{atom:2s} {coord[0]:12.6f} 1 {coord[1]:12.6f} 1 {coord[2]:12.6f} 1\n")

        return input_file

    def run_mopac_calculation(self, input_file):
        """Run MOPAC calculation."""
        try:
            # Run MOPAC
            result = subprocess.run(
                ["mopac", input_file],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                print(f"⚠ MOPAC failed with return code {result.returncode}")
                print(f"Error: {result.stderr}")
                return False

            return True

        except subprocess.TimeoutExpired:
            print("⚠ MOPAC calculation timed out")
            return False
        except Exception as e:
            print(f"⚠ Error running MOPAC: {e}")
            return False

    def parse_mopac_output(self, output_file):
        """Parse MOPAC output file for properties with improved parsing."""
        properties = {}
    
        try:
            with open(output_file, 'r') as f:
                content = f.read()
                lines = content.split('\n')
    
            print(f"📊 Parsing MOPAC output: {output_file}")
    
            # 1. Parse heat of formation
            hof_pattern = r"FINAL\s+HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.\d+)\s*KCAL/MOL"
            hof_match = re.search(hof_pattern, content, re.IGNORECASE)
            if hof_match:
                properties['heat_of_formation'] = float(hof_match.group(1))
                print(f"    Heat of Formation: {properties['heat_of_formation']:.3f} kcal/mol")
            else:
                print("    Heat of Formation: Not found")
    
            # 2. Parse dipole moment
            dipole_pattern = r"SUM\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)"
            dipole_match = re.search(dipole_pattern, content)
            if dipole_match:
                properties['dipole_moment'] = float(dipole_match.group(4))
                properties['dipole_x'] = float(dipole_match.group(1))
                properties['dipole_y'] = float(dipole_match.group(2))
                properties['dipole_z'] = float(dipole_match.group(3))
                print(f"    Dipole Moment: {properties['dipole_moment']:.3f} Debye")
            else:
                print("    Dipole Moment: Not found")
    
            # 3. Parse HOMO/LUMO energies - Handle both closed and open shell systems
            homo_lumo_pattern = r"HOMO\s+LUMO\s+ENERGIES\s*\(EV\)\s*=\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)"
            homo_lumo_match = re.search(homo_lumo_pattern, content)
            
            if homo_lumo_match:
                # Closed shell system
                properties['homo_ev'] = float(homo_lumo_match.group(1))
                properties['lumo_ev'] = float(homo_lumo_match.group(2))
                properties['gap_ev'] = properties['lumo_ev'] - properties['homo_ev']
                properties['spin_state'] = 'closed_shell'
                print(f"    HOMO: {properties['homo_ev']:.3f} eV")
                print(f"    LUMO: {properties['lumo_ev']:.3f} eV")
                print(f"    HOMO-LUMO Gap: {properties['gap_ev']:.3f} eV")
            else:
                # Try SOMO/LUMO pattern (open shell)
                alpha_somo_pattern = r"ALPHA\s+SOMO\s+LUMO\s*\(EV\)\s*=\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)"
                beta_somo_pattern = r"BETA\s+SOMO\s+LUMO\s*\(EV\)\s*=\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)"
                
                alpha_match = re.search(alpha_somo_pattern, content)
                beta_match = re.search(beta_somo_pattern, content)
                
                if alpha_match and beta_match:
                    # Open shell system
                    properties['alpha_somo_ev'] = float(alpha_match.group(1))
                    properties['alpha_lumo_ev'] = float(alpha_match.group(2))
                    properties['beta_somo_ev'] = float(beta_match.group(1))
                    properties['beta_lumo_ev'] = float(beta_match.group(2))
                    
                    # For compatibility, set HOMO to the higher SOMO energy
                    properties['homo_ev'] = max(properties['alpha_somo_ev'], properties['beta_somo_ev'])
                    properties['lumo_ev'] = min(properties['alpha_lumo_ev'], properties['beta_lumo_ev'])
                    properties['gap_ev'] = properties['lumo_ev'] - properties['homo_ev']
                    properties['spin_state'] = 'open_shell'
                    
                    print(f"    Alpha SOMO: {properties['alpha_somo_ev']:.3f} eV")
                    print(f"    Beta SOMO: {properties['beta_somo_ev']:.3f} eV")
                    print(f"    Effective HOMO: {properties['homo_ev']:.3f} eV")
                    print(f"    Effective LUMO: {properties['lumo_ev']:.3f} eV")
                    print(f"    SOMO-LUMO Gap: {properties['gap_ev']:.3f} eV")
                else:
                    print("    HOMO/LUMO/SOMO energies: Not found")
    
            # 4. Parse ionization potential
            ip_pattern = r"IONIZATION\s+POTENTIAL\s*=\s*([-+]?\d+\.\d+)\s*EV"
            ip_match = re.search(ip_pattern, content, re.IGNORECASE)
            if ip_match:
                properties['ionization_potential'] = float(ip_match.group(1))
                print(f"    Ionization Potential: {properties['ionization_potential']:.3f} eV")
            else:
                print("    Ionization Potential: Not found")
    
            # 5. Parse COSMO area and volume
            cosmo_area_pattern = r"COSMO\s+AREA\s*=\s*([-+]?\d+\.\d+)\s*SQUARE\s+ANGSTROMS"
            cosmo_area_match = re.search(cosmo_area_pattern, content, re.IGNORECASE)
            if cosmo_area_match:
                properties['cosmo_area'] = float(cosmo_area_match.group(1))
                print(f"    COSMO Area: {properties['cosmo_area']:.2f} Ų")
    
            cosmo_volume_pattern = r"COSMO\s+VOLUME\s*=\s*([-+]?\d+\.\d+)\s*CUBIC\s+ANGSTROMS"
            cosmo_volume_match = re.search(cosmo_volume_pattern, content, re.IGNORECASE)
            if cosmo_volume_match:
                properties['cosmo_volume'] = float(cosmo_volume_match.group(1))
                print(f"    COSMO Volume: {properties['cosmo_volume']:.2f} ų")
    
            # 6. Parse molecular weight
            mw_pattern = r"MOLECULAR\s+WEIGHT\s*=\s*([-+]?\d+\.\d+)"
            mw_match = re.search(mw_pattern, content, re.IGNORECASE)
            if mw_match:
                properties['molecular_weight'] = float(mw_match.group(1))
                print(f"    Molecular Weight: {properties['molecular_weight']:.2f} g/mol")
    
            # 7. Parse point group
            pg_pattern = r"POINT\s+GROUP:\s*([A-Za-z0-9]+)"
            pg_match = re.search(pg_pattern, content, re.IGNORECASE)
            if pg_match:
                properties['point_group'] = pg_match.group(1)
                print(f"    Point Group: {properties['point_group']}")
    
            # 8. Parse number of filled levels
            filled_levels_pattern = r"NO\.\s+OF\s+FILLED\s+LEVELS\s*=\s*(\d+)"
            filled_levels_match = re.search(filled_levels_pattern, content, re.IGNORECASE)
            if filled_levels_match:
                properties['filled_levels'] = int(filled_levels_match.group(1))
                print(f"    Filled Levels: {properties['filled_levels']}")
    
            # 9. Parse electron counts and spin properties for open shell systems
            alpha_electrons_pattern = r"NO\.\s+OF\s+ALPHA\s+ELECTRONS\s*=\s*(\d+)"
            beta_electrons_pattern = r"NO\.\s+OF\s+BETA\s+ELECTRONS\s*=\s*(\d+)"
            
            alpha_match = re.search(alpha_electrons_pattern, content)
            beta_match = re.search(beta_electrons_pattern, content)
            
            if alpha_match and beta_match:
                alpha_electrons = int(alpha_match.group(1))
                beta_electrons = int(beta_match.group(1))
                unpaired_electrons = abs(alpha_electrons - beta_electrons)
                multiplicity = unpaired_electrons + 1
                
                properties['alpha_electrons'] = alpha_electrons
                properties['beta_electrons'] = beta_electrons
                properties['unpaired_electrons'] = unpaired_electrons
                properties['multiplicity'] = multiplicity
                
                print(f"    Alpha electrons: {alpha_electrons}")
                print(f"    Beta electrons: {beta_electrons}")
                print(f"    Multiplicity: {multiplicity}")
                if unpaired_electrons > 0:
                    print(f"    Unpaired electrons: {unpaired_electrons}")
    
            # 10. Parse spin contamination (S**2)
            spin_pattern = r"\(S\*\*2\)\s*=\s*([\d.]+)"
            spin_match = re.search(spin_pattern, content)
            if spin_match:
                properties['s_squared'] = float(spin_match.group(1))
                print(f"    S²: {properties['s_squared']:.3f}")
    
            # 11. Parse molecular charge
            charge_pattern = r"CHARGE\s+ON\s+SYSTEM\s*=\s*([-+]?\d+)"
            charge_match = re.search(charge_pattern, content, re.IGNORECASE)
            if charge_match:
                properties['charge'] = int(charge_match.group(1))
                print(f"    Molecular Charge: {properties['charge']:+d}")
    
            # 12. Calculate total energy if heat of formation is available
            if 'heat_of_formation' in properties:
                # Convert kcal/mol to eV (1 kcal/mol ≈ 0.043363 eV)
                properties['total_energy_ev'] = properties['heat_of_formation'] * 0.043363
                properties['total_energy_kcal_mol'] = properties['heat_of_formation']
                print(f"    Total Energy: {properties['total_energy_kcal_mol']:.3f} kcal/mol")
                print(f"    Total Energy: {properties['total_energy_ev']:.3f} eV")
    
            # 13. Parse computation time
            comp_time_pattern = r"COMPUTATION\s+TIME\s*=\s*([\d.]+)\s*SECONDS"
            comp_time_match = re.search(comp_time_pattern, content, re.IGNORECASE)
            if comp_time_match:
                properties['computation_time'] = float(comp_time_match.group(1))
                print(f"    Computation Time: {properties['computation_time']:.3f} seconds")
    
            # Summary
            found_properties = len([k for k, v in properties.items() if v is not None])
            print(f"    Successfully parsed {found_properties} properties")
    
        except Exception as e:
            print(f"✗ Error parsing MOPAC output: {e}")
            import traceback
            traceback.print_exc()
    
        return properties

    def get_temp_files(self, label):
        """Get list of temporary files for a given label."""
        extensions = ['.mop', '.out', '.arc', '.aux', '.log', '.end']
        temp_files = []
        for ext in extensions:
            temp_file = os.path.join(self.temp_dir, f"{label}{ext}")
            if os.path.exists(temp_file):
                temp_files.append(temp_file)
        return temp_files

    def cleanup_files(self, label):
        """Clean up temporary files."""
        extensions = ['.mop', '.out', '.arc', '.aux', '.log', '.end']
        cleaned_files = []
        for ext in extensions:
            temp_file = os.path.join(self.temp_dir, f"{label}{ext}")
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    cleaned_files.append(temp_file)
            except Exception as e:
                print(f"Warning: Could not remove {temp_file}: {e}")

        if cleaned_files:
            print(f"    🗑 Cleaned up {len(cleaned_files)} temporary files")

        return cleaned_files

    def calculate_properties(self, smiles, charge=0, multiplicity=None, cleanup=True):
        """
        Calculate PM7 properties for a SMILES string with charge support.

        Args:
            smiles: SMILES string
            charge: Molecular charge (default: 0, use 'auto' to detect from SMILES)
            multiplicity: Spin multiplicity (default: None, auto-determined)
            cleanup: Whether to clean up temporary files (default: True)

        Returns:
            dict: Properties dictionary with file information
        """
        label = f"mol_{uuid.uuid4().hex[:8]}"

        try:
            print(f"🔬 Processing: {smiles} (charge={charge if charge != 'auto' else 'auto'})")

            # Generate 3D structure
            atoms, coords, detected_charge = self.smiles_to_3d(smiles, charge)
            if atoms is None:
                return {'success': False, 'error': 'Failed to generate 3D structure', 'smiles': smiles}

            # Use detected charge if auto
            final_charge = detected_charge if charge == 'auto' else charge

            print(f"    Generated 3D structure ({len(atoms)} atoms, charge={final_charge:+d})")

            # Write MOPAC input
            input_file = self.write_mopac_input(atoms, coords, label, final_charge, multiplicity)
            print(f"    Created input file: {input_file}")

            # Run MOPAC calculation
            success = self.run_mopac_calculation(input_file)
            if not success:
                return {'success': False, 'error': 'MOPAC calculation failed', 'smiles': smiles, 'charge': final_charge}

            print(f"   ⚡ MOPAC calculation completed")

            # Parse output
            output_file = os.path.join(self.temp_dir, f"{label}.out")
            properties = self.parse_mopac_output(output_file)

            if not properties:
                return {'success': False, 'error': 'Failed to parse properties', 'smiles': smiles, 'charge': final_charge}

            # Add metadata
            properties['success'] = True
            properties['smiles'] = smiles
            properties['num_atoms'] = len(atoms)
            properties['label'] = label
            properties['input_charge'] = final_charge

            # Add file information
            temp_files = self.get_temp_files(label)
            properties['temp_files'] = temp_files
            properties['files_kept'] = not cleanup

            if cleanup:
                cleaned_files = self.cleanup_files(label)
                properties['cleaned_files'] = cleaned_files
            else:
                print(f"    Keeping {len(temp_files)} temporary files:")
                for temp_file in temp_files:
                    print(f"      - {temp_file}")

            print(f"    ✅ Properties calculated successfully")
            return properties

        except Exception as e:
            return {'success': False, 'error': str(e), 'smiles': smiles, 'charge': charge if charge != 'auto' else 0}

        finally:
            # Additional safety cleanup only if requested
            if cleanup:
                try:
                    self.cleanup_files(label)
                except:
                    pass

    def calculate_proton_affinity(self, smiles_neutral, smiles_protonated=None, cleanup=True):
        """
        Calculate gas-phase proton affinity using PM7.
        
        PA = -ΔH for:  AH⁺ → A + H⁺
        PA = HOF(A) + HOF(H⁺) - HOF(AH⁺)
        
        Args:
            smiles_neutral: SMILES of neutral molecule
            smiles_protonated: SMILES of protonated form (if None, will attempt auto-protonation)
            cleanup: Whether to clean up temporary files
            
        Returns:
            dict: Results including PA and component energies
        """
        print("=" * 70)
        print("🧪 PROTON AFFINITY CALCULATION")
        print("=" * 70)
        
        results = {
            'success': False,
            'smiles_neutral': smiles_neutral,
            'smiles_protonated': smiles_protonated
        }
        
        # Calculate neutral species
        print("\n1️⃣  Calculating NEUTRAL species:")
        props_neutral = self.calculate_properties(smiles_neutral, charge=0, cleanup=cleanup)
        
        if not props_neutral['success']:
            results['error'] = f"Failed to calculate neutral species: {props_neutral.get('error', 'Unknown')}"
            return results
        
        hof_neutral = props_neutral.get('heat_of_formation')
        if hof_neutral is None:
            results['error'] = "Heat of formation not found for neutral species"
            return results
        
        results['hof_neutral'] = hof_neutral
        results['properties_neutral'] = props_neutral
        
        # Handle protonated form
        if smiles_protonated is None:
            print("\n⚠️  No protonated SMILES provided. You need to provide the protonated structure.")
            results['error'] = "Protonated SMILES required. Please provide explicitly."
            return results
        
        # Calculate protonated species
        print(f"\n2️⃣  Calculating PROTONATED species:")
        props_protonated = self.calculate_properties(smiles_protonated, charge=1, cleanup=cleanup)
        
        if not props_protonated['success']:
            results['error'] = f"Failed to calculate protonated species: {props_protonated.get('error', 'Unknown')}"
            return results
        
        hof_protonated = props_protonated.get('heat_of_formation')
        if hof_protonated is None:
            results['error'] = "Heat of formation not found for protonated species"
            return results
        
        results['hof_protonated'] = hof_protonated
        results['properties_protonated'] = props_protonated
        
        # Calculate proton affinity
        hof_proton = self.proton_hof
        pa = hof_neutral + hof_proton - hof_protonated
        
        results['hof_proton'] = hof_proton
        results['proton_affinity_kcal_mol'] = pa
        results['proton_affinity_kj_mol'] = pa * 4.184  # Convert to kJ/mol
        results['proton_affinity_ev'] = pa * 0.043363  # Convert to eV
        results['success'] = True
        
        # Display results
        print("\n" + "=" * 70)
        print("📊 PROTON AFFINITY RESULTS")
        print("=" * 70)
        print(f"Neutral molecule:     {smiles_neutral}")
        print(f"  HOF(A):            {hof_neutral:>10.3f} kcal/mol")
        print(f"\nProtonated molecule:  {smiles_protonated}")
        print(f"  HOF(AH⁺):          {hof_protonated:>10.3f} kcal/mol")
        print(f"\nProton reference:")
        print(f"  HOF(H⁺):           {hof_proton:>10.3f} kcal/mol (gas phase)")
        print(f"\n" + "-" * 70)
        print(f"PROTON AFFINITY (PA): {pa:>10.2f} kcal/mol")
        print(f"                      {pa * 4.184:>10.2f} kJ/mol")
        print(f"                      {pa * 0.043363:>10.3f} eV")
        print("=" * 70)
        
        return results

# ============================================================================
# ENHANCED ONE-LINE FUNCTIONS FOR COLAB WITH CHARGE SUPPORT
# ============================================================================

# Global calculator instance
_colab_calculator = None

def calculate_pm7_properties_colab(smiles, charge=0, multiplicity=None, method="PM7", cleanup=True):
    """
    One-line function to calculate PM7 properties in Colab with charge support.

    Args:
        smiles: SMILES string
        charge: Molecular charge (default: 0, use 'auto' to detect from SMILES)
        multiplicity: Spin multiplicity (default: None)
        method: Semi-empirical method
        cleanup: Whether to remove temporary files (default: True)

    Returns:
        dict: Properties dictionary
    """
    global _colab_calculator

    if _colab_calculator is None:
        _colab_calculator = ColabPM7Calculator(method=method)

    return _colab_calculator.calculate_properties(smiles, charge=charge, multiplicity=multiplicity, cleanup=cleanup)

def calculate_pm7_batch_colab(smiles_list, charges=None, method="PM7", max_molecules=None, cleanup=True):
    """
    Calculate PM7 properties for multiple SMILES in Colab with charge support.

    Args:
        smiles_list: List of SMILES strings
        charges: List of charges (same length as smiles_list) or single charge for all, or None for neutral
        method: Semi-empirical method
        max_molecules: Maximum number to process (None = all)
        cleanup: Whether to remove temporary files (default: True)

    Returns:
        list: List of property dictionaries
    """
    calculator = ColabPM7Calculator(method=method)
    results = []

    # Limit number if specified
    if max_molecules:
        smiles_list = smiles_list[:max_molecules]

    # Handle charges
    if charges is None:
        charges = [0] * len(smiles_list)
    elif isinstance(charges, int):
        charges = [charges] * len(smiles_list)
    elif len(charges) != len(smiles_list):
        raise ValueError("Length of charges must match length of smiles_list")

    cleanup_msg = "with file cleanup" if cleanup else "keeping temporary files"
    print(f"🚀 Processing {len(smiles_list)} molecules {cleanup_msg}...")

    for i, (smiles, charge) in enumerate(zip(smiles_list, charges)):
        print(f"\n📍 Molecule {i+1}/{len(smiles_list)}")

        props = calculator.calculate_properties(smiles, charge=charge, cleanup=cleanup)
        results.append(props)

        # Show progress
        if props['success']:
            hof = props.get('heat_of_formation', 'N/A')
            dipole = props.get('dipole_moment', 'N/A')
            homo = props.get('homo_ev', 'N/A')
            chg = props.get('input_charge', 'N/A')
            print(f"    Charge: {chg:+d}, HOF: {hof}, Dipole: {dipole}, HOMO: {homo}")
        else:
            print(f"    ❌ Failed: {props['error']}")

    # Summary
    successful = sum(1 for r in results if r['success'])
    total_files = sum(len(r.get('temp_files', [])) for r in results if r['success'])

    print(f"\n✅ Summary: {successful}/{len(results)} successful calculations")
    if not cleanup and total_files > 0:
        print(f"📁 Total temporary files kept: {total_files}")

    return results

def calculate_proton_affinity_colab(smiles_neutral, smiles_protonated, method="PM7", cleanup=True):
    """
    One-line function to calculate proton affinity in Colab.
    
    Args:
        smiles_neutral: SMILES of neutral molecule
        smiles_protonated: SMILES of protonated form
        method: Semi-empirical method
        cleanup: Whether to remove temporary files
        
    Returns:
        dict: Results including PA and component energies
    """
    calculator = ColabPM7Calculator(method=method)
    return calculator.calculate_proton_affinity(smiles_neutral, smiles_protonated, cleanup=cleanup)

def display_properties_enhanced(props):
    """Display PM7 properties with enhanced formatting including charge."""
    if not props.get('success', False):
        print(f"❌ Calculation failed: {props.get('error', 'Unknown error')}")
        return

    charge = props.get('input_charge', 0)
    charge_str = f" (charge={charge:+d})" if charge != 0 else ""
    
    print(f"✨ PM7 Properties for {props.get('smiles', 'Unknown')}{charge_str}:")
    print("=" * 60)

    # Charge information
    if 'charge' in props:
        print(f"⚡ Molecular Charge: {props['charge']:+d}")

    # Core properties
    if 'heat_of_formation' in props:
        print(f"🔥 Heat of Formation: {props['heat_of_formation']:.3f} kcal/mol")

    if 'dipole_moment' in props:
        print(f"🧲 Dipole Moment: {props['dipole_moment']:.3f} Debye")
        if 'dipole_x' in props:
            print(f"   Components: X={props['dipole_x']:.3f}, Y={props['dipole_y']:.3f}, Z={props['dipole_z']:.3f}")

    # Electronic properties
    if 'spin_state' in props:
        print(f"🌀 Spin State: {props['spin_state']}")
    
    if 'multiplicity' in props:
        print(f"   Multiplicity: {props['multiplicity']}")
    
    if 'unpaired_electrons' in props and props['unpaired_electrons'] > 0:
        print(f"   Unpaired electrons: {props['unpaired_electrons']}")

    if 'homo_ev' in props and 'lumo_ev' in props:
        print(f"⚛️  HOMO Energy: {props['homo_ev']:.3f} eV")
        print(f"⚛️  LUMO Energy: {props['lumo_ev']:.3f} eV")
        print(f"📊 HOMO-LUMO Gap: {props['gap_ev']:.3f} eV")

    if 'ionization_potential' in props:
        print(f"⚡ Ionization Potential: {props['ionization_potential']:.3f} eV")

    # Molecular properties
    if 'molecular_weight' in props:
        print(f"⚖️  Molecular Weight: {props['molecular_weight']:.2f} g/mol")

    if 'point_group' in props:
        print(f"🔷 Point Group: {props['point_group']}")

    if 'filled_levels' in props:
        print(f"📈 Filled Levels: {props['filled_levels']}")

    # COSMO properties
    if 'cosmo_area' in props:
        print(f"🌊 COSMO Area: {props['cosmo_area']:.2f} Ų")

    if 'cosmo_volume' in props:
        print(f"📦 COSMO Volume: {props['cosmo_volume']:.2f} ų")

    # Computational info
    if 'computation_time' in props:
        print(f"⏱️  Computation Time: {props['computation_time']:.3f} seconds")

    print(f"🔬 Number of Atoms: {props.get('num_atoms', 'N/A')}")

    # File information
    if 'temp_files' in props:
        temp_files = props['temp_files']
        files_kept = props.get('files_kept', False)

        if files_kept and temp_files:
            print(f"\n📁 Temporary files kept ({len(temp_files)}):")
            for temp_file in temp_files:
                file_size = os.path.getsize(temp_file) if os.path.exists(temp_file) else 0
                print(f"  - {os.path.basename(temp_file)} ({file_size} bytes)")

    print()

def calculate_pm7_dataframe_colab(df, smiles_column='smiles', charge_column=None, cleanup=True):
    """
    Calculate PM7 properties for SMILES in a pandas DataFrame with charge support.
    
    Args:
        df: pandas DataFrame containing SMILES
        smiles_column: Column name containing SMILES strings
        charge_column: Column name containing charges (optional, defaults to 0)
        cleanup: Whether to remove temporary files
        
    Returns:
        pandas DataFrame with added property columns
    """
    if smiles_column not in df.columns:
        raise ValueError(f"Column '{smiles_column}' not found in DataFrame")
    
    # Get SMILES list
    smiles_list = df[smiles_column].tolist()
    
    # Get charges list
    if charge_column and charge_column in df.columns:
        charges = df[charge_column].tolist()
    else:
        charges = None
    
    # Calculate properties
    results = calculate_pm7_batch_colab(smiles_list, charges=charges, cleanup=cleanup)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Merge with original DataFrame
    df_with_props = df.copy()
    
    # Add property columns
    property_columns = ['heat_of_formation', 'dipole_moment', 'homo_ev', 'lumo_ev', 
                       'gap_ev', 'ionization_potential', 'molecular_weight', 
                       'point_group', 'charge', 'multiplicity', 'spin_state', 'success']
    
    for col in property_columns:
        if col in results_df.columns:
            df_with_props[f'pm7_{col}'] = results_df[col]
    
    return df_with_props
