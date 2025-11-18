"""
PM7 Calculator Usage Examples with Charge Support
==================================================

This script demonstrates how to use the enhanced PM7 calculator with:
1. Neutral molecules
2. Charged species (cations and anions)
3. Proton affinity calculations
"""

from pm7calculator import (
    calculate_pm7_properties_colab,
    calculate_pm7_batch_colab,
    calculate_proton_affinity_colab,
    display_properties_enhanced,
    ColabPM7Calculator
)

# ============================================================================
# EXAMPLE 1: Basic Neutral Molecule Calculation
# ============================================================================

def example_neutral_molecule():
    """Calculate properties for a neutral molecule"""
    print("\n" + "="*70)
    print("EXAMPLE 1: Neutral Molecule (Water)")
    print("="*70)
    
    smiles = "O"  # Water
    props = calculate_pm7_properties_colab(smiles, charge=0)
    display_properties_enhanced(props)
    
    return props

# ============================================================================
# EXAMPLE 2: Cation Calculation
# ============================================================================

def example_cation():
    """Calculate properties for a positively charged species"""
    print("\n" + "="*70)
    print("EXAMPLE 2: Cation (Protonated Water, H3O+)")
    print("="*70)
    
    smiles = "[OH3+]"  # Hydronium ion
    props = calculate_pm7_properties_colab(smiles, charge=1)
    display_properties_enhanced(props)
    
    return props

# ============================================================================
# EXAMPLE 3: Anion Calculation
# ============================================================================

def example_anion():
    """Calculate properties for a negatively charged species"""
    print("\n" + "="*70)
    print("EXAMPLE 3: Anion (Hydroxide, OH-)")
    print("="*70)
    
    smiles = "[OH-]"  # Hydroxide ion
    props = calculate_pm7_properties_colab(smiles, charge=-1)
    display_properties_enhanced(props)
    
    return props

# ============================================================================
# EXAMPLE 4: Batch Calculation with Mixed Charges
# ============================================================================

def example_batch_mixed_charges():
    """Calculate properties for multiple molecules with different charges"""
    print("\n" + "="*70)
    print("EXAMPLE 4: Batch Calculation with Mixed Charges")
    print("="*70)
    
    smiles_list = [
        "O",           # Water (neutral)
        "[OH3+]",      # Hydronium (cation)
        "[OH-]",       # Hydroxide (anion)
        "N",           # Ammonia (neutral)
        "[NH4+]",      # Ammonium (cation)
    ]
    
    charges = [0, 1, -1, 0, 1]
    
    results = calculate_pm7_batch_colab(smiles_list, charges=charges)
    
    # Display summary
    print("\n" + "="*70)
    print("SUMMARY OF RESULTS")
    print("="*70)
    for i, result in enumerate(results):
        if result['success']:
            hof = result.get('heat_of_formation', 'N/A')
            chg = result.get('input_charge', 'N/A')
            print(f"{i+1}. {smiles_list[i]:20s} | Charge: {chg:+2d} | HOF: {hof:8.2f} kcal/mol")
        else:
            print(f"{i+1}. {smiles_list[i]:20s} | FAILED: {result.get('error')}")
    
    return results

# ============================================================================
# EXAMPLE 5: Proton Affinity Calculation - Water
# ============================================================================

def example_proton_affinity_water():
    """Calculate proton affinity of water"""
    print("\n" + "="*70)
    print("EXAMPLE 5: Proton Affinity of Water")
    print("="*70)
    print("\nReaction: H2O + H+ → H3O+")
    print("Expected PA (experimental): ~165 kcal/mol")
    
    smiles_neutral = "O"       # Water
    smiles_protonated = "[OH3+]"  # Hydronium
    
    result = calculate_proton_affinity_colab(smiles_neutral, smiles_protonated)
    
    if result['success']:
        pa = result['proton_affinity_kcal_mol']
        print(f"\n🎯 Calculated PA: {pa:.2f} kcal/mol")
        print(f"📊 Typical PM7 error: ±3-8 kcal/mol")
    
    return result

# ============================================================================
# EXAMPLE 6: Proton Affinity Calculation - Ammonia
# ============================================================================

def example_proton_affinity_ammonia():
    """Calculate proton affinity of ammonia"""
    print("\n" + "="*70)
    print("EXAMPLE 6: Proton Affinity of Ammonia")
    print("="*70)
    print("\nReaction: NH3 + H+ → NH4+")
    print("Expected PA (experimental): ~204 kcal/mol")
    
    smiles_neutral = "N"       # Ammonia
    smiles_protonated = "[NH4+]"  # Ammonium
    
    result = calculate_proton_affinity_colab(smiles_neutral, smiles_protonated)
    
    if result['success']:
        pa = result['proton_affinity_kcal_mol']
        print(f"\n🎯 Calculated PA: {pa:.2f} kcal/mol")
        print(f"📊 Typical PM7 error: ±3-8 kcal/mol")
    
    return result

# ============================================================================
# EXAMPLE 7: Proton Affinity for Organic Bases
# ============================================================================

def example_proton_affinity_organic_bases():
    """Calculate proton affinities for a series of organic bases"""
    print("\n" + "="*70)
    print("EXAMPLE 7: Proton Affinities of Organic Bases")
    print("="*70)
    
    molecules = [
        ("Methylamine", "CN", "C[NH3+]"),
        ("Dimethylamine", "CNC", "C[NH2+]C"),
        ("Trimethylamine", "CN(C)C", "C[NH+](C)C"),
        ("Pyridine", "c1cccnc1", "c1ccc[nH+]c1"),
    ]
    
    calculator = ColabPM7Calculator()
    results = []
    
    print("\nCalculating proton affinities...")
    for name, neutral_smiles, protonated_smiles in molecules:
        print(f"\n{name}:")
        result = calculator.calculate_proton_affinity(neutral_smiles, protonated_smiles, cleanup=True)
        results.append((name, result))
    
    # Display comparison
    print("\n" + "="*70)
    print("PROTON AFFINITY COMPARISON")
    print("="*70)
    print(f"{'Molecule':<20} {'PA (kcal/mol)':<15} {'PA (kJ/mol)':<15}")
    print("-"*70)
    
    for name, result in results:
        if result['success']:
            pa_kcal = result['proton_affinity_kcal_mol']
            pa_kj = result['proton_affinity_kj_mol']
            print(f"{name:<20} {pa_kcal:>12.2f}   {pa_kj:>12.2f}")
        else:
            print(f"{name:<20} FAILED")
    
    return results

# ============================================================================
# EXAMPLE 8: DataFrame with Charged Species
# ============================================================================

def example_dataframe_with_charges():
    """Calculate properties for molecules in a DataFrame with charges"""
    import pandas as pd
    from pm7calculator import calculate_pm7_dataframe_colab
    
    print("\n" + "="*70)
    print("EXAMPLE 8: DataFrame with Charged Species")
    print("="*70)
    
    # Create example DataFrame
    data = {
        'name': ['Water', 'Hydronium', 'Hydroxide', 'Ammonia', 'Ammonium'],
        'smiles': ['O', '[OH3+]', '[OH-]', 'N', '[NH4+]'],
        'charge': [0, 1, -1, 0, 1],
        'type': ['neutral', 'cation', 'anion', 'neutral', 'cation']
    }
    
    df = pd.DataFrame(data)
    
    print("\nInput DataFrame:")
    print(df)
    
    # Calculate properties
    df_with_props = calculate_pm7_dataframe_colab(df, smiles_column='smiles', charge_column='charge')
    
    print("\n\nResults DataFrame:")
    print(df_with_props[['name', 'charge', 'pm7_heat_of_formation', 'pm7_dipole_moment', 'pm7_success']])
    
    return df_with_props

# ============================================================================
# EXAMPLE 9: Comparing HOF for Neutral vs Charged
# ============================================================================

def example_hof_comparison():
    """Compare heat of formation for neutral molecule and its ions"""
    print("\n" + "="*70)
    print("EXAMPLE 9: HOF Comparison - Neutral vs Charged")
    print("="*70)
    
    molecule = "Methanol"
    smiles_neutral = "CO"
    smiles_cation = "C[OH2+]"     # Protonated methanol
    smiles_anion = "C[O-]"         # Methoxide anion
    
    calculator = ColabPM7Calculator()
    
    # Calculate all three
    print(f"\nCalculating properties for {molecule}...")
    props_neutral = calculator.calculate_properties(smiles_neutral, charge=0)
    props_cation = calculator.calculate_properties(smiles_cation, charge=1)
    props_anion = calculator.calculate_properties(smiles_anion, charge=-1)
    
    # Display comparison
    print("\n" + "="*70)
    print(f"HEAT OF FORMATION COMPARISON - {molecule}")
    print("="*70)
    
    species_data = [
        ("Neutral (CH3OH)", props_neutral),
        ("Cation (CH3OH2+)", props_cation),
        ("Anion (CH3O-)", props_anion),
    ]
    
    print(f"{'Species':<25} {'HOF (kcal/mol)':<20} {'Charge'}")
    print("-"*70)
    
    for name, props in species_data:
        if props['success']:
            hof = props.get('heat_of_formation', 'N/A')
            chg = props.get('input_charge', 'N/A')
            print(f"{name:<25} {hof:>15.2f}      {chg:>5}")
        else:
            print(f"{name:<25} FAILED")
    
    # Calculate ionization energy (approximate)
    if all(p['success'] for _, p in species_data):
        hof_neutral = props_neutral['heat_of_formation']
        hof_cation = props_cation['heat_of_formation']
        
        # Approximate ionization energy (not exact, needs electron affinity correction)
        print(f"\n📊 Energy differences:")
        print(f"   Protonation energy: {hof_cation - hof_neutral - 365.7:.2f} kcal/mol")
        print(f"   (PA would be negative of this)")
    
    return species_data

# ============================================================================
# EXAMPLE 10: Validation Against Known Values
# ============================================================================

def example_validation():
    """Validate PM7 results against known experimental values"""
    print("\n" + "="*70)
    print("EXAMPLE 10: Validation Against Experimental Data")
    print("="*70)
    
    # Known experimental proton affinities (kcal/mol)
    experimental_data = [
        ("Water", "O", "[OH3+]", 165.0),
        ("Ammonia", "N", "[NH4+]", 204.0),
        ("Methanol", "CO", "C[OH2+]", 180.3),
    ]
    
    calculator = ColabPM7Calculator()
    
    print("\nCalculating and comparing with experimental values...")
    print("\n" + "="*70)
    print(f"{'Molecule':<15} {'PM7 (kcal/mol)':<18} {'Exp (kcal/mol)':<18} {'Error'}")
    print("-"*70)
    
    errors = []
    
    for name, neutral, protonated, exp_pa in experimental_data:
        result = calculator.calculate_proton_affinity(neutral, protonated, cleanup=True)
        
        if result['success']:
            calc_pa = result['proton_affinity_kcal_mol']
            error = calc_pa - exp_pa
            errors.append(error)
            
            print(f"{name:<15} {calc_pa:>15.2f}   {exp_pa:>15.2f}   {error:>+8.2f}")
        else:
            print(f"{name:<15} FAILED")
    
    if errors:
        mean_error = sum(errors) / len(errors)
        mae = sum(abs(e) for e in errors) / len(errors)
        
        print("\n" + "="*70)
        print(f"Mean Error (ME):            {mean_error:>8.2f} kcal/mol")
        print(f"Mean Absolute Error (MAE):  {mae:>8.2f} kcal/mol")
        print(f"\n📊 PM7 performance: MAE typically 3-8 kcal/mol for PA")
        print("="*70)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║  PM7 Calculator - Usage Examples with Charge Support             ║
    ║  Version 0.2.0                                                   ║
    ╚══════════════════════════════════════════════════════════════════╝
    
    This script demonstrates all features of the enhanced PM7 calculator.
    Uncomment the examples you want to run below.
    """)
    
    # Run examples (uncomment as needed)
    
    # Basic calculations
    # example_neutral_molecule()
    # example_cation()
    # example_anion()
    # example_batch_mixed_charges()
    
    # Proton affinity calculations
    # example_proton_affinity_water()
    # example_proton_affinity_ammonia()
    # example_proton_affinity_organic_bases()
    
    # Advanced examples
    # example_dataframe_with_charges()
    # example_hof_comparison()
    # example_validation()
    
    print("\n✅ Examples completed! Uncomment specific examples to run them.")
