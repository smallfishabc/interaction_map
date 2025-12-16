# Examples

This directory contains example scripts demonstrating how to use the IDP Interaction Map package.

## 📋 Available Examples

### 1. Basic Usage (`basic_usage.py`)
**What it does**: Complete workflow for analyzing MD trajectories and generating interaction maps

**Use this example to learn**:
- How to load trajectory files
- How to generate contact maps
- How to visualize interaction networks
- Basic configuration options

**Run**:
```bash
python examples/basic_usage.py
```

### 2. Mutation Scanner (`generate_mutations.py`)
**What it does**: Generate targeted mutations based on interaction analysis

**Use this example to learn**:
- How to scan for mutation candidates
- How to generate single, pair, and chunk mutations
- How to filter by interaction strength
- How to protect critical regions (forbidden regions)
- How to save mutations for simulations

**Run**:
```bash
python examples/generate_mutations.py
```

### 3. Variant Comparison (`compare_variants.py`)
**What it does**: Compare interaction patterns between protein variants

**Use this example to learn**:
- How to analyze multiple variants
- How to compare interaction differences
- How to identify variant-specific interactions
- How to visualize comparative results

**Run**:
```bash
python examples/compare_variants.py
```

## 🚀 Quick Start

1. **Install the package** (if not already installed):
   ```bash
   pip install -e .
   ```

2. **Prepare your data**:
   - MD trajectory file (`.xtc`, `.dcd`, `.trr`, etc.)
   - Topology file (`.pdb`, `.gro`, `.psf`, etc.)
   - Protein sequence (one-letter amino acid codes)

3. **Modify an example**:
   - Copy one of the example scripts
   - Update file paths to your data
   - Adjust parameters as needed
   - Run the script

## 📚 Documentation

For more detailed information:
- **[Main README](../README.md)** - Package overview and installation
- **[Quick Reference](../QUICK_REFERENCE.md)** - Command cheat sheet
- **[Mutation Scanner Guide](../MUTATION_SCANNER.md)** - Detailed mutation scanner documentation
- **[API Documentation](../docs/README.md)** - Complete API reference

## 💡 Tips

### For Beginners
- Start with `basic_usage.py` to understand the workflow
- Read the inline comments in each script
- Try running with the default parameters first
- Experiment with different visualization options

### For Advanced Users
- Combine multiple analyses in one script
- Use the Python API for custom workflows
- Integrate with your existing analysis pipeline
- Create batch processing scripts

### Common Modifications

**Change input files**:
```python
# In any example script, update these paths
trajectory_file = "path/to/your/trajectory.xtc"
topology_file = "path/to/your/topology.pdb"
sequence = "YOUR_PROTEIN_SEQUENCE"
```

**Adjust visualization**:
```python
# Modify plot parameters
plot_interaction_network(
    df, 
    sequence, 
    threshold=0.05,  # Change interaction threshold
    output_file="custom_output.png"
)
```

**Change mutation parameters**:
```python
# In generate_mutations.py
results = scanner.generate_full_scan(
    output_dir=output_dir,
    interaction_type='attractive',  # or 'repulsive'
    min_chunk_strength=1.0,  # Adjust threshold
    generate_chunks=True,  # Enable/disable chunk mutations
    chunk_size=3  # Change chunk size
)
```

## 🐛 Troubleshooting

**Script fails to run?**
- Check that all dependencies are installed: `pip install -e .`
- Verify your input files exist and paths are correct
- Make sure your trajectory and topology are compatible

**Output looks wrong?**
- Check that your sequence matches your structure
- Verify you're using the correct mode (CG or all-atom)
- Try adjusting the `cutoff` and `threshold` parameters

**Need more help?**
- Read the [Quick Reference](../QUICK_REFERENCE.md)
- Check the [full documentation](../docs/README.md)
- Review the [Mutation Scanner Guide](../MUTATION_SCANNER.md)

## 🔬 Scientific Use

These examples demonstrate workflows suitable for:
- **Protein structure analysis**: Understanding IDP conformational preferences
- **Mutation design**: Engineering specific interaction properties
- **Comparative studies**: Analyzing differences between variants
- **High-throughput screening**: Batch analysis of multiple proteins

Remember to cite the appropriate methods and force fields used in your simulations!

## 📝 Creating Your Own Scripts

Feel free to use these examples as templates for your own analysis scripts. The basic pattern is:

1. **Import the package**:
   ```python
   from idp_interaction_map import analyze_contacts, plot_interaction_network
   ```

2. **Load and analyze data**:
   ```python
   df = analyze_contacts(
       trajectory_file="your_trajectory.xtc",
       topology_file="your_topology.pdb",
       sequence="YOURSEQUENCE"
   )
   ```

3. **Visualize or export results**:
   ```python
   plot_interaction_network(df, sequence)
   df.to_csv("results.csv")
   ```

4. **(Optional) Generate mutations**:
   ```python
   from idp_interaction_map import MutationScanner
   scanner = MutationScanner(df, sequence, "ProteinName")
   mutations = scanner.generate_full_scan(output_dir="mutations/")
   ```

Happy analyzing! 🎉
