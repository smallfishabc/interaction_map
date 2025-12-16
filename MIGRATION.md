# Migration Guide: Legacy → Modern IDP Interaction Map

This guide helps you migrate from the old codebase to the modernized package structure.

## Quick Start

### For First-Time Users

If you're new to the package, simply:

```bash
# Install
pip install -e .

# Run analysis
idp-interaction-map -d ./data -n my_protein -r 5
```

### For Existing Users

The old files remain in the root directory for backward compatibility. You can:

1. **Continue using old code** (no changes needed)
2. **Gradually migrate** to the new API
3. **Switch completely** to the modern version

## What Changed?

### File Mapping

| Old File | New Location | Purpose |
|----------|--------------|---------|
| `contact_map_generation.py` | `src/idp_interaction_map/contact_map.py` | Contact probability calculation |
| `normalization.py` | `src/idp_interaction_map/normalization.py` | Interaction strength normalization |
| `interaction_plot.py` | `src/idp_interaction_map/plotting.py` | Visualization |
| `default_function.py` | `src/idp_interaction_map/core.py` | Main workflow |
| `readpath.py` | `src/idp_interaction_map/utils.py` | File I/O utilities |
| `main.py` | `src/idp_interaction_map/cli.py` | Command-line interface |
| `showoff.py` | Integrated into `cli.py` | Welcome banner |

### Command-Line Usage

#### Old Way
```bash
python main.py --pdb protein.pdb --xtc trajectory.xtc -dir ./data -name my_protein
```

#### New Way
```bash
idp-interaction-map -p protein.pdb -x trajectory.xtc -d ./data -n my_protein
```

**Key Changes:**
- ✅ No need to call `python main.py`
- ✅ Cleaner flag names (`-d` instead of `-dir`)
- ✅ Consistent flag format (`-n` for name)
- ✅ Better help documentation (`--help`)

### Python API Usage

#### Old Way
```python
import os
import default_function
import readpath

# Change to data directory
os.chdir(path)

# Read sequence
seq = readpath.readsequence(path)

# Run analysis
default_function.interaction_map_pairwise(
    name, traj_path, seq, output_dir,
    pdb_top='__START_0.pdb',
    xtc_input=5,
    read_from_file=False
)
```

#### New Way
```python
from pathlib import Path
from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt

# Read sequence
sequence = read_sequence_from_txt(data_dir)

# Run analysis (no need to change directory!)
interaction_df = analyze_interaction_map(
    name="my_protein",
    trajectory_path=data_dir,
    sequence=sequence,
    output_dir=output_dir,
    pdb_top="__START_0.pdb",
    xtc_input=5,
    read_from_file=False
)

# Returns a DataFrame for further analysis
print(interaction_df.head())
```

**Key Improvements:**
- ✅ No directory changes needed
- ✅ Returns DataFrame for analysis
- ✅ Type hints for IDE support
- ✅ Better error messages
- ✅ Proper logging instead of prints

## Detailed Migration Examples

### Example 1: Single Trajectory Analysis

**Old Code:**
```python
import os
import readpath
import default_function

test = 1
path = 'F:/DATA/protein-summary'
name = path.split("\\")[-1].split("-")[0]
psi = '0'
residue = 'BB'

os.chdir(path)
seq = readpath.readsequence(path)
default_function.interaction_map_pairwise(
    name, path, seq, path,
    pdb_name='__START_0.pdb',
    xtc_name='__traj_0.xtc'
)
```

**New Code:**
```python
from pathlib import Path
from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt

# Define paths
data_dir = Path("./DATA/protein-summary")
name = data_dir.stem.split("-")[0]

# Read sequence
sequence = read_sequence_from_txt(data_dir)

# Run analysis
result = analyze_interaction_map(
    name=name,
    trajectory_path=data_dir,
    sequence=sequence,
    output_dir=data_dir,
    pdb_top="__START_0.pdb",
    xtc_input="__traj_0.xtc"
)
```

### Example 2: Multiple Trajectories

**Old Code:**
```python
def multi_traj_pre(name, path, psi, residue):
    _traj_path = os.path.join(path, residue, '_'.join(['S', psi]))
    _seq = readpath.readsequence(path)
    return _traj_path, _seq

traj_p, seq = multi_traj_pre(name, path, psi, residue)
default_function.interaction_map_pairwise(
    name, traj_p, seq, traj_p,
    xtc_input=5,
    read_from_file=False
)
```

**New Code:**
```python
from pathlib import Path
from idp_interaction_map import analyze_interaction_map
from idp_interaction_map.utils import read_sequence_from_txt

# Define paths
data_dir = Path("./data") / residue / f"S_{psi}"
sequence = read_sequence_from_txt(data_dir)

# Run analysis with 5 trajectories
result = analyze_interaction_map(
    name=name,
    trajectory_path=data_dir,
    sequence=sequence,
    output_dir=data_dir,
    xtc_input=5  # Loads __traj_0.xtc through __traj_4.xtc
)
```

### Example 3: Using Individual Components

If you need more control, use individual modules:

```python
from idp_interaction_map.contact_map import generate_contact
from idp_interaction_map.normalization import normalize_interaction_map
from idp_interaction_map.plotting import create_interaction_map

# Step 1: Generate contact map
contact_data = generate_contact(
    protein_name="my_protein",
    pdb_top="topology.pdb",
    xtc_input=5,
    cutoff=1.2
)

# Step 2: Normalize interactions
interaction_df = normalize_interaction_map(
    contact_data.contact,
    a1=13.12,
    b1=-2.32
)

# Step 3: Create visualization
create_interaction_map(
    seq=sequence,
    length=len(sequence),
    interaction_df=interaction_df,
    output_name="output_figure"
)
```

## Breaking Changes

### 1. Function Names
- `interaction_map_pairwise()` → `analyze_interaction_map()`
- `readsequence()` → `read_sequence_from_txt()` or `read_sequence_from_fasta()`

### 2. Return Values
- Old: Functions mostly had side effects (saved files)
- New: Functions return DataFrames for further analysis

### 3. Import Paths
- Old: `import default_function`
- New: `from idp_interaction_map import analyze_interaction_map`

### 4. Working Directory
- Old: Functions changed working directory with `os.chdir()`
- New: Functions use explicit paths (no directory changes)

## New Features

### 1. Logging
```python
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

from idp_interaction_map import analyze_interaction_map
# Will now see detailed progress messages
```

### 2. Type Checking
```python
# Modern IDEs will provide:
# - Autocomplete
# - Type hints
# - Error detection
from idp_interaction_map import analyze_interaction_map

# Your IDE knows the return type is pd.DataFrame!
result = analyze_interaction_map(...)
```

### 3. Better Error Messages
```python
# Old: Generic errors
# New: Helpful error messages with suggestions

# Example error:
# FileNotFoundError: Sequence file not found: ./data/seq.txt
# Please create seq.txt in the data directory with the protein sequence
```

## Testing Your Migration

### 1. Side-by-Side Comparison
Keep both implementations and compare outputs:

```python
# Run old version
old_result = old_function(...)

# Run new version
new_result = analyze_interaction_map(...)

# Compare results
assert np.allclose(old_result, new_result)
```

### 2. Unit Tests
The new package includes comprehensive tests:

```bash
pytest tests/ -v
```

### 3. Visual Inspection
Compare the generated plots to ensure consistency.

## Troubleshooting

### Issue: Import errors
```python
# Error: ModuleNotFoundError: No module named 'idp_interaction_map'

# Solution: Install the package
pip install -e .
```

### Issue: Can't find sequence file
```python
# Error: FileNotFoundError: Sequence file not found

# Solution: Ensure seq.txt exists in data directory
echo "MDEYK..." > data/seq.txt
```

### Issue: Old hardcoded paths
```python
# Old code had Windows paths like:
# path = 'F:\\DATA\\protein-summary'

# Update to cross-platform Path:
from pathlib import Path
path = Path("./DATA/protein-summary")
```

## Getting Help

1. **Documentation**: See `README.md` for full documentation
2. **Examples**: Check `examples/` directory for working code
3. **Tests**: Look at `tests/` to see usage patterns
4. **Issues**: The old code is still available if needed

## Gradual Migration Strategy

### Phase 1: Installation
```bash
pip install -e .
```

### Phase 2: Test CLI
```bash
idp-interaction-map --help
# Test on one protein
idp-interaction-map -d ./test_data -n test_protein -r 1
```

### Phase 3: Update Scripts
Replace one script at a time, keeping backups

### Phase 4: Validate Results
Compare outputs between old and new versions

### Phase 5: Full Transition
Once validated, remove old code dependencies

## Benefits of Migration

✅ **Better Performance**: Optimized algorithms  
✅ **Type Safety**: Catch errors before runtime  
✅ **IDE Support**: Autocomplete and hints  
✅ **Logging**: Better debugging  
✅ **Testing**: Comprehensive test suite  
✅ **Documentation**: Clear, complete docs  
✅ **Maintenance**: Modern code structure  
✅ **Extensibility**: Easy to add features  

## Still Need the Old Code?

The original files are preserved in the root directory. You can:

1. Continue using them as-is
2. Reference them during migration
3. Keep them as backup

No functionality has been removed—only improved and reorganized!
