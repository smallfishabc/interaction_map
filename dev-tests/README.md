# Development Tests

This directory contains development test scripts and test outputs that were used during development and debugging. These are not part of the main test suite.

## Contents

### Test Scripts
- **test_all_atom_mode.py** - Tests for all-atom mode functionality
- **test_dual_mode_comparison.py** - Comparing coarse-grained and all-atom modes
- **test_multi_traj.py** - Multi-trajectory processing tests
- **test_corrected_parameters.py** - Parameter validation tests
- **debug_traj.py** - Trajectory debugging utilities

### Test Output Directories
- **test_all_atom_output/** - Output from all-atom tests
- **test_cli_output/** - Output from CLI testing
- **test_output_multi_new/** - Multi-trajectory test outputs (new)
- **test_output_multi_old/** - Multi-trajectory test outputs (old)

### Test Data
- **test_0.6_contact_df_1201.csv** - Sample contact data for testing

## Purpose

These scripts were used for:
- Development and debugging
- Manual testing of features
- Validation of parameters
- Comparing different modes
- Generating test outputs

## Main Test Suite

**For running tests**, use the main test suite in the `tests/` directory:

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=src/idp_interaction_map

# Run specific test file
pytest tests/test_mutation_scanner.py
```

## For Developers

These development scripts may be useful for:
- Debugging specific issues
- Manual testing during development
- Comparing outputs between versions
- Reproducing specific test scenarios

**Note**: These are not maintained as part of the CI/CD pipeline. For automated testing, always use the `tests/` directory.
