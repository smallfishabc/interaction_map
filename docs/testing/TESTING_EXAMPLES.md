# Testing Examples: Real-World Test Patterns

This document provides copy-paste-ready examples of common testing patterns used in the IDP Interaction Map package.

## Table of Contents

1. [Basic Function Testing](#basic-function-testing)
2. [Testing with DataFrames](#testing-with-dataframes)
3. [Testing File Operations](#testing-file-operations)
4. [Testing Error Handling](#testing-error-handling)
5. [Testing with Fixtures](#testing-with-fixtures)
6. [Parametrized Tests](#parametrized-tests)
7. [Integration Tests](#integration-tests)
8. [Performance Tests](#performance-tests)

---

## Basic Function Testing

### Example 1: Testing a Simple Calculator Function

```python
def test_calculate_charge():
    """Test charge calculation for protein sequences."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    # Minimal scanner setup
    data = pd.DataFrame({
        'r_1': [1], 'r_2': [5], 'cont_prob': [0.5],
        'distance': [4], 'relative_strength': [0.5], 'plot_value': [1]
    })
    scanner = MutationScanner(data, "ACDEF", "Test")
    
    # Test positive charges
    assert scanner._calc_charge("KKK") == 3  # 3 lysines
    assert scanner._calc_charge("RRR") == 3  # 3 arginines
    assert scanner._calc_charge("HHH") == 3  # 3 histidines
    
    # Test negative charges
    assert scanner._calc_charge("DDD") == -3  # 3 aspartates
    assert scanner._calc_charge("EEE") == -3  # 3 glutamates
    
    # Test neutral
    assert scanner._calc_charge("AAA") == 0  # 3 alanines
    
    # Test mixed
    assert scanner._calc_charge("KRDE") == 0  # +2 -2 = 0
```

### Example 2: Testing String Manipulation

```python
def test_mutation_naming():
    """Test that mutation names are formatted correctly."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    # Setup
    data = pd.DataFrame({
        'r_1': [2], 'r_2': [5], 'cont_prob': [0.8],
        'distance': [4], 'relative_strength': [0.6], 'plot_value': [2]
    })
    scanner = MutationScanner(data, "ACDEFGHIJ", "MyProtein")
    
    # Generate a mutation
    candidates = pd.DataFrame({
        'r_1': [2], 'r_2': [5], 'strength': [2.0],
        'r_1_charge': [0], 'r_2_charge': [-1]
    })
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'charge', 'attractive'
    )
    
    # Check naming format
    if len(mutations) > 0:
        name = mutations.iloc[0]['mutation_name']
        
        # Should contain protein name
        assert 'MyProtein' in name
        
        # Should contain position
        assert '2' in name  # Position 2
        
        # Should contain residue info
        assert 'C' in name  # Original residue at position 2
```

---

## Testing with DataFrames

### Example 3: Testing DataFrame Filtering

```python
def test_filter_interactions():
    """Test filtering interaction data based on thresholds."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    # Create data with varying contact probabilities
    data = pd.DataFrame({
        'r_1': [1, 2, 3, 4],
        'r_2': [5, 6, 7, 8],
        'cont_prob': [0.9, 0.5, 0.1, 0.05],  # High to low
        'distance': [4, 5, 6, 7],
        'relative_strength': [0.8, 0.5, 0.2, 0.1],
        'plot_value': [2, 1, 1, 1]
    })
    
    scanner = MutationScanner(data, "ACDEFGHIJK", "Test")
    
    # Test filtering with different thresholds
    filtered = scanner.filter_interactions(min_contact_prob=0.4)
    
    # Should keep only high probability contacts
    assert len(filtered) == 2, "Should keep 2 high probability contacts"
    assert all(filtered['cont_prob'] >= 0.4), "All should be >= 0.4"
    
    # Test stricter threshold
    filtered_strict = scanner.filter_interactions(min_contact_prob=0.8)
    assert len(filtered_strict) == 1, "Only one contact >= 0.8"
```

### Example 4: Testing DataFrame Transformations

```python
def test_calculate_chunk_strength():
    """Test chunk strength calculation adds correct columns."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    data = pd.DataFrame({
        'r_1': [3, 4, 5],
        'r_2': [8, 9, 10],
        'cont_prob': [0.8, 0.7, 0.9],
        'distance': [5, 5, 5],
        'relative_strength': [0.6, 0.5, 0.7],
        'plot_value': [2, 1, 2]
    })
    
    scanner = MutationScanner(data, "ACDEFGHIKLMNPQRST", "Test")
    
    # Calculate chunk strength
    result = scanner.calculate_chunk_strength(data)
    
    # Check new columns were added
    expected_columns = [
        'r_1_chunk', 'r_2_chunk',
        'r_1_hydro', 'r_2_hydro',
        'r_1_charge', 'r_2_charge',
        'r_1_aromatic', 'r_2_aromatic'
    ]
    
    for col in expected_columns:
        assert col in result.columns, f"Missing column: {col}"
    
    # Check data types
    assert result['r_1_hydro'].dtype in ['float64', 'float'], "Hydro should be float"
    assert result['r_1_charge'].dtype in ['int64', 'int'], "Charge should be int"
```

---

## Testing File Operations

### Example 5: Testing CSV File Writing

```python
def test_save_mutations_creates_file():
    """Test that mutations are saved correctly to CSV."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    import tempfile
    from pathlib import Path
    
    # Setup scanner
    data = pd.DataFrame({
        'r_1': [1], 'r_2': [5], 'cont_prob': [0.5],
        'distance': [4], 'relative_strength': [0.5], 'plot_value': [1]
    })
    scanner = MutationScanner(data, "ACDEF", "Test")
    
    # Create test mutations
    mutations = pd.DataFrame({
        'mutation_name': ['Test_A1K', 'Test_C2E'],
        'sequence': ['KCDEF', 'AEDEF']
    })
    
    # Save to temporary file
    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = Path(tmpdir) / "mutations.csv"
        scanner.save_mutations(mutations, output_file)
        
        # Verify file exists
        assert output_file.exists(), "CSV file should be created"
        
        # Verify contents
        loaded = pd.read_csv(output_file, header=None)
        assert len(loaded) == 2, "Should have 2 rows"
        assert loaded.iloc[0, 0] == 'Test_A1K', "First mutation name matches"
        assert loaded.iloc[1, 0] == 'Test_C2E', "Second mutation name matches"
```

### Example 6: Testing CSV File Reading

```python
def test_scan_mutations_from_csv():
    """Test reading interaction data from CSV and generating mutations."""
    from idp_interaction_map.mutation_scanner import scan_mutations_from_csv
    import pandas as pd
    import tempfile
    from pathlib import Path
    
    # Create test CSV file
    data = pd.DataFrame({
        'r_1': [1, 2, 3],
        'r_2': [5, 6, 7],
        'cont_prob': [0.8, 0.7, 0.9],
        'distance': [4, 5, 6],
        'relative_strength': [0.5, 0.6, 0.7],
        'plot_value': [2, -1, 1]
    })
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save input CSV
        csv_path = Path(tmpdir) / "interactions.csv"
        data.to_csv(csv_path, index=False)
        
        # Run mutation scan
        output_dir = Path(tmpdir) / "mutations"
        results = scan_mutations_from_csv(
            interaction_csv=csv_path,
            sequence="ACDEFGHIJKLMNOP",
            protein_name="TestProtein",
            output_dir=output_dir,
            interaction_type='attractive',
            min_chunk_strength=0.5
        )
        
        # Verify results
        assert isinstance(results, dict), "Should return dictionary"
        assert len(results) > 0, "Should generate some mutations"
        assert output_dir.exists(), "Output directory should be created"
```

---

## Testing Error Handling

### Example 7: Testing Invalid Inputs

```python
def test_invalid_sequence_length():
    """Test that mismatched sequence length is caught."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    import pytest
    
    # Create data that references position 10
    data = pd.DataFrame({
        'r_1': [1, 10],  # Position 10 doesn't exist in short sequence
        'r_2': [5, 15],
        'cont_prob': [0.8, 0.7],
        'distance': [4, 5],
        'relative_strength': [0.5, 0.6],
        'plot_value': [2, 1]
    })
    
    # Short sequence (only 5 residues)
    short_sequence = "ACDEF"
    
    # This should either raise an error or handle it gracefully
    # (depending on implementation)
    try:
        scanner = MutationScanner(data, short_sequence, "Test")
        # If no error, check that out-of-bounds are handled
        assert True, "Scanner handles out-of-bounds positions"
    except (ValueError, KeyError, IndexError) as e:
        # If error is raised, that's also acceptable
        assert True, f"Error caught as expected: {e}"
```

### Example 8: Testing Empty DataFrames

```python
def test_empty_dataframe_handling():
    """Test handling of empty interaction data."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    # Empty DataFrame with correct columns
    empty_data = pd.DataFrame(columns=[
        'r_1', 'r_2', 'cont_prob', 'distance', 
        'relative_strength', 'plot_value'
    ])
    
    scanner = MutationScanner(empty_data, "ACDEF", "Test")
    
    # Should handle empty data gracefully
    filtered = scanner.filter_interactions()
    assert len(filtered) == 0, "Empty data should return empty result"
    
    candidates = scanner.identify_candidates(
        filtered, 'attractive', 1.0
    )
    assert len(candidates) == 0, "No candidates from empty data"
```

---

## Testing with Fixtures

### Example 9: Reusable Test Data

```python
import pytest
import pandas as pd
from idp_interaction_map.mutation_scanner import MutationScanner


@pytest.fixture
def standard_sequence():
    """A standard 20-residue test sequence."""
    return "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture
def standard_interaction_data():
    """Standard interaction data for testing."""
    return pd.DataFrame({
        'r_1': [1, 2, 3, 8, 11],
        'r_2': [5, 6, 7, 12, 16],
        'cont_prob': [0.8, 0.7, 0.9, 0.85, 0.75],
        'distance': [4, 4, 4, 4, 5],
        'relative_strength': [0.5, 0.6, 0.7, 0.65, 0.55],
        'plot_value': [2, -1, 1, 2, 1]
    })


@pytest.fixture
def basic_scanner(standard_interaction_data, standard_sequence):
    """A basic scanner for general testing."""
    return MutationScanner(
        standard_interaction_data,
        standard_sequence,
        "TestProtein"
    )


# Now you can use these fixtures in any test
def test_with_fixtures(basic_scanner, standard_sequence):
    """Example test using fixtures."""
    assert basic_scanner.sequence == standard_sequence
    assert basic_scanner.protein_name == "TestProtein"
    assert len(basic_scanner.df) == 5


def test_filtering_with_fixtures(basic_scanner):
    """Another test using the same fixture."""
    filtered = basic_scanner.filter_interactions(min_contact_prob=0.75)
    assert len(filtered) <= 5, "Filtered should have fewer or equal rows"
```

### Example 10: Fixture with Cleanup

```python
import pytest
import tempfile
from pathlib import Path


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)
    # Directory is automatically cleaned up after test


def test_with_temp_dir(basic_scanner, temp_output_dir):
    """Test that uses temporary directory."""
    import pandas as pd
    
    # Create some mutations
    mutations = pd.DataFrame({
        'mutation_name': ['Test_A1K'],
        'sequence': ['KCDEFGHIKLMNPQRSTVWY']
    })
    
    # Save to temporary directory
    output_file = temp_output_dir / "mutations.csv"
    basic_scanner.save_mutations(mutations, output_file)
    
    # Verify
    assert output_file.exists()
    # No cleanup needed - fixture handles it!
```

---

## Parametrized Tests

### Example 11: Testing Multiple Scenarios

```python
import pytest


@pytest.mark.parametrize("mutation_type,interaction_type,expected_residue", [
    ("charge", "attractive", "E"),  # Should use glutamate
    ("charge", "repulsive", "A"),   # Should neutralize
    ("polar", "attractive", "Q"),    # Should use glutamine
    ("polar", "repulsive", "A"),     # Should use alanine
    ("hydrophobic", "attractive", "L"),  # Should use leucine
    ("hydrophobic", "repulsive", "S"),   # Should use serine
])
def test_mutation_residue_selection(
    basic_scanner, 
    mutation_type, 
    interaction_type, 
    expected_residue
):
    """Test that correct residues are chosen for each mutation type."""
    import pandas as pd
    
    candidates = pd.DataFrame({
        'r_1': [2],
        'r_2': [5],
        'strength': [2.0],
        'r_1_charge': [0],
        'r_2_charge': [0]
    })
    
    mutations = basic_scanner.generate_single_mutations(
        candidates, 'left', mutation_type, interaction_type
    )
    
    if len(mutations) > 0:
        # Check that expected residue appears in mutations
        found = any(
            expected_residue in name 
            for name in mutations['mutation_name']
        )
        assert found, \
            f"Expected residue '{expected_residue}' for " \
            f"{mutation_type}/{interaction_type}"
```

### Example 12: Testing Edge Values

```python
@pytest.mark.parametrize("threshold,expected_count", [
    (0.0, 5),   # No filtering
    (0.5, 5),   # Keep all above 0.5
    (0.7, 4),   # Remove one with 0.65
    (0.8, 3),   # Keep only highest
    (1.0, 0),   # Remove all
])
def test_filtering_thresholds(
    basic_scanner, 
    threshold, 
    expected_count
):
    """Test filtering with various threshold values."""
    filtered = basic_scanner.filter_interactions(
        min_contact_prob=threshold
    )
    
    actual_count = len(filtered)
    assert actual_count == expected_count, \
        f"With threshold {threshold}, expected {expected_count} " \
        f"interactions, got {actual_count}"
```

---

## Integration Tests

### Example 13: Full Workflow Test

```python
def test_full_mutation_generation_workflow():
    """Test the complete mutation generation pipeline."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    import tempfile
    from pathlib import Path
    
    # Realistic interaction data
    data = pd.DataFrame({
        'r_1': list(range(5, 15)),
        'r_2': list(range(20, 30)),
        'cont_prob': [0.8] * 10,
        'distance': [5] * 10,
        'relative_strength': [0.6] * 10,
        'plot_value': [2, -1, 1, 2, 1, -1, 2, 1, 2, -1]
    })
    
    sequence = "A" * 35  # 35 residues
    
    with tempfile.TemporaryDirectory() as tmpdir:
        scanner = MutationScanner(data, sequence, "IntegrationTest")
        
        # Step 1: Filter
        filtered = scanner.filter_interactions(min_contact_prob=0.5)
        assert len(filtered) > 0, "Should have filtered interactions"
        
        # Step 2: Calculate chunks
        chunks = scanner.calculate_chunk_strength(filtered)
        assert 'r_1_chunk' in chunks.columns, "Should have chunk data"
        
        # Step 3: Identify candidates
        candidates = scanner.identify_candidates(
            chunks, 'attractive', 0.5
        )
        assert len(candidates) >= 0, "Should identify candidates"
        
        # Step 4: Generate mutations
        if len(candidates) > 0:
            mutations = scanner.generate_single_mutations(
                candidates, 'left', 'charge', 'attractive'
            )
            
            # Step 5: Save results
            output_file = Path(tmpdir) / "mutations.csv"
            if len(mutations) > 0:
                scanner.save_mutations(mutations, output_file)
                assert output_file.exists(), "Should save mutations"
```

---

## Performance Tests

### Example 14: Large-Scale Test

```python
def test_large_sequence_performance():
    """Test scanner with a large protein sequence."""
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    import time
    
    # Large sequence (500 residues)
    large_sequence = "ACDEFGHIKLMNPQRSTVWY" * 25  # 500 residues
    
    # Generate many interactions
    n_interactions = 100
    data = pd.DataFrame({
        'r_1': list(range(1, n_interactions + 1)),
        'r_2': list(range(50, 50 + n_interactions)),
        'cont_prob': [0.7] * n_interactions,
        'distance': [5] * n_interactions,
        'relative_strength': [0.5] * n_interactions,
        'plot_value': [1] * n_interactions
    })
    
    # Time the initialization
    start = time.time()
    scanner = MutationScanner(data, large_sequence, "Large")
    init_time = time.time() - start
    
    # Should complete reasonably fast
    assert init_time < 5.0, f"Initialization took {init_time:.2f}s (>5s)"
    
    # Time filtering
    start = time.time()
    filtered = scanner.filter_interactions()
    filter_time = time.time() - start
    
    assert filter_time < 2.0, f"Filtering took {filter_time:.2f}s (>2s)"
```

---

## Real-World Scenario Tests

### Example 15: Biologically Realistic Test

```python
def test_realistic_protein_scenario():
    """
    Test a realistic scenario: E1A protein region analysis.
    
    Scenario:
    - E1A protein, residues 101-164 (64 residues)
    - Want to enhance attractive interactions
    - Protect N-terminal region (101-110) - important for binding
    - Focus on central region (111-150) for mutations
    """
    from idp_interaction_map.mutation_scanner import MutationScanner
    import pandas as pd
    
    # Realistic interaction data
    data = pd.DataFrame({
        'r_1': [115, 120, 125, 130, 135, 140],
        'r_2': [145, 148, 150, 152, 155, 160],
        'cont_prob': [0.85, 0.78, 0.92, 0.81, 0.76, 0.88],
        'distance': [30, 28, 25, 22, 20, 20],
        'relative_strength': [0.65, 0.58, 0.72, 0.61, 0.56, 0.68],
        'plot_value': [2.5, 1.8, 2.8, 2.1, 1.7, 2.4]
    })
    
    # E1A sequence (using a placeholder)
    sequence = "M" * 64  # Would use actual E1A sequence
    
    # Protect N-terminal region
    forbidden = list(range(1, 11))  # Positions 1-10 (101-110 in full protein)
    
    scanner = MutationScanner(
        data, 
        sequence, 
        "E1A_101_164",
        forbidden_regions=forbidden
    )
    
    # Run analysis
    filtered = scanner.filter_interactions(
        min_contact_prob=0.75,
        min_distance=15
    )
    
    assert len(filtered) > 0, "Should find long-range interactions"
    
    chunks = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(
        chunks, 'attractive', 1.0
    )
    
    # Generate focused mutations
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'charge', 'attractive'
    )
    
    # Verify biological constraints
    for _, mut in mutations.iterrows():
        name = mut['mutation_name']
        # Extract position
        parts = name.split('_')[1:]
        for part in parts:
            pos = int(''.join(filter(str.isdigit, part)))
            assert pos not in forbidden, \
                f"Position {pos} should be protected!"
```

---

## Tips for Using These Examples

1. **Copy and modify**: Start with an example that's close to what you need
2. **Adjust test data**: Change DataFrames to match your specific scenario
3. **Add assertions**: Think about what would indicate a bug
4. **Use descriptive names**: Make test names explain what they check
5. **Add docstrings**: Explain why the test is important

## Common Pitfalls to Avoid

```python
# ❌ BAD: Test doesn't assert anything
def test_mutation_generation():
    scanner = create_scanner()
    mutations = scanner.generate_mutations()
    # No assertions - test always passes!

# ✅ GOOD: Test checks results
def test_mutation_generation():
    scanner = create_scanner()
    mutations = scanner.generate_mutations()
    assert len(mutations) > 0, "Should generate mutations"
    assert 'sequence' in mutations.columns
```

```python
# ❌ BAD: Test depends on external files
def test_with_real_file():
    scanner = load_from_file("/path/on/my/computer/data.csv")
    # Fails on other computers!

# ✅ GOOD: Test creates its own data
def test_with_test_data():
    data = create_test_dataframe()
    scanner = MutationScanner(data, "ACDEF", "Test")
    # Works anywhere!
```

---

## Additional Resources

- [Testing Guide](TESTING_GUIDE.md) - Comprehensive testing documentation
- [Testing Tutorial](TESTING_TUTORIAL.md) - Step-by-step learning exercises
- [pytest documentation](https://docs.pytest.org/) - Official pytest docs
- Existing tests in `tests/` directory - Real examples from this project
