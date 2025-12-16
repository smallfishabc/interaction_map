# Testing Guide for Scientists

## Introduction

This guide is designed for scientists who may not have extensive programming experience but want to understand, create, and modify tests for the IDP Interaction Map package. Tests help ensure that code works correctly and continues to work as changes are made.

## Table of Contents

1. [What Are Unit Tests?](#what-are-unit-tests)
2. [Why Testing Matters](#why-testing-matters)
3. [Setting Up Your Testing Environment](#setting-up-your-testing-environment)
4. [Understanding Test Structure](#understanding-test-structure)
5. [Writing Your First Test](#writing-your-first-test)
6. [Testing the Mutation Scanner](#testing-the-mutation-scanner)
7. [Common Testing Patterns](#common-testing-patterns)
8. [Running Tests](#running-tests)
9. [Debugging Failed Tests](#debugging-failed-tests)
10. [Best Practices](#best-practices)

---

## What Are Unit Tests?

**Unit tests** are small pieces of code that check if your main code works correctly. Think of them like quality control experiments in a lab:

- **Your Code** = Your experimental protocol
- **Unit Test** = A control experiment that verifies the protocol works
- **Test Suite** = A collection of all your control experiments

### Example Analogy

Imagine you have a function that calculates protein charge:

```python
def calculate_charge(sequence):
    """Calculate net charge of a protein sequence."""
    charge = 0
    for amino_acid in sequence:
        if amino_acid in ['K', 'R', 'H']:  # Positive
            charge += 1
        elif amino_acid in ['D', 'E']:  # Negative
            charge -= 1
    return charge
```

A unit test would verify this works:

```python
def test_calculate_charge():
    """Test that charge calculation is correct."""
    # Test a sequence with net positive charge
    result = calculate_charge("KRH")  # 3 positive residues
    assert result == 3, "Should have +3 charge"
    
    # Test a sequence with net negative charge
    result = calculate_charge("DE")  # 2 negative residues
    assert result == -2, "Should have -2 charge"
    
    # Test a neutral sequence
    result = calculate_charge("KDE")  # +1, -1, -1 = -1
    assert result == -1, "Should have -1 charge"
```

---

## Why Testing Matters

### For Scientists

1. **Reproducibility**: Tests document exactly what your code should do
2. **Confidence**: Know that your analysis code works correctly
3. **Catch Errors Early**: Find bugs before they affect your research
4. **Safe Refactoring**: Improve code without breaking functionality
5. **Documentation**: Tests show examples of how to use your code

### Real-World Example

Without tests:
```
You: "My mutation scanner gave weird results last month..."
Colleague: "Did you change anything?"
You: "I'm not sure... maybe?"
Result: Hours of debugging
```

With tests:
```
You: Make a change → Run tests → Tests fail
You: "Ah, I broke the charge calculation. Let me fix that."
Result: Bug caught in 2 minutes
```

---

## Setting Up Your Testing Environment

### 1. Install pytest

```bash
pip install pytest pytest-cov
```

### 2. Verify Installation

```bash
cd /path/to/interaction_map-0519_CG
pytest --version
```

You should see something like:
```
pytest 7.4.0
```

### 3. Run Existing Tests

```bash
pytest tests/
```

This will run all tests and show you the results.

---

## Understanding Test Structure

### Basic Test Anatomy

```python
# 1. Import the function you want to test
from idp_interaction_map.mutation_scanner import MutationScanner

# 2. Define a test function (must start with "test_")
def test_something():
    """
    3. Write a docstring explaining what you're testing
    """
    # 4. Set up test data (ARRANGE)
    sequence = "ACDEFG"
    
    # 5. Call the function being tested (ACT)
    result = do_something(sequence)
    
    # 6. Check the result is correct (ASSERT)
    assert result == expected_value, "Helpful error message"
```

### The AAA Pattern

Tests follow the **Arrange-Act-Assert** pattern:

```python
def test_filter_interactions():
    # ARRANGE: Set up test data
    interaction_data = pd.DataFrame({
        'r_1': [1, 2],
        'r_2': [5, 6],
        'cont_prob': [0.8, 0.05],  # One high, one low
        'distance': [4, 5],
        'relative_strength': [0.5, 0.6],
        'plot_value': [2, -1]
    })
    scanner = MutationScanner(interaction_data, "ACDEFGHIJ", "Test")
    
    # ACT: Call the function
    filtered = scanner.filter_interactions(min_contact_prob=0.1)
    
    # ASSERT: Check the results
    assert len(filtered) == 1, "Should filter out low probability interaction"
    assert filtered.iloc[0]['cont_prob'] == 0.8, "Should keep high probability"
```

---

## Writing Your First Test

Let's write a test step-by-step.

### Step 1: Choose What to Test

Let's test the `_calc_charge()` helper method from `MutationScanner`.

### Step 2: Create a Test File

Tests live in the `tests/` directory. If testing `mutation_scanner.py`, the test file is `tests/test_mutation_scanner.py`.

### Step 3: Import What You Need

```python
import pytest
import pandas as pd
from idp_interaction_map.mutation_scanner import MutationScanner
```

### Step 4: Write the Test

```python
def test_calc_charge_basic():
    """Test that charge calculation works for simple sequences."""
    # ARRANGE: Create a scanner object with minimal data
    test_data = pd.DataFrame({
        'r_1': [1], 
        'r_2': [5],
        'cont_prob': [0.5],
        'distance': [4],
        'relative_strength': [0.5],
        'plot_value': [1]
    })
    scanner = MutationScanner(test_data, "ACDEFG", "Test")
    
    # ACT & ASSERT: Test different charge scenarios
    
    # All positive residues: K, R, H
    result = scanner._calc_charge("KRH")
    assert result == 3, "KRH should have +3 charge"
    
    # All negative residues: D, E
    result = scanner._calc_charge("DE")
    assert result == -2, "DE should have -2 charge"
    
    # Mixed residues
    result = scanner._calc_charge("KRDE")
    assert result == 0, "KRDE should be neutral (2+ and 2-)"
    
    # Neutral residues: A, C, etc.
    result = scanner._calc_charge("ACF")
    assert result == 0, "ACF should have 0 charge"
```

### Step 5: Run Your Test

```bash
pytest tests/test_mutation_scanner.py::test_calc_charge_basic -v
```

Output:
```
tests/test_mutation_scanner.py::test_calc_charge_basic PASSED [100%]
```

✅ Success!

---

## Testing the Mutation Scanner

### Example 1: Testing Single Mutations

This tests whether the scanner can generate single-point mutations correctly.

```python
def test_generate_single_mutations_example():
    """
    Test generating single mutations for attractive interactions.
    
    What we're testing:
    - Can the scanner generate charge mutations?
    - Do we get the expected mutations for each position?
    - Are the sequences modified correctly?
    """
    # ARRANGE: Create test data
    interaction_data = pd.DataFrame({
        'r_1': [2, 3],           # Positions to mutate
        'r_2': [5, 6],           # Interaction partners
        'cont_prob': [0.8, 0.7],
        'distance': [4, 5],
        'relative_strength': [0.5, 0.6],
        'plot_value': [2, 1]
    })
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"  # 20 amino acids
    scanner = MutationScanner(interaction_data, sequence, "MyProtein")
    
    # Prepare candidate data
    candidates = pd.DataFrame({
        'r_1': [2],           # Position 2 (residue C)
        'r_2': [5],           # Position 5 (residue F)
        'strength': [2.0],
        'r_1_charge': [0],    # Neutral
        'r_2_charge': [-1]    # Negative
    })
    
    # ACT: Generate mutations
    mutations = scanner.generate_single_mutations(
        candidates, 
        position='left',        # Mutate the left position (r_1)
        mutation_type='charge', # Change charge
        interaction_type='attractive'  # Enhance attraction
    )
    
    # ASSERT: Check results
    assert len(mutations) > 0, "Should generate at least one mutation"
    assert 'mutation_name' in mutations.columns, "Should have mutation names"
    assert 'sequence' in mutations.columns, "Should have sequences"
    
    # Check that mutation names are formatted correctly
    first_mutation = mutations.iloc[0]
    assert 'MyProtein' in first_mutation['mutation_name'], "Should include protein name"
    assert 'C2' in first_mutation['mutation_name'], "Should show original residue and position"
    
    # Check that sequence was modified
    original_sequence = "ACDEFGHIKLMNPQRSTVWY"
    mutated_sequence = first_mutation['sequence']
    assert mutated_sequence != original_sequence, "Sequence should be changed"
    assert len(mutated_sequence) == len(original_sequence), "Length should be same"
```

### Example 2: Testing with Forbidden Regions

Scientists often want to protect certain regions (like active sites) from mutations.

```python
def test_forbidden_regions_example():
    """
    Test that forbidden regions are not mutated.
    
    Scenario: You have a protein with an active site at positions 5-10
    that you don't want to mutate.
    """
    # ARRANGE: Create interaction data
    interaction_data = pd.DataFrame({
        'r_1': [2, 5, 8],      # Position 5 and 8 are in forbidden region
        'r_2': [5, 8, 12],
        'cont_prob': [0.8, 0.7, 0.9],
        'distance': [4, 5, 6],
        'relative_strength': [0.5, 0.6, 0.7],
        'plot_value': [2, 1, 2]
    })
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    
    # Create scanner with forbidden region (positions 5-10 protected)
    scanner = MutationScanner(
        interaction_data, 
        sequence, 
        "MyProtein",
        forbidden_regions=[5, 6, 7, 8, 9, 10]  # Active site
    )
    
    # Prepare candidates
    candidates = pd.DataFrame({
        'r_1': [2, 5, 8],      # Try to mutate positions 2, 5, and 8
        'r_2': [5, 8, 12],
        'strength': [2.0, 1.5, 2.5],
        'r_1_charge': [0, 0, 0],
        'r_2_charge': [0, 0, 0]
    })
    
    # ACT: Generate mutations
    mutations = scanner.generate_single_mutations(
        candidates,
        position='left',
        mutation_type='charge',
        interaction_type='attractive'
    )
    
    # ASSERT: Check that forbidden positions were NOT mutated
    for _, mutation in mutations.iterrows():
        mutation_name = mutation['mutation_name']
        
        # Extract position from name (format: MyProtein_C2E)
        parts = mutation_name.split('_')[1:]  # Skip protein name
        for part in parts:
            # Position is the number in the middle (e.g., "C2E" → 2)
            position = int(''.join(filter(str.isdigit, part)))
            assert position not in [5, 6, 7, 8, 9, 10], \
                f"Position {position} should be protected!"
    
    # Position 2 should be mutated (not forbidden)
    assert len(mutations) > 0, "Should mutate non-forbidden positions"
```

### Example 3: Testing Edge Cases

Edge cases are unusual situations that might break your code.

```python
def test_boundary_conditions_example():
    """
    Test mutations at the edges of the sequence.
    
    What could go wrong:
    - Trying to create a chunk at position 1 might try to access position 0 or -1
    - Trying to mutate near the end might go past the sequence length
    """
    # ARRANGE: Small sequence
    interaction_data = pd.DataFrame({
        'r_1': [1, 10],        # First and last positions
        'r_2': [5, 5],
        'cont_prob': [0.8, 0.7],
        'distance': [4, 5],
        'relative_strength': [0.5, 0.6],
        'plot_value': [2, 1]
    })
    
    sequence = "ACDEFGHIJK"  # Only 10 residues
    scanner = MutationScanner(interaction_data, sequence, "Test")
    
    chunk_data = pd.DataFrame({
        'r_1': [1, 10],        # Edges
        'r_2': [5, 5],
        'strength': [2.0, 2.0],
        'r_1_charge': [1, 1],
        'r_2_charge': [-1, -1],
        'r_1_aromatic': [0, 0],
        'r_2_aromatic': [0, 0]
    })
    
    # ACT: Try to create chunk mutations at edges
    # Chunk size 3 at position 1 would need positions -1, 0, 1 (invalid!)
    mutations = scanner.generate_chunk_mutations(
        chunk_data,
        position='left',
        mutation_type='charge',
        interaction_type='attractive',
        chunk_size=3
    )
    
    # ASSERT: Should handle boundaries gracefully
    # No crashes = good!
    assert isinstance(mutations, pd.DataFrame), "Should return DataFrame"
    
    # Should skip invalid chunks
    for _, mutation in mutations.iterrows():
        seq = mutation['sequence']
        assert len(seq) == len(sequence), "Should maintain sequence length"
```

---

## Common Testing Patterns

### Pattern 1: Testing With Fixtures

Fixtures are reusable test data. They prevent code duplication.

```python
import pytest

@pytest.fixture
def sample_sequence():
    """Reusable test sequence."""
    return "ACDEFGHIKLMNPQRSTVWY"

@pytest.fixture
def sample_interaction_data():
    """Reusable interaction data."""
    return pd.DataFrame({
        'r_1': [1, 2, 3, 8],
        'r_2': [5, 6, 7, 12],
        'cont_prob': [0.8, 0.7, 0.9, 0.85],
        'distance': [4, 5, 6, 4],
        'relative_strength': [0.5, 0.6, 0.7, 0.65],
        'plot_value': [2, -1, 1, 2]
    })

@pytest.fixture
def scanner(sample_interaction_data, sample_sequence):
    """Reusable scanner object."""
    return MutationScanner(
        sample_interaction_data, 
        sample_sequence, 
        "TestProtein"
    )

# Now you can use these fixtures in any test
def test_with_fixture(scanner, sample_sequence):
    """This test automatically gets the fixtures."""
    assert scanner.sequence == sample_sequence
    assert scanner.protein_name == "TestProtein"
```

### Pattern 2: Testing Multiple Scenarios

Use `pytest.mark.parametrize` to test many scenarios efficiently.

```python
@pytest.mark.parametrize("sequence,expected_charge", [
    ("KKK", 3),      # Three positive
    ("EEE", -3),     # Three negative
    ("KE", 0),       # Neutral
    ("KRKR", 4),     # Four positive
    ("ACFG", 0),     # No charged residues
])
def test_charge_calculation(sequence, expected_charge):
    """Test charge calculation for various sequences."""
    scanner = create_test_scanner()  # Helper function
    result = scanner._calc_charge(sequence)
    assert result == expected_charge
```

### Pattern 3: Testing Error Handling

Make sure your code fails gracefully.

```python
def test_invalid_sequence():
    """Test that invalid sequences are caught."""
    interaction_data = create_test_data()
    
    # Try to create scanner with invalid sequence
    with pytest.raises(ValueError):
        scanner = MutationScanner(
            interaction_data,
            "XYZABC",  # Contains invalid amino acid X, Z
            "Test"
        )
```

### Pattern 4: Testing File I/O

Test functions that read or write files.

```python
import tempfile
from pathlib import Path

def test_save_mutations(scanner):
    """Test saving mutations to CSV file."""
    # Create mutations
    mutations = pd.DataFrame({
        'mutation_name': ['Test_A1E', 'Test_C2K'],
        'sequence': ['ECDEF', 'AKDEF']
    })
    
    # Use temporary directory (automatically cleaned up)
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "mutations.csv"
        
        # ACT: Save file
        scanner.save_mutations(mutations, output_path)
        
        # ASSERT: Check file exists and contents are correct
        assert output_path.exists(), "File should be created"
        
        loaded = pd.read_csv(output_path, header=None)
        assert len(loaded) == 2, "Should have 2 mutations"
        assert loaded.iloc[0, 0] == 'Test_A1E', "First mutation name correct"
```

---

## Running Tests

### Run All Tests

```bash
cd /path/to/interaction_map-0519_CG
pytest tests/
```

### Run Tests in One File

```bash
pytest tests/test_mutation_scanner.py
```

### Run One Specific Test

```bash
pytest tests/test_mutation_scanner.py::test_generate_single_mutations
```

### Run Tests with Verbose Output

```bash
pytest tests/ -v
```

Output shows each test:
```
tests/test_mutation_scanner.py::test_init PASSED                    [1%]
tests/test_mutation_scanner.py::test_filter_interactions PASSED     [2%]
...
```

### Run Tests with Coverage Report

```bash
pytest tests/ --cov=src/idp_interaction_map --cov-report=term-missing
```

Shows which lines of code are tested:
```
Name                                  Stmts   Miss  Cover   Missing
-------------------------------------------------------------------
src/idp_interaction_map/mutation_scanner.py   259      0   100%
```

### Run Tests Matching a Pattern

```bash
pytest tests/ -k "charge" -v
```

This runs only tests with "charge" in the name.

---

## Debugging Failed Tests

### Understanding Failure Messages

When a test fails, pytest shows detailed information:

```python
def test_example():
    result = calculate_charge("KKE")
    assert result == 5  # Wrong! Should be 1

# Output:
# >       assert result == 5
# E       assert 1 == 5
```

The `E` line shows what went wrong: the result was 1, but we expected 5.

### Adding Debug Information

```python
def test_with_debug_info():
    result = calculate_charge("KKE")
    print(f"DEBUG: result = {result}")  # Shows in output if test fails
    assert result == 1, f"Expected 1 but got {result}"
```

Run with `-s` to see print statements:

```bash
pytest tests/test_example.py::test_with_debug_info -s
```

### Using pytest's Built-in Debugger

```bash
pytest tests/ --pdb
```

When a test fails, drops into an interactive debugger where you can inspect variables.

### Common Failure Patterns

**1. Assertion Error**
```python
assert len(mutations) > 0
# AssertionError

# Fix: Add message
assert len(mutations) > 0, f"No mutations generated! Got: {mutations}"
```

**2. Attribute Error**
```python
result = scanner.nonexistent_method()
# AttributeError: 'MutationScanner' object has no attribute 'nonexistent_method'

# Fix: Check you're calling the right method
```

**3. Type Error**
```python
scanner = MutationScanner("ACDEF", interaction_data)  # Wrong order!
# TypeError: ...

# Fix: Check function signature
scanner = MutationScanner(interaction_data, "ACDEF", "Test")
```

---

## Best Practices

### 1. Test Names Should Be Descriptive

❌ Bad:
```python
def test_1():
    ...
```

✅ Good:
```python
def test_filter_interactions_removes_low_probability_contacts():
    ...
```

### 2. One Concept Per Test

❌ Bad: Testing everything in one test
```python
def test_everything():
    # Test filtering
    # Test mutations
    # Test saving
    # Test loading
    # ... 100 lines later
```

✅ Good: Separate tests
```python
def test_filter_interactions():
    ...

def test_generate_mutations():
    ...

def test_save_mutations():
    ...
```

### 3. Use Descriptive Assertions

❌ Bad:
```python
assert result == 5
```

✅ Good:
```python
assert result == 5, f"Expected 5 positive charges for 'KKKKK', got {result}"
```

### 4. Test Edge Cases

Always test:
- Empty inputs
- Boundary values (first/last positions)
- Invalid inputs
- Large inputs
- Typical inputs

```python
def test_charge_edge_cases():
    scanner = create_test_scanner()
    
    # Empty sequence
    assert scanner._calc_charge("") == 0
    
    # Single residue
    assert scanner._calc_charge("K") == 1
    
    # Very long sequence
    long_seq = "K" * 1000
    assert scanner._calc_charge(long_seq) == 1000
```

### 5. Keep Tests Independent

Each test should work on its own. Don't rely on test order.

❌ Bad:
```python
# test_a.py
global_scanner = None

def test_create_scanner():
    global global_scanner
    global_scanner = MutationScanner(...)

def test_use_scanner():
    # Relies on test_create_scanner running first!
    result = global_scanner.do_something()
```

✅ Good:
```python
def test_create_scanner():
    scanner = MutationScanner(...)
    assert scanner is not None

def test_use_scanner():
    scanner = MutationScanner(...)  # Create fresh scanner
    result = scanner.do_something()
```

### 6. Use Meaningful Test Data

❌ Bad: Mystery numbers
```python
data = pd.DataFrame({
    'r_1': [1, 2, 3],
    'r_2': [5, 6, 7],
    'cont_prob': [0.8, 0.7, 0.9],
    'distance': [4, 5, 6],
    'relative_strength': [0.5, 0.6, 0.7],
    'plot_value': [2, -1, 1]
})
```

✅ Good: Explained values
```python
# Create test data representing a small protein interaction network
# - Position 1 and 5: Strong attractive interaction (plot_value=2)
# - Position 2 and 6: Repulsive interaction (plot_value=-1)
# - Position 3 and 7: Weak attractive interaction (plot_value=1)
data = pd.DataFrame({
    'r_1': [1, 2, 3],           # First residue positions
    'r_2': [5, 6, 7],           # Partner residue positions
    'cont_prob': [0.8, 0.7, 0.9],  # High probability contacts
    'distance': [4, 5, 6],      # Sequence separation
    'relative_strength': [0.5, 0.6, 0.7],
    'plot_value': [2, -1, 1]    # Positive = attractive, negative = repulsive
})
```

---

## Quick Reference Card

### Creating a New Test

```python
# 1. Import
from idp_interaction_map.mutation_scanner import MutationScanner
import pandas as pd

# 2. Define test function
def test_my_feature():
    """What I'm testing"""
    
    # 3. ARRANGE: Set up
    data = create_test_data()
    scanner = MutationScanner(data, "ACDEF", "Test")
    
    # 4. ACT: Do something
    result = scanner.do_something()
    
    # 5. ASSERT: Check result
    assert result == expected_value, "Why this should be true"
```

### Running Tests

```bash
pytest tests/                              # All tests
pytest tests/test_file.py                  # One file
pytest tests/test_file.py::test_name       # One test
pytest tests/ -v                           # Verbose
pytest tests/ -k "pattern"                 # Match pattern
pytest tests/ --cov=src                    # With coverage
pytest tests/ -s                           # Show prints
pytest tests/ --pdb                        # Debug on failure
```

### Common Assertions

```python
assert x == y                              # Equality
assert x != y                              # Inequality
assert x > y                               # Greater than
assert x in [1, 2, 3]                      # Membership
assert len(x) == 5                         # Length
assert x is not None                       # Not None
assert isinstance(x, DataFrame)            # Type check
assert 'substring' in string               # Contains
```

---

## Getting Help

1. **Read the error message carefully** - It often tells you exactly what's wrong
2. **Check existing tests** - Look at `tests/test_mutation_scanner.py` for examples
3. **Use pytest documentation** - https://docs.pytest.org/
4. **Run tests with `-v`** - Verbose mode shows more information
5. **Add print statements** - Run with `-s` to see output
6. **Ask for help** - Share the error message and test code

---

## Next Steps

1. Read through the existing tests in `tests/test_mutation_scanner.py`
2. Try modifying a simple test and re-running it
3. Write a test for a new feature you're adding
4. Check the [Testing Tutorial](TESTING_TUTORIAL.md) for hands-on exercises
5. Review the [Test Examples](TESTING_EXAMPLES.md) for more patterns

Remember: **Writing tests is a skill that improves with practice!** Start simple and gradually build up to more complex tests.
