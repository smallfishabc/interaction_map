# Testing Tutorial: Hands-On Exercises

This tutorial provides step-by-step exercises to help you learn testing by doing. Complete these exercises in order.

## Prerequisites

Make sure you have:
- Python environment set up
- pytest installed: `pip install pytest pytest-cov`
- The interaction_map package code

## Exercise 1: Your First Test

**Goal**: Create a simple test and run it successfully.

### Step 1: Create a Test File

Create a new file: `tests/test_tutorial.py`

```python
"""Tutorial tests for learning pytest."""
import pytest


def test_simple_math():
    """Test that basic math works."""
    result = 2 + 2
    assert result == 4, "2 + 2 should equal 4"
```

### Step 2: Run the Test

```bash
cd /path/to/interaction_map-0519_CG
pytest tests/test_tutorial.py -v
```

**Expected Output**:
```
tests/test_tutorial.py::test_simple_math PASSED                [100%]
```

✅ If you see PASSED, congratulations! You've run your first test.

### Step 3: Make It Fail

Change the assertion to see what a failure looks like:

```python
def test_simple_math():
    """Test that basic math works."""
    result = 2 + 2
    assert result == 5, "2 + 2 should equal 5"  # This is wrong!
```

Run again:
```bash
pytest tests/test_tutorial.py -v
```

**Expected Output**:
```
>       assert result == 5, "2 + 2 should equal 5"
E       AssertionError: 2 + 2 should equal 5
E       assert 4 == 5
```

**Learn**: The `E` lines show you what went wrong. Fix the test back to `assert result == 4`.

---

## Exercise 2: Testing a Real Function

**Goal**: Test one of the mutation scanner helper functions.

### Step 1: Import the Code

```python
"""Tutorial tests for mutation scanner."""
import pytest
import pandas as pd
from idp_interaction_map.mutation_scanner import MutationScanner


# Helper function to create a minimal scanner for testing
def create_test_scanner():
    """Create a scanner with minimal test data."""
    test_data = pd.DataFrame({
        'r_1': [1],
        'r_2': [5],
        'cont_prob': [0.5],
        'distance': [4],
        'relative_strength': [0.5],
        'plot_value': [1]
    })
    return MutationScanner(test_data, "ACDEFGHIJK", "Tutorial")
```

### Step 2: Test Hydrophobicity Calculation

The Kyte-Doolittle hydrophobicity scale assigns values to amino acids:
- Positive values = hydrophobic (water-repelling)
- Negative values = hydrophilic (water-loving)

```python
def test_hydrophobicity_calculation():
    """
    Test that hydrophobicity is calculated correctly.
    
    Hydrophobicity scale (Kyte-Doolittle):
    - Isoleucine (I) = 4.5 (very hydrophobic)
    - Aspartate (D) = -3.5 (very hydrophilic)
    - Alanine (A) = 1.8 (slightly hydrophobic)
    """
    scanner = create_test_scanner()
    
    # Test 1: Single hydrophobic residue
    result = scanner._calc_hydrophobicity("I")
    assert result == 4.5, "Isoleucine should be 4.5"
    
    # Test 2: Single hydrophilic residue
    result = scanner._calc_hydrophobicity("D")
    assert result == -3.5, "Aspartate should be -3.5"
    
    # Test 3: Multiple residues (sum them up)
    result = scanner._calc_hydrophobicity("III")
    assert result == 13.5, "Three I's should be 3 × 4.5 = 13.5"
    
    # Test 4: Mixed residues
    result = scanner._calc_hydrophobicity("IA")
    expected = 4.5 + 1.8  # I + A
    assert result == expected, f"I+A should be {expected}"
```

### Step 3: Run and Verify

```bash
pytest tests/test_tutorial.py::test_hydrophobicity_calculation -v
```

**Challenge**: Add one more test case with different amino acids. Look up their values in the `HYDROPHOBICITY` dictionary in `mutation_scanner.py`.

---

## Exercise 3: Testing with Different Inputs

**Goal**: Learn to test multiple scenarios using parametrize.

### The Scenario

We want to test charge calculation for many different sequences without writing separate tests for each.

```python
@pytest.mark.parametrize("sequence,expected_charge,description", [
    # (input_sequence, expected_result, what_we're_testing)
    ("KKK", 3, "three lysines"),
    ("EEE", -3, "three glutamates"),
    ("KE", 0, "one positive, one negative"),
    ("RRH", 3, "different positive residues"),
    ("DD", -2, "two aspartates"),
    ("ACF", 0, "neutral residues"),
    ("KRDE", 0, "balanced charges"),
])
def test_charge_calculation_comprehensive(sequence, expected_charge, description):
    """Test charge calculation for various sequences."""
    scanner = create_test_scanner()
    
    result = scanner._calc_charge(sequence)
    
    assert result == expected_charge, \
        f"Failed for {description}: '{sequence}' should have charge {expected_charge}, got {result}"
```

### Run This Test

```bash
pytest tests/test_tutorial.py::test_charge_calculation_comprehensive -v
```

You should see **7 tests** run (one for each parameter set):
```
test_charge_calculation_comprehensive[KKK-3-three lysines] PASSED
test_charge_calculation_comprehensive[EEE--3-three glutamates] PASSED
...
```

**Challenge**: Add 2 more test cases to the parametrize list.

---

## Exercise 4: Testing the Full Workflow

**Goal**: Test a complete mutation generation workflow.

### The Scenario

You want to test if the scanner can:
1. Take interaction data
2. Identify candidates for mutation
3. Generate single mutations
4. Return properly formatted results

```python
def test_single_mutation_workflow():
    """Test the complete workflow for generating single mutations."""
    
    # ARRANGE: Set up test data
    # This represents a protein with some interactions
    interaction_data = pd.DataFrame({
        'r_1': [2, 3, 4],           # Positions to analyze
        'r_2': [7, 8, 9],           # Their interaction partners
        'cont_prob': [0.8, 0.7, 0.9],
        'distance': [5, 5, 5],
        'relative_strength': [0.6, 0.5, 0.7],
        'plot_value': [2, 1, 2]     # Positive = attractive
    })
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    protein_name = "TutorialProtein"
    
    scanner = MutationScanner(interaction_data, sequence, protein_name)
    
    # ACT: Run the workflow
    # Step 1: Filter interactions
    filtered = scanner.filter_interactions(
        min_contact_prob=0.1,
        min_distance=4
    )
    
    # Step 2: Calculate chunk strengths
    chunk_data = scanner.calculate_chunk_strength(filtered)
    
    # Step 3: Identify candidates
    candidates = scanner.identify_candidates(
        chunk_data,
        interaction_type='attractive',
        min_chunk_strength=0.5
    )
    
    # Step 4: Generate mutations
    mutations = scanner.generate_single_mutations(
        candidates,
        position='left',
        mutation_type='charge',
        interaction_type='attractive'
    )
    
    # ASSERT: Check the results
    # 1. Should return a DataFrame
    assert isinstance(mutations, pd.DataFrame), \
        "Result should be a DataFrame"
    
    # 2. Should have the right columns
    assert 'mutation_name' in mutations.columns, \
        "Should have mutation_name column"
    assert 'sequence' in mutations.columns, \
        "Should have sequence column"
    
    # 3. Should have some mutations
    assert len(mutations) >= 0, \
        "Should return at least empty DataFrame"
    
    # 4. If mutations were generated, check their format
    if len(mutations) > 0:
        first_mutation = mutations.iloc[0]
        
        # Check mutation name format
        assert protein_name in first_mutation['mutation_name'], \
            f"Mutation name should include protein name '{protein_name}'"
        
        # Check sequence length unchanged
        assert len(first_mutation['sequence']) == len(sequence), \
            "Mutated sequence should have same length as original"
        
        # Check sequence was actually changed
        # (at least one mutation should change the sequence)
        sequences_changed = any(
            mut_seq != sequence 
            for mut_seq in mutations['sequence']
        )
        assert sequences_changed, \
            "At least one mutation should change the sequence"
    
    print(f"\n✓ Workflow complete! Generated {len(mutations)} mutations.")
```

### Run with Output

```bash
pytest tests/test_tutorial.py::test_single_mutation_workflow -v -s
```

The `-s` flag shows the print statement.

---

## Exercise 5: Testing Edge Cases

**Goal**: Learn to test unusual or extreme conditions.

### What Are Edge Cases?

Edge cases are situations that are at the "edges" of normal behavior:
- Empty inputs
- Very large inputs
- Boundary values (first/last items)
- Invalid inputs

```python
def test_edge_cases():
    """Test unusual but valid conditions."""
    scanner = create_test_scanner()
    
    # Edge Case 1: Empty sequence
    result = scanner._calc_charge("")
    assert result == 0, "Empty sequence should have 0 charge"
    
    # Edge Case 2: Single character
    result = scanner._calc_charge("K")
    assert result == 1, "Single K should have +1 charge"
    
    # Edge Case 3: Very long sequence
    long_sequence = "K" * 1000  # 1000 lysines
    result = scanner._calc_charge(long_sequence)
    assert result == 1000, "1000 K's should have +1000 charge"
    
    # Edge Case 4: All neutral residues
    result = scanner._calc_charge("AAAAA")
    assert result == 0, "All alanines should be neutral"
```

**Challenge**: Add a test for a sequence with ALL the positive residues (K, R, H) mixed together.

---

## Exercise 6: Testing with Temporary Files

**Goal**: Learn to test functions that read/write files safely.

### Why Temporary Files?

When testing file operations, we don't want to:
- Create files that clutter the project
- Risk overwriting important files
- Leave files behind if tests crash

Solution: Use Python's `tempfile` module.

```python
import tempfile
from pathlib import Path


def test_save_mutations_to_file():
    """Test that mutations can be saved to a CSV file."""
    
    # ARRANGE: Create scanner and mutations
    scanner = create_test_scanner()
    
    mutations = pd.DataFrame({
        'mutation_name': ['Tutorial_A1E', 'Tutorial_C2K', 'Tutorial_D3R'],
        'sequence': ['ECDEFGHIJK', 'AKDEFGHIJK', 'ACRЕФGHIJK']
    })
    
    # Use temporary directory (automatically cleaned up!)
    with tempfile.TemporaryDirectory() as tmpdir:
        # ACT: Save to file
        output_path = Path(tmpdir) / "test_mutations.csv"
        scanner.save_mutations(mutations, output_path)
        
        # ASSERT: Check file exists
        assert output_path.exists(), \
            f"File should be created at {output_path}"
        
        # Check file contents
        loaded_data = pd.read_csv(output_path, header=None)
        
        assert len(loaded_data) == 3, \
            f"Should have 3 rows, found {len(loaded_data)}"
        
        assert loaded_data.iloc[0, 0] == 'Tutorial_A1E', \
            "First mutation name should match"
        
        print(f"\n✓ File created and verified: {output_path}")
        print(f"✓ File will be auto-deleted when test completes")
    
    # File is automatically deleted when we exit the 'with' block
```

### Run and See

```bash
pytest tests/test_tutorial.py::test_save_mutations_to_file -v -s
```

---

## Exercise 7: Testing Forbidden Regions

**Goal**: Test that protected regions are respected.

### The Biological Context

Sometimes you have regions in a protein you don't want to mutate:
- Active sites
- Binding domains
- Structural motifs
- Post-translational modification sites

```python
def test_forbidden_regions_protection():
    """
    Test that forbidden regions are protected from mutations.
    
    Biological scenario: 
    - Protein has 20 residues
    - Positions 5-10 are an active site (don't mutate!)
    - Positions 15-17 are a binding site (don't mutate!)
    """
    
    # ARRANGE: Create interaction data
    interaction_data = pd.DataFrame({
        'r_1': [2, 5, 8, 15],       # Including forbidden positions
        'r_2': [12, 13, 14, 18],
        'cont_prob': [0.8, 0.7, 0.9, 0.85],
        'distance': [10, 8, 6, 3],
        'relative_strength': [0.6, 0.5, 0.7, 0.6],
        'plot_value': [2, 1, 2, 1]
    })
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    
    # Define forbidden regions
    active_site = [5, 6, 7, 8, 9, 10]
    binding_site = [15, 16, 17]
    forbidden = active_site + binding_site
    
    scanner = MutationScanner(
        interaction_data, 
        sequence, 
        "ProteinWithActiveSite",
        forbidden_regions=forbidden
    )
    
    # ACT: Try to generate mutations
    filtered = scanner.filter_interactions(min_contact_prob=0.1)
    chunk_data = scanner.calculate_chunk_strength(filtered)
    candidates = scanner.identify_candidates(
        chunk_data, 'attractive', 0.5
    )
    
    mutations = scanner.generate_single_mutations(
        candidates, 'left', 'charge', 'attractive'
    )
    
    # ASSERT: Check no forbidden positions were mutated
    for _, mutation in mutations.iterrows():
        mutation_name = mutation['mutation_name']
        
        # Extract position from mutation name
        # Format is: ProteinName_OriginalResidue Position MutatedResidue
        # Example: "ProteinWithActiveSite_C2K" means position 2
        parts = mutation_name.split('_')[1:]  # Skip protein name
        
        for part in parts:
            # Extract the number (position)
            position = int(''.join(filter(str.isdigit, part)))
            
            assert position not in forbidden, \
                f"Position {position} is forbidden but was mutated! " \
                f"Forbidden regions: active site {active_site}, " \
                f"binding site {binding_site}"
    
    print(f"\n✓ Tested {len(mutations)} mutations")
    print(f"✓ All mutations respect forbidden regions")
    print(f"✓ Active site (positions {active_site}) protected")
    print(f"✓ Binding site (positions {binding_site}) protected")
```

### Challenge

Modify the test to:
1. Add another forbidden region
2. Verify that at least one mutation WAS generated (not everything is forbidden)

---

## Exercise 8: Debugging a Failing Test

**Goal**: Learn to debug when tests don't pass.

### Create a Broken Test

```python
def test_intentionally_broken():
    """This test has a bug. Can you find it?"""
    scanner = create_test_scanner()
    
    # Calculate charge for lysines (K = positive)
    result = scanner._calc_charge("KKKK")
    
    # Bug: Expected value is wrong!
    assert result == 3, f"Four K's should give +3 charge, got {result}"
```

### Run It

```bash
pytest tests/test_tutorial.py::test_intentionally_broken -v
```

### Debug Steps

1. **Read the error message**:
   ```
   AssertionError: Four K's should give +3 charge, got 4
   assert 4 == 3
   ```

2. **The problem**: We expected 3, but got 4. Four K's = +4 charge (correct!)

3. **The fix**: Change assertion to `assert result == 4`

### Practice Debugging

Create this test and fix it:

```python
def test_debug_practice():
    """Find and fix the bug in this test."""
    scanner = create_test_scanner()
    
    # This should calculate hydrophobicity for "AAA"
    # Alanine (A) = 1.8 on Kyte-Doolittle scale
    result = scanner._calc_hydrophobicity("AAA")
    
    # Bug is here - can you spot it?
    assert result == 1.8, f"Expected 1.8, got {result}"
```

**Hint**: How many alanines are there?

---

## Exercise 9: Create Your Own Test

**Goal**: Write a test from scratch for a feature you care about.

### Choose a Function

Pick one function from `mutation_scanner.py` that you want to test. Examples:
- `filter_interactions()`
- `generate_pair_mutations()`
- `generate_chunk_mutations()`

### Write the Test

Follow the AAA pattern:

```python
def test_my_chosen_function():
    """
    Test [FUNCTION NAME] does [WHAT IT SHOULD DO].
    
    Why this is important:
    [EXPLAIN WHY THIS FUNCTION MATTERS FOR YOUR RESEARCH]
    """
    
    # ARRANGE: Set up test data
    # ... your setup code ...
    
    # ACT: Call the function
    # ... your function call ...
    
    # ASSERT: Check results
    # ... your assertions ...
    
    print("\n✓ My test passed!")
```

### Get Feedback

Run your test and make sure it:
1. Passes
2. Tests something meaningful
3. Would catch a real bug if the code broke

---

## Exercise 10: Test Coverage

**Goal**: Understand what code your tests cover.

### Run with Coverage

```bash
pytest tests/test_tutorial.py --cov=src/idp_interaction_map/mutation_scanner --cov-report=term-missing -v
```

### Understand the Report

```
Name                                    Stmts   Miss  Cover   Missing
---------------------------------------------------------------------
src/.../mutation_scanner.py               259     50    81%   45-52, 78-85
```

- **Stmts**: Total lines of code
- **Miss**: Lines not tested
- **Cover**: Percentage tested
- **Missing**: Which line numbers aren't tested

### Goal

Try to increase coverage by testing the "Missing" lines!

---

## Final Challenge: Comprehensive Test Suite

Create a test that exercises a complete mutation workflow:

```python
def test_complete_mutation_workflow():
    """
    Test a realistic mutation generation workflow.
    
    Scenario: You have a disordered protein with 30 residues.
    You want to:
    1. Identify strong attractive interactions
    2. Protect an important region (positions 10-15)
    3. Generate charge mutations to enhance attraction
    4. Save results to a file
    """
    
    # TODO: Implement this test!
    # 
    # Steps:
    # 1. Create realistic interaction data (at least 10 interactions)
    # 2. Create a scanner with forbidden regions
    # 3. Run the full workflow
    # 4. Verify results make biological sense
    # 5. Save to a temporary file
    # 6. Verify file contents
    
    pass  # Remove this and write your test!
```

---

## Checking Your Progress

### Beginner Level ✓

You should be able to:
- [ ] Run existing tests
- [ ] Understand test output
- [ ] Write a simple assertion
- [ ] Create a basic test function

### Intermediate Level ✓

You should be able to:
- [ ] Use pytest.mark.parametrize
- [ ] Test with temporary files
- [ ] Understand test fixtures
- [ ] Debug failing tests
- [ ] Check test coverage

### Advanced Level ✓

You should be able to:
- [ ] Design test scenarios
- [ ] Test edge cases
- [ ] Create reusable test fixtures
- [ ] Write comprehensive test suites
- [ ] Interpret coverage reports

---

## Next Steps

1. **Complete all exercises** in this tutorial
2. **Read existing tests** in `tests/test_mutation_scanner.py`
3. **Write tests for your own features**
4. **Review the Testing Guide** for more patterns
5. **Share your tests** with colleagues for feedback

## Getting Help

If you get stuck:
1. Read the error message carefully
2. Check the Testing Guide
3. Look at similar tests in the test suite
4. Run with `-v` for verbose output
5. Add print statements and run with `-s`
6. Ask a colleague or create an issue on GitHub

Remember: **Every expert was once a beginner!** Keep practicing and you'll get better at testing.
