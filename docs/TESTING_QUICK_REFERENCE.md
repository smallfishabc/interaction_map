# Testing Quick Reference Card

**Print this page and keep it handy while writing tests!**

---

## Basic Test Template

```python
def test_my_feature():
    """Describe what you're testing."""
    # ARRANGE: Set up test data
    data = create_test_data()
    
    # ACT: Do something
    result = do_something(data)
    
    # ASSERT: Check result
    assert result == expected, "Helpful error message"
```

---

## Common Commands

```bash
# Run all tests
pytest tests/

# Run one file
pytest tests/test_mutation_scanner.py

# Run one test
pytest tests/test_file.py::test_name

# Verbose output
pytest tests/ -v

# Show print statements
pytest tests/ -s

# Check coverage
pytest tests/ --cov=src/idp_interaction_map

# Debug on failure
pytest tests/ --pdb
```

---

## Common Assertions

```python
assert x == y               # Equal
assert x != y               # Not equal
assert x > y                # Greater than
assert x < y                # Less than
assert x >= y               # Greater or equal
assert x <= y               # Less or equal
assert x in [1, 2, 3]       # In list
assert x not in [1, 2]      # Not in list
assert len(x) == 5          # Length
assert x is None            # Is None
assert x is not None        # Is not None
assert isinstance(x, str)   # Type check
assert 'sub' in string      # Contains substring
assert x                    # Truthy
assert not x                # Falsy
```

---

## Creating Test Data

```python
import pandas as pd

# Minimal interaction data
data = pd.DataFrame({
    'r_1': [1, 2, 3],
    'r_2': [5, 6, 7],
    'cont_prob': [0.8, 0.7, 0.9],
    'distance': [4, 5, 6],
    'relative_strength': [0.5, 0.6, 0.7],
    'plot_value': [2, -1, 1]
})

# Create scanner
from idp_interaction_map.mutation_scanner import MutationScanner
scanner = MutationScanner(data, "ACDEFGHIJ", "Test")
```

---

## Testing with Temporary Files

```python
import tempfile
from pathlib import Path

def test_file_operation():
    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = Path(tmpdir) / "test.csv"
        # Do file operations
        save_data(output_file)
        # Check results
        assert output_file.exists()
    # File auto-deleted here!
```

---

## Using Fixtures (Reusable Data)

```python
import pytest

@pytest.fixture
def test_sequence():
    return "ACDEFGHIKLMNPQRSTVWY"

@pytest.fixture  
def test_scanner(test_sequence):
    data = create_test_data()
    return MutationScanner(data, test_sequence, "Test")

# Use in test
def test_something(test_scanner):
    result = test_scanner.do_something()
    assert result is not None
```

---

## Parametrized Tests (Multiple Scenarios)

```python
@pytest.mark.parametrize("input,expected", [
    ("KKK", 3),      # 3 positive
    ("EEE", -3),     # 3 negative
    ("KE", 0),       # Neutral
])
def test_charge(input, expected):
    scanner = create_scanner()
    assert scanner._calc_charge(input) == expected
```

---

## Testing Errors

```python
import pytest

def test_invalid_input():
    with pytest.raises(ValueError):
        do_something_invalid()
```

---

## Checking Test Coverage

```bash
# Run with coverage report
pytest tests/ --cov=src/idp_interaction_map --cov-report=term-missing

# Output shows:
# Name                    Stmts   Miss  Cover   Missing
# -----------------------------------------------------
# module.py                 100     10    90%   45-52
#                          ^^^^   ^^^^  ^^^^    ^^^^^
#                         total  missed  %    line numbers
```

---

## Debugging Failed Tests

1. **Read the error message**
   ```
   >       assert result == 5
   E       assert 3 == 5
   ```
   
2. **Add debug info**
   ```python
   print(f"DEBUG: result = {result}")
   assert result == expected, f"Got {result}, expected {expected}"
   ```

3. **Run with -s to see prints**
   ```bash
   pytest tests/test_file.py::test_name -s
   ```

4. **Use debugger**
   ```bash
   pytest tests/ --pdb
   ```

---

## AAA Pattern (Arrange-Act-Assert)

```python
def test_example():
    # ARRANGE: Set up
    scanner = create_scanner()
    data = prepare_data()
    
    # ACT: Do the thing
    result = scanner.process(data)
    
    # ASSERT: Check it worked
    assert result.is_valid()
    assert len(result) > 0
```

---

## Best Practices

✅ **DO:**
- Write descriptive test names
- Test one thing per test
- Use meaningful test data
- Add helpful assertion messages
- Keep tests independent
- Test edge cases (empty, boundary values)

❌ **DON'T:**
- Test multiple unrelated things
- Use mysterious magic numbers
- Depend on test order
- Skip assertion messages
- Leave tests commented out
- Hardcode file paths

---

## Common Test Patterns

### Testing a calculation
```python
def test_calculation():
    result = calculate_charge("KKE")  # 2 pos, 1 neg
    assert result == 1, f"Expected +1, got {result}"
```

### Testing data filtering
```python
def test_filter():
    data = create_data_with_high_and_low_values()
    filtered = filter_data(data, threshold=0.5)
    assert all(filtered['value'] >= 0.5)
```

### Testing file I/O
```python
def test_save_load():
    with tempfile.TemporaryDirectory() as tmpdir:
        file_path = Path(tmpdir) / "data.csv"
        save_data(data, file_path)
        loaded = load_data(file_path)
        assert len(loaded) == len(data)
```

### Testing with forbidden regions
```python
def test_forbidden():
    scanner = MutationScanner(
        data, sequence, "Test",
        forbidden_regions=[5, 6, 7]
    )
    mutations = generate_mutations(scanner)
    for mut in mutations:
        pos = extract_position(mut['mutation_name'])
        assert pos not in [5, 6, 7]
```

---

## Getting Help

1. Read error messages carefully
2. Check [Testing Guide](TESTING_GUIDE.md)
3. Look at existing tests in `tests/`
4. Use `-v` for verbose output
5. Add print statements and use `-s`
6. Ask a colleague!

---

## Test File Structure

```
tests/
├── __init__.py
├── test_mutation_scanner.py    # Mutation scanner tests
├── test_core.py                # Core functionality tests  
├── test_plotting.py            # Visualization tests
└── test_utils.py               # Utility function tests
```

---

## Quick Diagnostic

**Test won't run?**
- Check it starts with `test_`
- Check file starts with `test_`
- Check imports are correct

**Test always passes?**
- Add an assertion!
- Check assertion is reachable

**Test fails unexpectedly?**
- Read full error message
- Check test data is correct
- Use `print()` and `-s` flag
- Verify function signature

---

**Remember**: Tests are **experiments** that verify your code works correctly!

**For detailed explanations and tutorials, see:**
- [Testing Guide](TESTING_GUIDE.md) - Complete introduction
- [Testing Tutorial](TESTING_TUTORIAL.md) - Hands-on exercises  
- [Testing Examples](TESTING_EXAMPLES.md) - Copy-paste examples
