# Documentation Index

Welcome to the IDP Interaction Map documentation!

## Getting Started

- **[README](../README.md)** - Main project overview, installation, and quick start
- **[Mutation Scanner Guide](MUTATION_SCANNER.md)** - Complete guide to the mutation scanner module

## Testing Documentation (For Scientists)

Our testing documentation is designed specifically for scientists who may not have extensive programming experience. Start here if you're new to testing:

### 📚 Learning Path

1. **[Testing Guide](TESTING_GUIDE.md)** - START HERE!
   - What are unit tests and why they matter
   - Setting up your testing environment
   - Understanding test structure
   - Writing your first test
   - Common testing patterns
   - Running and debugging tests
   - Best practices
   - **Audience**: Complete beginners to testing
   - **Time**: 30-60 minutes to read, longer to practice

2. **[Testing Tutorial](TESTING_TUTORIAL.md)** - Hands-on practice
   - 10 step-by-step exercises
   - Real examples from the codebase
   - Progressive difficulty
   - Practice problems with solutions
   - **Audience**: After reading the Testing Guide
   - **Time**: 2-4 hours with all exercises

3. **[Testing Examples](TESTING_EXAMPLES.md)** - Reference library
   - 15+ copy-paste-ready test patterns
   - Real-world scenarios
   - Organized by test type
   - Commented examples
   - **Audience**: When writing your own tests
   - **Time**: Use as needed (reference)

4. **[Testing Quick Reference](TESTING_QUICK_REFERENCE.md)** - Cheat sheet
   - One-page reference card (printable!)
   - Common commands
   - Common assertions
   - Quick troubleshooting
   - **Audience**: Everyone (keep handy!)
   - **Time**: 5 minutes to review

### 🎯 Choose Your Path

**I've never written a test before:**
1. Read [Testing Guide](TESTING_GUIDE.md) (focus on sections 1-5)
2. Do Exercise 1 from [Testing Tutorial](TESTING_TUTORIAL.md)
3. Keep [Quick Reference](TESTING_QUICK_REFERENCE.md) open
4. Practice with Exercise 2-3

**I understand the basics:**
1. Skim [Testing Guide](TESTING_GUIDE.md) sections 6-9
2. Do [Testing Tutorial](TESTING_TUTORIAL.md) exercises 4-7
3. Use [Testing Examples](TESTING_EXAMPLES.md) as templates

**I need to write tests now:**
1. Check [Testing Examples](TESTING_EXAMPLES.md) for similar cases
2. Use [Quick Reference](TESTING_QUICK_REFERENCE.md) for syntax
3. Refer to existing tests in `tests/` directory
4. Consult [Testing Guide](TESTING_GUIDE.md) for patterns

**I'm debugging a test:**
1. Check [Quick Reference](TESTING_QUICK_REFERENCE.md) "Debugging" section
2. Review [Testing Guide](TESTING_GUIDE.md) "Debugging Failed Tests"
3. Look at similar tests in `tests/` directory

## Feature Documentation

### Core Features

- **Contact Map Generation** - Generate contact maps from MD trajectories
- **Interaction Analysis** - Analyze residue-residue interactions
- **Visualization** - Create publication-quality figures
- **Dual Mode Support** - Both coarse-grained and all-atom simulations

### Mutation Scanner

- **[Mutation Scanner Guide](MUTATION_SCANNER.md)** - Complete documentation
  - Overview and installation
  - Quick start examples
  - Detailed usage guide
  - CLI and Python API
  - Mutation types and strategies
  - Forbidden regions
  - Output format
  - Advanced features
  - Troubleshooting

## API Reference

The package is organized into the following modules:

### Main Modules

- `idp_interaction_map.core` - Core functionality for contact analysis
- `idp_interaction_map.contact_map` - Contact map generation
- `idp_interaction_map.normalization` - Normalization functions
- `idp_interaction_map.plotting` - Visualization functions
- `idp_interaction_map.mutation_scanner` - Mutation generation (NEW!)
- `idp_interaction_map.utils` - Utility functions

### Command-Line Tools

- `idp-interaction-map` - Main CLI tool for interaction analysis
- `idp-mutation-scan` - CLI tool for mutation generation (NEW!)

## Development

### For Contributors

- **Testing Philosophy** - We emphasize testing for research reproducibility
- **Code Quality** - Black formatting, isort for imports, ruff for linting
- **Coverage Goals** - Core modules at 100%, overall >80%

### Running Tests

```bash
# All tests
pytest tests/

# Specific module
pytest tests/test_mutation_scanner.py

# With coverage
pytest tests/ --cov=src/idp_interaction_map
```

### Current Test Status

- **Total Tests**: 83 (all passing)
- **Overall Coverage**: 82%
- **Mutation Scanner Coverage**: 100% ⭐

## Additional Resources

### External Links

- **GitHub Repository**: https://github.com/smallfishabc/interaction_map
- **pytest Documentation**: https://docs.pytest.org/
- **MDTraj Documentation**: http://mdtraj.org/

### Examples

The `examples/` directory contains:
- `generate_mutations.py` - Mutation scanner examples
- Additional usage examples (coming soon)

## FAQ

### Testing

**Q: I've never written tests before. Where do I start?**  
A: Start with the [Testing Guide](TESTING_GUIDE.md), then do Exercise 1 in the [Tutorial](TESTING_TUTORIAL.md).

**Q: How do I run tests?**  
A: `pytest tests/` from the project root directory. See [Quick Reference](TESTING_QUICK_REFERENCE.md) for more commands.

**Q: My test is failing. How do I debug it?**  
A: Check the "Debugging Failed Tests" section in the [Testing Guide](TESTING_GUIDE.md) or [Quick Reference](TESTING_QUICK_REFERENCE.md).

**Q: How do I test file I/O?**  
A: Use `tempfile.TemporaryDirectory()` - see examples in [Testing Examples](TESTING_EXAMPLES.md) section "Testing File Operations".

### Mutation Scanner

**Q: How do I generate mutations?**  
A: See [Mutation Scanner Guide](MUTATION_SCANNER.md) Quick Start section.

**Q: What's the difference between attractive and repulsive mutations?**  
A: Attractive mutations enhance favorable interactions, repulsive mutations disrupt unfavorable ones. Details in [Mutation Scanner Guide](MUTATION_SCANNER.md).

**Q: How do I protect important regions?**  
A: Use the `forbidden_regions` parameter. See [Mutation Scanner Guide](MUTATION_SCANNER.md) "Forbidden Regions" section.

### General

**Q: What Python version do I need?**  
A: Python 3.8 or higher.

**Q: How do I report a bug?**  
A: Create an issue on GitHub with error messages and minimal example.

**Q: Can I contribute?**  
A: Yes! Read the testing documentation, write tests for your changes, and submit a pull request.

## Getting Help

1. **Check this documentation** - Start with relevant guide
2. **Read error messages** - They often tell you what's wrong
3. **Look at examples** - Check `tests/` and `examples/` directories
4. **GitHub Issues** - Search or create an issue
5. **Ask colleagues** - Share code and error messages

## Document Overview

| Document | Purpose | Length | Best For |
|----------|---------|--------|----------|
| [Testing Guide](TESTING_GUIDE.md) | Comprehensive introduction | ~60 min read | Learning testing from scratch |
| [Testing Tutorial](TESTING_TUTORIAL.md) | Hands-on exercises | 2-4 hours | Practice and skill building |
| [Testing Examples](TESTING_EXAMPLES.md) | Code templates | Reference | Writing your own tests |
| [Quick Reference](TESTING_QUICK_REFERENCE.md) | Cheat sheet | 5 min | Quick lookups |
| [Mutation Scanner](MUTATION_SCANNER.md) | Feature guide | 30 min read | Using mutation scanner |

## Updates

- **2024-12-15**: Added comprehensive testing documentation for scientists
- **2024-12-15**: Achieved 100% test coverage for mutation scanner
- **2024-12-15**: Added mutation scanner module
- **2024-12**: Initial documentation release

---

**Need help?** Start with the [Testing Guide](TESTING_GUIDE.md) if you're new to testing, or jump to the [Quick Reference](TESTING_QUICK_REFERENCE.md) if you need a quick syntax reminder!
