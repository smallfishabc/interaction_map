# Testing Documentation

This directory contains comprehensive testing documentation designed specifically for scientists who may not have extensive programming experience.

## 📚 Documentation Files

### Getting Started
1. **[TESTING_GUIDE.md](TESTING_GUIDE.md)** - Comprehensive introduction (START HERE!)
   - What are unit tests and why they matter
   - Setting up your environment
   - Writing your first test
   - Best practices
   - ~60 pages, 30-60 minutes to read

2. **[TESTING_TUTORIAL.md](TESTING_TUTORIAL.md)** - Hands-on exercises
   - 10 progressive exercises
   - Real examples from the codebase
   - Practice problems with solutions
   - 2-4 hours with all exercises

### Reference Materials
3. **[TESTING_EXAMPLES.md](TESTING_EXAMPLES.md)** - Copy-paste-ready templates
   - 15+ common testing patterns
   - Real-world scenarios
   - Commented examples
   - Use when writing your own tests

4. **[TESTING_QUICK_REFERENCE.md](TESTING_QUICK_REFERENCE.md)** - One-page cheat sheet
   - Common commands
   - Common assertions
   - Quick troubleshooting
   - Print and keep handy!

5. **[TEST_ANNOTATIONS_SUMMARY.md](TEST_ANNOTATIONS_SUMMARY.md)** - Code walkthrough
   - Explanation of test file annotations
   - Line-by-line comment guide
   - Scientific context
   - How to read the annotated tests

## 🎯 Quick Start Paths

### Never Written a Test?
1. Read [TESTING_GUIDE.md](TESTING_GUIDE.md) sections 1-5
2. Do Exercise 1 from [TESTING_TUTORIAL.md](TESTING_TUTORIAL.md)
3. Keep [TESTING_QUICK_REFERENCE.md](TESTING_QUICK_REFERENCE.md) handy

### Need to Write Tests Now?
1. Check [TESTING_EXAMPLES.md](TESTING_EXAMPLES.md) for similar cases
2. Use [TESTING_QUICK_REFERENCE.md](TESTING_QUICK_REFERENCE.md) for syntax
3. Look at `../tests/test_mutation_scanner.py` for annotated examples

### Learning by Example?
1. Read [TEST_ANNOTATIONS_SUMMARY.md](TEST_ANNOTATIONS_SUMMARY.md)
2. Open `../tests/test_mutation_scanner.py`
3. Study the line-by-line explanations

## 🔗 Related Resources

- **Main Tests**: `../../tests/` - The actual test suite
- **Main Docs**: `../` - General documentation
- **Project Root**: `../../README.md` - Project overview

## 💡 Tips

- Start with the guide, practice with the tutorial
- Use examples as templates when writing new tests
- Keep the quick reference handy for common tasks
- Read the annotated test file to see concepts in action
- Tests help ensure your research is reproducible!

## 📊 Test Coverage

Current status:
- **Total Tests**: 83 (all passing)
- **Overall Coverage**: 82%
- **Mutation Scanner**: 100% ⭐

Run tests with:
```bash
pytest tests/
pytest tests/ --cov=src/idp_interaction_map
```
