# Repository Structure

This document describes the organization of the IDP Interaction Map repository after the December 2025 cleanup.

## 📁 Directory Structure

```
interaction_map/
├── README.md                      # Main project documentation
├── LICENSE                        # MIT License
├── CHANGELOG.md                   # Version history
├── QUICK_REFERENCE.md             # Command cheat sheet
├── MUTATION_SCANNER.md            # Mutation scanner guide
├── pyproject.toml                 # Package configuration
├── setup.py                       # Installation script
│
├── src/                           # Source code
│   └── idp_interaction_map/       # Main package
│       ├── __init__.py
│       ├── core.py                # Core analysis functions
│       ├── contact_map.py         # Contact map generation
│       ├── normalization.py       # Normalization functions
│       ├── plotting.py            # Visualization
│       ├── mutation_scanner.py    # Mutation generation
│       ├── utils.py               # Utilities
│       ├── cli.py                 # Main CLI
│       └── cli_mutations.py       # Mutation CLI
│
├── tests/                         # Main test suite
│   ├── test_core.py               # Tests for core module
│   ├── test_contact_map.py        # Tests for contact maps
│   ├── test_normalization.py     # Tests for normalization
│   ├── test_plotting.py           # Tests for plotting
│   ├── test_mutation_scanner.py   # Tests for mutation scanner (fully annotated!)
│   └── test_utils.py              # Tests for utilities
│
├── examples/                      # Example scripts
│   ├── README.md                  # Examples guide
│   ├── basic_usage.py             # Basic workflow example
│   ├── generate_mutations.py     # Mutation scanner example
│   └── compare_variants.py       # Variant comparison example
│
├── docs/                          # Documentation
│   ├── README.md                  # Documentation index
│   └── testing/                   # Testing documentation
│       ├── README.md              # Testing docs index
│       ├── TESTING_GUIDE.md       # Complete testing guide (~2000 lines)
│       ├── TESTING_TUTORIAL.md    # Hands-on exercises
│       ├── TESTING_EXAMPLES.md    # Copy-paste templates
│       ├── TESTING_QUICK_REFERENCE.md  # Cheat sheet
│       └── TEST_ANNOTATIONS_SUMMARY.md # Code annotations guide
│
├── archive/                       # Historical documents
│   ├── README.md                  # Archive index
│   ├── MIGRATION.md               # Migration notes
│   ├── MODERNIZATION.md           # Modernization process
│   ├── IMPLEMENTATION_NOTES.md    # Technical notes
│   ├── PROJECT_COMPLETE.md        # Completion docs
│   ├── VALIDATION_RESULTS.md      # Validation data
│   └── ... (6 more historical docs)
│
├── dev-tests/                     # Development tests
│   ├── README.md                  # Dev tests guide
│   ├── test_all_atom_mode.py      # All-atom testing
│   ├── test_dual_mode_comparison.py  # Mode comparison
│   ├── test_multi_traj.py         # Multi-trajectory tests
│   ├── debug_traj.py              # Debug utilities
│   └── ... (test output directories)
│
└── legacy_code_backup/            # Original code (reference)
```

## 🎯 Quick Navigation

### For First-Time Users
1. Start with [README.md](README.md) - Installation and overview
2. Check [examples/README.md](examples/README.md) - Run example scripts
3. Read [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Command cheatsheet

### For Developers
1. Browse [tests/](tests/) - Main test suite
2. Read [docs/testing/TESTING_GUIDE.md](docs/testing/TESTING_GUIDE.md) - Learn testing
3. Check [src/idp_interaction_map/](src/idp_interaction_map/) - Source code

### For Contributors
1. Review [CHANGELOG.md](CHANGELOG.md) - Version history
2. Check [archive/](archive/) - Historical context
3. Read [dev-tests/README.md](dev-tests/README.md) - Development testing

## 📂 Directory Purposes

### Root Directory
**Purpose**: Essential files only - what users need first  
**Contents**: README, LICENSE, configuration files, main documentation  
**Why clean**: Makes first impression clear and professional

### `src/idp_interaction_map/`
**Purpose**: All source code  
**Contents**: Python modules implementing the package functionality  
**Unchanged**: Standard Python package structure

### `tests/`
**Purpose**: Automated test suite  
**Contents**: pytest tests for all modules  
**Note**: `test_mutation_scanner.py` has comprehensive line-by-line annotations

### `examples/`
**Purpose**: Example scripts demonstrating usage  
**Contents**: 3 example scripts + README explaining each  
**Usage**: Copy and modify for your own analyses

### `docs/`
**Purpose**: Comprehensive documentation  
**Contents**: Documentation index + testing subdirectory  
**Organization**: Main docs in root, testing docs in `testing/`

### `docs/testing/`
**Purpose**: Testing documentation for scientists  
**Contents**: 5 comprehensive guides (60+ pages total)  
**Audience**: Scientists new to testing  
**Highlight**: Complete learning path from basics to advanced

### `archive/`
**Purpose**: Historical development documents  
**Contents**: 10 documents covering migration, validation, project management  
**Usage**: Reference for understanding design decisions  
**Not needed**: For normal package usage

### `dev-tests/`
**Purpose**: Development and debug scripts  
**Contents**: Manual test scripts, test outputs, debug utilities  
**Usage**: Development and troubleshooting  
**Not needed**: For normal package usage (use `tests/` instead)

### `legacy_code_backup/`
**Purpose**: Original code before modernization  
**Contents**: Backup of pre-refactor code  
**Usage**: Historical reference only

## 🗂️ File Categories

### Essential Files (Always Read These)
- `README.md` - Start here!
- `QUICK_REFERENCE.md` - Commands
- `MUTATION_SCANNER.md` - Mutation guide
- `examples/README.md` - Examples guide
- `LICENSE` - Usage terms

### Documentation (Read When Needed)
- `docs/README.md` - Documentation index
- `docs/testing/*.md` - Testing guides (for writing tests)
- `CHANGELOG.md` - Version history

### Reference (Optional)
- `archive/*.md` - Historical documents
- `dev-tests/README.md` - Development testing info

### Configuration (Don't Usually Need to Modify)
- `pyproject.toml` - Package metadata and dependencies
- `setup.py` - Installation configuration

## 📊 Documentation Statistics

| Directory | Files | Total Lines | Purpose |
|-----------|-------|-------------|---------|
| Root | 5 docs | ~2,000 | Main guides |
| `docs/testing/` | 5 docs | ~10,000 | Testing education |
| `examples/` | 3 scripts + README | ~500 | Usage examples |
| `archive/` | 10 docs + README | ~5,000 | Historical reference |
| `dev-tests/` | 5 scripts + README | ~1,000 | Development tools |
| `tests/` | 6 test files | ~3,000 | Automated tests |

## 🎓 Learning Paths

### Path 1: Quick Start User
```
README.md
  ↓
examples/README.md
  ↓
examples/basic_usage.py
  ↓
QUICK_REFERENCE.md
```

### Path 2: Mutation Scanner User
```
README.md
  ↓
MUTATION_SCANNER.md
  ↓
examples/generate_mutations.py
  ↓
QUICK_REFERENCE.md
```

### Path 3: Developer/Contributor
```
README.md
  ↓
docs/README.md
  ↓
docs/testing/TESTING_GUIDE.md
  ↓
tests/test_mutation_scanner.py (annotated)
  ↓
src/idp_interaction_map/
```

### Path 4: Understanding History
```
README.md
  ↓
CHANGELOG.md
  ↓
archive/README.md
  ↓
archive/MIGRATION.md
```

## ✅ Design Principles

### 1. Progressive Disclosure
- Essential information in root
- Detailed docs in subdirectories
- Historical/development info in archives

### 2. Clear Purpose
- Each directory has a README explaining its purpose
- File names indicate content
- Related files grouped together

### 3. User-Focused
- First-time users see clean structure
- Common tasks easy to find
- Advanced features documented but not prominent

### 4. Maintainability
- Clear separation of concerns
- Easy to add new examples
- Development artifacts separated

## 🔄 Recent Changes (Dec 2025)

### Moved to `archive/`
- All historical project documents
- Validation results
- Migration notes
- Modernization documentation

### Moved to `dev-tests/`
- Development test scripts
- Debug utilities
- Test output directories
- Sample test data

### Moved to `docs/testing/`
- All testing guides
- Tutorial exercises
- Example templates
- Quick reference
- Test annotation summary

### Added READMEs
- `archive/README.md` - Explains archived documents
- `dev-tests/README.md` - Explains dev testing
- `docs/testing/README.md` - Testing docs index
- `examples/README.md` - Examples guide

### Updated Links
- Main README testing section
- Documentation index
- All cross-references

## 🎉 Result

**Before**: 14 files in root (cluttered, overwhelming)  
**After**: 7 files in root (clean, focused)

**Before**: Testing docs scattered  
**After**: All in `docs/testing/` with clear index

**Before**: Development artifacts mixed with user docs  
**After**: Separated into `archive/` and `dev-tests/`

**Result**: Professional, user-friendly repository structure! ✨
