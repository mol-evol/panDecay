# panDecay Code Style Guide

## Overview

This document outlines the code style standards for the panDecay project. The goal is to maintain readable, consistent, and maintainable code following Python best practices.

## Code Style Standards

### 1. Line Length
- **Maximum line length**: 100 characters
- **Preferred line length**: 88 characters (Black formatter default)
- **Long lines should be broken using**:
  - Parentheses for expressions
  - Backslashes for simple continuations
  - Triple quotes for multi-line strings

### 2. String Formatting
- **Prefer f-strings** for string formatting
- **Break long f-strings** into multiple lines:
  ```python
  # Good
  message = (f"Processing {filename} with {count} sequences "
            f"using {model} model")
  
  # Avoid
  message = f"Processing {filename} with {count} sequences using {model} model"
  ```

### 3. Function Signatures
- **Break long function signatures** after the opening parenthesis:
  ```python
  # Good
  def process_alignment(self, alignment_file: str, 
                       output_dir: Path, model: str = "GTR",
                       threads: int = 1) -> Dict[str, Any]:
  
  # Avoid
  def process_alignment(self, alignment_file: str, output_dir: Path, model: str = "GTR", threads: int = 1) -> Dict[str, Any]:
  ```

### 4. Complex Expressions
- **Break complex expressions** using parentheses:
  ```python
  # Good
  result = (some_long_variable_name + another_long_variable_name *
           third_variable_name / fourth_variable_name)
  
  # Avoid
  result = some_long_variable_name + another_long_variable_name * third_variable_name / fourth_variable_name
  ```

### 5. Logging Statements
- **Break long logging statements**:
  ```python
  # Good
  logger.info(f"Analysis completed: {analysis_type} with {count} sequences "
             f"in {duration:.2f} seconds")
  
  # Avoid
  logger.info(f"Analysis completed: {analysis_type} with {count} sequences in {duration:.2f} seconds")
  ```

## Current Code Quality Metrics

### Line Length Analysis
- **Total lines**: 7,617
- **Lines > 100 characters**: 509 (6.7%)
- **Lines > 120 characters**: ~200 (2.6%)

### Common Long Line Patterns
1. **Table formatting** (35% of long lines)
2. **Logging statements** (25% of long lines)
3. **Complex expressions** (20% of long lines)
4. **Function signatures** (15% of long lines)
5. **String operations** (5% of long lines)

## Improvement Priorities

### High Priority (Affecting readability)
1. **Function signatures** - Break long parameter lists
2. **Logging statements** - Split informational messages
3. **Complex expressions** - Use parentheses for clarity

### Medium Priority (Affecting maintainability)
1. **Error messages** - Split long error descriptions
2. **Configuration parsing** - Break long conditional statements
3. **Mathematical operations** - Clarify complex calculations

### Low Priority (Cosmetic)
1. **Table formatting** - ASCII art tables (acceptable to exceed 100 chars)
2. **Constants and headers** - Documentation and metadata
3. **Import statements** - Already well-formatted

## Automated Tools

### Recommended Tools
- **Black**: Python code formatter (88 character line length)
- **isort**: Import sorting and organization
- **flake8**: Style guide enforcement
- **mypy**: Type checking

### Configuration
```ini
# setup.cfg
[tool:isort]
line_length = 88
multi_line_output = 3

[flake8]
max-line-length = 100
extend-ignore = E203, W503
```

## Implementation Strategy

### Phase 1: Critical Lines (Completed)
- [x] Function signatures with type hints
- [x] Long logging statements  
- [x] Complex mathematical expressions
- [x] Error message formatting

### Phase 2: Systematic Improvements
- [ ] Break remaining long function signatures
- [ ] Split complex conditional statements
- [ ] Improve string formatting consistency
- [ ] Add docstring line breaks

### Phase 3: Automated Formatting
- [ ] Apply Black formatter with 88-character limit
- [ ] Run isort for import organization
- [ ] Validate with flake8

## Exception Cases

### Acceptable Long Lines
1. **URLs and file paths** - Keep on single line for clarity
2. **ASCII art and tables** - Preserve visual formatting
3. **Regular expressions** - Keep pattern on single line
4. **Documentation strings** - Some flexibility for readability

### Example Exceptions
```python
# Acceptable: URL
HELP_URL = "https://github.com/user/repo/wiki/very-long-documentation-page-name"

# Acceptable: Table formatting
f.write(f"│ {name:<20} │ {value:>10} │ {status:^15} │")

# Acceptable: Regex pattern
pattern = r"^(\w+)\s*=\s*([+-]?\d*\.?\d+([eE][+-]?\d+)?)\s*$"
```

## Code Review Guidelines

### Before Committing
1. Check line length with `python -c "max(len(line) for line in open('file.py'))"`
2. Verify readability of broken lines
3. Ensure logical grouping of related parameters
4. Test that code still functions correctly

### Review Checklist
- [ ] No lines exceed 100 characters unnecessarily
- [ ] Long lines are broken at logical points
- [ ] Broken lines are properly indented
- [ ] Complex expressions are clarified with parentheses
- [ ] String formatting is consistent