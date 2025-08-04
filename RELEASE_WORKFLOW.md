# panDecay Release Workflow

This document describes the process for releasing new versions of panDecay to PyPI.

## Prerequisites

### One-time Setup

1. **PyPI Account**: Create account at https://pypi.org/account/register/
2. **TestPyPI Account**: Create account at https://test.pypi.org/account/register/
3. **API Tokens**: Generate API tokens for both PyPI and TestPyPI
4. **Tools Installation**: 
   ```bash
   pipx install build
   pipx install twine
   ```

### Configure Credentials

Create or update `~/.pypirc`:
```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-YOUR_PYPI_TOKEN_HERE

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-YOUR_TESTPYPI_TOKEN_HERE
```

## Release Process

### 1. Prepare Release

1. **Update Version**: Edit `pandecay/core/constants.py`
   ```python
   VERSION = "1.2.0"  # New version number
   ```

2. **Update Version in pyproject.toml**:
   ```toml
   [project]
   version = "1.2.0"
   ```

3. **Update CHANGELOG.md**: Add new release section with changes

4. **Test Locally**: Ensure all tests pass and functionality works

### 2. Build Package

```bash
# Clean previous builds
rm -rf dist/ build/ *.egg-info/

# Build package
pyproject-build
```

This creates:
- `dist/pandecay-X.Y.Z.tar.gz` (source distribution)
- `dist/pandecay-X.Y.Z-py3-none-any.whl` (wheel)

### 3. Test on TestPyPI

```bash
# Upload to TestPyPI
twine upload --repository testpypi dist/*

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pandecay

# Test the installed package
pandecay --help
```

### 4. Create Git Tag and Push

```bash
# Commit version changes
git add .
git commit -m "Release version 1.2.0"

# Create and push tag
git tag v1.2.0
git push origin main
git push origin v1.2.0
```

### 5. Upload to PyPI

```bash
# Upload to real PyPI
twine upload dist/*
```

### 6. Verify Release

```bash
# Test installation from PyPI
pip install --upgrade pandecay

# Verify version
pandecay --version
```

### 7. Post-Release

1. **GitHub Release**: Create release on GitHub with changelog
2. **Update Development Version**: Bump to next dev version (e.g., 1.2.1-dev)

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):

- **Major (X.0.0)**: Breaking changes
- **Minor (1.X.0)**: New features, backward compatible
- **Patch (1.1.X)**: Bug fixes, backward compatible

## Common Issues

### Build Failures
- Check `pyproject.toml` syntax
- Ensure all imports use absolute paths
- Verify `MANIFEST.in` includes necessary files

### Upload Failures
- Check API token permissions
- Verify version number is unique (can't overwrite)
- Ensure package name availability

### Import Errors
- Test package installation in clean environment
- Check entry point configuration
- Verify dependencies are properly specified

## Automation (Future)

Consider setting up GitHub Actions for automatic releases:

1. **Trigger**: Push to `v*` tags
2. **Build**: Automatic package building
3. **Test**: Run test suite
4. **Upload**: Automatic PyPI upload
5. **Release**: Create GitHub release

## Emergency Procedures

### Yanking a Release
If a critical issue is discovered:

```bash
# Yank the problematic version (keeps it installed but hidden)
twine upload --repository pypi --skip-existing --verbose dist/*
# Then use PyPI web interface to yank the version
```

### Hot Fix Release
For critical bugs:

1. Create hotfix branch from tag
2. Apply minimal fix
3. Release as patch version (e.g., 1.1.1)
4. Merge back to main

## Contact

For questions about the release process:
- GitHub Issues: https://github.com/mol-evol/panDecay/issues
- Maintainer: James McInerney