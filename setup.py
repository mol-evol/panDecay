#!/usr/bin/env python3
"""
Setup script for panDecay package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# Read requirements from requirements.txt
requirements = []
requirements_path = this_directory / "requirements.txt"
if requirements_path.exists():
    requirements = requirements_path.read_text().strip().split('\n')
    requirements = [req.strip() for req in requirements if req.strip() and not req.startswith('#')]

# Read version from package
import pandecay
version = pandecay.__version__

setup(
    name="pandecay",
    version=version,
    author="James McInerney",
    author_email="mcinerney.james+panGPT@gmail.com",
    description="Phylogenetic Analysis using Decay Indices",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mol-evol/panDecay",
    project_urls={
        "Bug Reports": "https://github.com/mol-evol/panDecay/issues",
        "Source": "https://github.com/mol-evol/panDecay",
        "Documentation": "https://github.com/mol-evol/panDecay#readme",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="phylogenetics bioinformatics decay-indices bremer-support maximum-likelihood bayesian parsimony",
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "visualization": ["matplotlib>=3.5.0", "seaborn>=0.11.0"],
        "dev": ["pytest>=6.0", "pytest-cov", "black", "flake8"],
    },
    entry_points={
        "console_scripts": [
            "pandecay=pandecay.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "pandecay": [
            "examples/data/*",
            "examples/*",
            "images/*",
        ],
    },
    zip_safe=False,
)