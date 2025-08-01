#!/usr/bin/env python3
"""
Parameter Audit Script for panDecay
Comprehensive analysis of CLI parameters vs config file mappings vs analysis engine parameters
"""

import re
import argparse
from pathlib import Path

def extract_cli_parameters():
    """Extract all CLI parameters from main.py"""
    main_py_path = Path("src/main.py")
    
    with open(main_py_path, 'r') as f:
        content = f.read()
    
    # Find all add_argument calls
    arg_pattern = r'add_argument\(["\'](-{1,2}[^"\']+)["\'].*?\)'
    matches = re.findall(arg_pattern, content, re.MULTILINE | re.DOTALL)
    
    cli_params = []
    for match in matches:
        # Clean up parameter name
        param = match.replace('--', '').replace('-', '_')
        if param not in ['h', 'help']:  # Skip help parameters
            cli_params.append(param)
    
    # Add positional argument
    if 'alignment' not in cli_params:
        cli_params.append('alignment')
    
    return sorted(set(cli_params))

def extract_config_parameters():
    """Extract all config file parameters from configuration.py"""
    config_py_path = Path("src/core/configuration.py")
    
    with open(config_py_path, 'r') as f:
        content = f.read()
    
    # Find the param_map dictionary
    param_map_match = re.search(r'param_map = \{(.*?)\}', content, re.MULTILINE | re.DOTALL)
    
    if not param_map_match:
        return []
    
    param_map_content = param_map_match.group(1)
    
    # Extract config parameter names (keys in the dictionary)
    config_params = re.findall(r"'([^']+)':\s*'[^']+'", param_map_content)
    
    return sorted(set(config_params))

def extract_analysis_engine_parameters():
    """Extract all analysis engine __init__ parameters"""
    analysis_engine_path = Path("src/core/analysis_engine.py")
    
    with open(analysis_engine_path, 'r') as f:
        content = f.read()
    
    # Find the __init__ method parameters
    init_match = re.search(r'def __init__\(self,([^)]+)\):', content, re.MULTILINE | re.DOTALL)
    
    if not init_match:
        return []
    
    params_content = init_match.group(1)
    
    # Extract parameter names (before = or before comma/newline)
    params = re.findall(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*(?:=|,|\n)', params_content)
    
    # Clean up and filter
    engine_params = []
    for param in params:
        param = param.strip()
        if param and param not in ['self']:
            engine_params.append(param)
    
    return sorted(set(engine_params))

def extract_config_template_parameters():
    """Extract all parameters from config template"""
    template_path = Path("src/core/config_template.txt")
    
    with open(template_path, 'r') as f:
        content = f.read()
    
    # Find all parameter assignments (param = value)
    param_pattern = r'^([a-zA-Z_][a-zA-Z0-9_]*)\s*='
    matches = re.findall(param_pattern, content, re.MULTILINE)
    
    return sorted(set(matches))

def normalize_parameter_name(param):
    """Normalize parameter name (handle dashes vs underscores)"""
    return param.replace('-', '_')

def main():
    print("=" * 80)
    print("PANDECAY PARAMETER AUDIT")
    print("=" * 80)
    print()
    
    # Extract parameters from all sources
    cli_params = extract_cli_parameters()
    config_params = extract_config_parameters()
    engine_params = extract_analysis_engine_parameters()
    template_params = extract_config_template_parameters()
    
    # Normalize all parameter names
    cli_params = [normalize_parameter_name(p) for p in cli_params]
    config_params = [normalize_parameter_name(p) for p in config_params]
    engine_params = [normalize_parameter_name(p) for p in engine_params]
    template_params = [normalize_parameter_name(p) for p in template_params]
    
    print(f"CLI Parameters Found: {len(cli_params)}")
    print(f"Config Map Parameters Found: {len(config_params)}")
    print(f"Analysis Engine Parameters Found: {len(engine_params)}")
    print(f"Config Template Parameters Found: {len(template_params)}")
    print()
    
    # Find missing mappings
    cli_set = set(cli_params)
    config_set = set(config_params)
    engine_set = set(engine_params)
    template_set = set(template_params)
    
    # Parameters in CLI but not in config mapping
    missing_in_config = cli_set - config_set
    
    # Parameters in CLI but not in template
    missing_in_template = cli_set - template_set
    
    # Parameters in config but not in CLI (shouldn't happen)
    extra_in_config = config_set - cli_set
    
    # Parameters in CLI/config but not in engine
    missing_in_engine = cli_set - engine_set
    
    print("MISSING CONFIG FILE MAPPINGS:")
    print("=" * 40)
    if missing_in_config:
        for param in sorted(missing_in_config):
            print(f"  ❌ --{param.replace('_', '-')} (missing from config parser)")
    else:
        print("  ✅ All CLI parameters have config mappings")
    print()
    
    print("MISSING CONFIG TEMPLATE ENTRIES:")
    print("=" * 40)
    if missing_in_template:
        for param in sorted(missing_in_template):
            print(f"  ❌ {param} (missing from config template)")
    else:
        print("  ✅ All CLI parameters have template entries")
    print()
    
    print("EXTRA CONFIG MAPPINGS:")
    print("=" * 40)
    if extra_in_config:
        for param in sorted(extra_in_config):
            print(f"  ⚠️  {param} (in config but not CLI)")
    else:
        print("  ✅ No extra config mappings")
    print()
    
    print("MISSING ANALYSIS ENGINE PARAMETERS:")
    print("=" * 40)
    if missing_in_engine:
        for param in sorted(missing_in_engine):
            print(f"  ❌ {param} (CLI parameter not passed to engine)")
    else:
        print("  ✅ All CLI parameters reach analysis engine")
    print()
    
    # Detailed parameter lists
    print("DETAILED PARAMETER LISTS:")
    print("=" * 40)
    
    print(f"\nCLI Parameters ({len(cli_params)}):")
    for i, param in enumerate(cli_params, 1):
        status = "✅" if param in config_set else "❌"
        print(f"  {i:2d}. {status} {param}")
    
    print(f"\nConfig Mapping Parameters ({len(config_params)}):")
    for i, param in enumerate(config_params, 1):
        print(f"  {i:2d}. {param}")
    
    print(f"\nAnalysis Engine Parameters ({len(engine_params)}):")
    for i, param in enumerate(engine_params, 1):
        print(f"  {i:2d}. {param}")
    
    # Summary
    print("\n" + "=" * 80)
    print("AUDIT SUMMARY:")
    print("=" * 80)
    
    total_issues = len(missing_in_config) + len(missing_in_template) + len(missing_in_engine)
    
    if total_issues == 0:
        print("✅ PERFECT CONSISTENCY - All parameters properly mapped across all interfaces")
    else:
        print(f"❌ ISSUES FOUND: {total_issues} parameter mapping problems detected")
        print("\nACTION ITEMS:")
        if missing_in_config:
            print(f"  1. Add {len(missing_in_config)} missing config mappings to configuration.py")
        if missing_in_template:
            print(f"  2. Add {len(missing_in_template)} missing parameters to config_template.txt")
        if missing_in_engine:
            print(f"  3. Review {len(missing_in_engine)} parameters not reaching analysis engine")

if __name__ == "__main__":
    main()