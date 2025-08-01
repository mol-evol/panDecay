#!/usr/bin/env python3
"""
Parameter Consistency Test Suite for panDecay
Tests CLI parameters, config file parameters, and runtime display consistency
"""

import tempfile
import argparse
from pathlib import Path
import sys
import os

# Add src to path to import modules
sys.path.insert(0, str(Path(__file__).parent / "src"))

# Import using absolute imports to avoid relative import issues
import importlib.util

def load_module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

# Load the required modules
src_path = Path(__file__).parent / "src"
main_module = load_module_from_path("main", src_path / "main.py")
utils_module = load_module_from_path("utils", src_path / "core" / "utils.py")
config_module = load_module_from_path("config", src_path / "core" / "configuration.py")

# Import the functions we need
build_effective_model_display = utils_module.build_effective_model_display
parse_config = config_module.parse_config
generate_config_template = config_module.generate_config_template
setup_argument_parser = main_module.setup_argument_parser

def test_model_parameter_consistency():
    """Test that model parameter overrides work correctly across CLI and config"""
    
    print("=" * 80)
    print("MODEL PARAMETER CONSISTENCY TESTS")
    print("=" * 80)
    
    test_cases = [
        {
            "name": "DNA NST=1 Override (Original Issue)",
            "cli_args": ["--nst", "1", "--base-freq", "empirical"],
            "expected_display": "JC (nst=1, basefreq=empirical)",
            "config_content": "nst = 1\nbase_freq = empirical\n"
        },
        {
            "name": "DNA NST=2 with Gamma", 
            "cli_args": ["--nst", "2", "--gamma", "--gamma-shape", "0.5"],
            "expected_display": "HKY+G (nst=2, shape=0.5)",
            "config_content": "nst = 2\ngamma = true\ngamma_shape = 0.5\n"
        },
        {
            "name": "Protein Model Override",
            "cli_args": ["--data-type", "protein", "--protein-model", "WAG", "--invariable"],
            "expected_display": "WAG+I (protein=WAG)",
            "config_content": "data_type = protein\nprotein_model = WAG\ninvariable = true\n"
        },
        {
            "name": "Discrete Data Type",
            "cli_args": ["--data-type", "discrete", "--base-freq", "equal"],
            "expected_display": "Mk (nst=1, basefreq=equal)",
            "config_content": "data_type = discrete\nbase_freq = equal\n"
        },
        {
            "name": "Rates Override",
            "cli_args": ["--rates", "equal", "--model", "GTR"],
            "expected_display": "GTR (rates=equal)", 
            "config_content": "rates = equal\nmodel = GTR\n"
        }
    ]
    
    parser = setup_argument_parser()
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\nTest {i}: {test_case['name']}")
        print("-" * 60)
        
        # Test CLI parameters
        cli_args = test_case['cli_args'] + ['dummy_alignment.fasta']
        args_cli = parser.parse_args(cli_args)
        display_cli = build_effective_model_display(args_cli)
        
        # Test config file parameters  
        with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
            f.write("alignment = dummy_alignment.fasta\n")
            f.write(test_case['config_content'])
            config_file = f.name
        
        try:
            args_config = parser.parse_args(['dummy_alignment.fasta'])
            args_config = parse_config(config_file, args_config)
            display_config = build_effective_model_display(args_config)
            
            # Compare results
            cli_match = display_cli == test_case['expected_display']
            config_match = display_config == test_case['expected_display']
            
            print(f"  CLI Display:      {display_cli}")
            print(f"  Config Display:   {display_config}")
            print(f"  Expected:         {test_case['expected_display']}")
            print(f"  CLI Result:       {'✅ PASS' if cli_match else '❌ FAIL'}")
            print(f"  Config Result:    {'✅ PASS' if config_match else '❌ FAIL'}")
            
            if not (cli_match and config_match):
                print(f"  ⚠️  CONSISTENCY ISSUE DETECTED")
                
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
        finally:
            os.unlink(config_file)

def test_parameter_precedence():
    """Test that CLI parameters override config file parameters"""
    
    print("\n\n" + "=" * 80)
    print("PARAMETER PRECEDENCE TESTS (CLI > Config)")
    print("=" * 80)
    
    # Create config file with one set of values
    config_content = """
alignment = config_alignment.fasta
model = HKY
nst = 2
base_freq = equal
gamma = false
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
        f.write(config_content)
        config_file = f.name
    
    try:
        parser = setup_argument_parser()
        
        # Parse with CLI overrides
        args = parser.parse_args([
            'cli_alignment.fasta',  # Override alignment
            '--model', 'GTR',       # Override model
            '--nst', '1',           # Override nst
            '--base-freq', 'empirical',  # Override base_freq
            '--gamma',              # Override gamma
            '--config', config_file
        ])
        
        # Parse config file
        args = parse_config(config_file, args)
        
        # Check that CLI values take precedence
        print("\nPrecedence Test Results:")
        print("-" * 40)
        
        tests = [
            ("alignment", "cli_alignment.fasta", args.alignment),
            ("model", "GTR", args.model),
            ("nst", 1, args.nst),
            ("base_freq", "empirical", args.base_freq),
            ("gamma", True, args.gamma),
        ]
        
        all_passed = True
        for param_name, expected, actual in tests:
            passed = actual == expected
            status = "✅ PASS" if passed else "❌ FAIL"
            print(f"  {param_name:12} = {actual:20} (expected: {expected}) {status}")
            if not passed:
                all_passed = False
        
        print(f"\nOverall Precedence Test: {'✅ PASS' if all_passed else '❌ FAIL'}")
        
        # Test model display with overrides
        display = build_effective_model_display(args)
        print(f"Final Model Display: {display}")
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
    finally:
        os.unlink(config_file)

def test_boolean_parameter_handling():
    """Test boolean parameter conversion from config files"""
    
    print("\n\n" + "=" * 80)
    print("BOOLEAN PARAMETER CONVERSION TESTS")
    print("=" * 80)
    
    boolean_tests = [
        ("true", True),
        ("false", False),
        ("yes", True),
        ("no", False),
        ("on", True),
        ("off", False),
        ("1", True),
        ("0", False),
    ]
    
    parser = setup_argument_parser()
    
    for bool_str, expected in boolean_tests:
        config_content = f"""
alignment = test.fasta
gamma = {bool_str}
invariable = {bool_str}
bootstrap = {bool_str}
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
            f.write(config_content)
            config_file = f.name
        
        try:
            args = parser.parse_args(['test.fasta'])
            args = parse_config(config_file, args)
            
            gamma_ok = args.gamma == expected
            invariable_ok = args.invariable == expected
            bootstrap_ok = args.bootstrap == expected
            
            all_ok = gamma_ok and invariable_ok and bootstrap_ok
            status = "✅ PASS" if all_ok else "❌ FAIL"
            
            print(f"  '{bool_str}' -> {expected}: {status}")
            if not all_ok:
                print(f"    gamma={args.gamma}, invariable={args.invariable}, bootstrap={args.bootstrap}")
                
        except Exception as e:
            print(f"  '{bool_str}' -> {expected}: ❌ ERROR - {e}")
        finally:
            os.unlink(config_file)

def test_config_generation():
    """Test config file generation functionality"""
    
    print("\n\n" + "=" * 80)
    print("CONFIG FILE GENERATION TEST")
    print("=" * 80)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
        config_file = f.name
    
    try:
        # Generate config file
        generate_config_template(config_file)
        
        # Check that file was created and has content
        config_path = Path(config_file)
        if not config_path.exists():
            print("❌ FAIL: Config file not created")
            return
            
        with open(config_file, 'r') as f:
            content = f.read()
        
        # Check for key sections and parameters
        required_sections = [
            "# BASIC INPUT/OUTPUT SETTINGS",
            "# MODEL SETTINGS",
            "# COMPUTATIONAL SETTINGS", 
            "# BAYESIAN-SPECIFIC SETTINGS",
            "# CONVERGENCE CHECKING",
            "# PARALLEL PROCESSING",
            "# CONSTRAINT SETTINGS"
        ]
        
        required_params = [
            "alignment =",
            "model =", 
            "nst =",
            "base_freq =",
            "gamma_shape =",
            "prop_invar =",
            "bayes_model =",
            "mrbayes_parse_timeout =",
            "output_style ="
        ]
        
        missing_sections = [sec for sec in required_sections if sec not in content]
        missing_params = [param for param in required_params if param not in content]
        
        if not missing_sections and not missing_params:
            print("✅ PASS: Config template generated successfully")
            print(f"  File size: {len(content)} characters")
            print(f"  Parameters found: {len([line for line in content.split('\n') if '=' in line and not line.strip().startswith('#')])}")
        else:
            print("❌ FAIL: Config template incomplete")
            if missing_sections:
                print(f"  Missing sections: {missing_sections}")
            if missing_params:
                print(f"  Missing parameters: {missing_params}")
                
    except Exception as e:
        print(f"❌ ERROR: {e}")
    finally:
        if Path(config_file).exists():
            os.unlink(config_file)

def main():
    """Run all parameter consistency tests"""
    
    print("PANDECAY PARAMETER CONSISTENCY TEST SUITE")
    print("Testing CLI, Config File, and Runtime Display Consistency")
    print("=" * 80)
    
    # Run all test suites
    test_model_parameter_consistency()
    test_parameter_precedence()
    test_boolean_parameter_handling()
    test_config_generation()
    
    print("\n\n" + "=" * 80)
    print("PARAMETER CONSISTENCY TESTS COMPLETED")
    print("=" * 80)
    print("\nIf any tests failed, review the parameter handling logic")
    print("All tests should pass for a fully consistent parameter system")

if __name__ == "__main__":
    main()