#!/usr/bin/env python3
"""
CLI Parameter Consistency Test Suite for panDecay
Tests CLI parameters vs config file parameters by running actual commands
"""

import subprocess
import tempfile
import os
from pathlib import Path

def create_test_alignment():
    """Create a minimal test alignment file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">seq1\nATCGATCGATCGATCG\n>seq2\nATCGATCGATCGATCG\n>seq3\nATCGATCGATCGATCG\n>seq4\nATCGATCGATCGATCG\n")
        return f.name

def run_pandecay_command(args, config_content=None):
    """Run panDecay command and capture runtime parameters display"""
    
    alignment_file = create_test_alignment()
    
    try:
        if config_content:
            # Create config file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
                f.write(f"alignment = {alignment_file}\n")
                f.write(config_content)
                config_file = f.name
            
            cmd = ["python3", "-m", "src.main", "--config", config_file]
        else:
            cmd = ["python3", "-m", "src.main"] + args + [alignment_file]
        
        # Run command and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        # Extract the runtime parameters box
        output_lines = result.stdout.split('\n') + result.stderr.split('\n')
        
        # Find the model line in the runtime parameters box
        model_line = None
        for line in output_lines:
            if "Model:" in line and "│" in line:
                model_line = line.strip()
                break
        
        return model_line, result.returncode
        
    except subprocess.TimeoutExpired:
        return "TIMEOUT", -1
    except Exception as e:
        return f"ERROR: {e}", -1
    finally:
        # Cleanup
        if os.path.exists(alignment_file):
            os.unlink(alignment_file)
        if config_content and 'config_file' in locals() and os.path.exists(config_file):
            os.unlink(config_file)

def extract_model_from_line(line):
    """Extract model string from runtime parameters line"""
    if not line or "Model:" not in line:
        return None
    
    # Extract text between "Model:" and the next │
    parts = line.split("Model:")
    if len(parts) < 2:
        return None
    
    model_part = parts[1].split("│")[0].strip()
    return model_part

def test_model_parameter_consistency():
    """Test model parameter display consistency between CLI and config"""
    
    print("=" * 80)
    print("CLI vs CONFIG PARAMETER CONSISTENCY TESTS")
    print("=" * 80)
    
    test_cases = [
        {
            "name": "NST=1 Override (Original Issue)",
            "cli_args": ["--nst", "1", "--base-freq", "empirical"],
            "config_content": "nst = 1\nbase_freq = empirical\n",
            "expected_contains": ["JC", "nst=1", "basefreq=empirical"]
        },
        {
            "name": "NST=2 with Gamma",
            "cli_args": ["--nst", "2", "--gamma"],
            "config_content": "nst = 2\ngamma = true\n", 
            "expected_contains": ["HKY", "nst=2"]
        },
        {
            "name": "Protein Model Override",
            "cli_args": ["--data-type", "protein", "--protein-model", "WAG"],
            "config_content": "data_type = protein\nprotein_model = WAG\n",
            "expected_contains": ["WAG", "protein=WAG"]
        },
        {
            "name": "Discrete Data Type", 
            "cli_args": ["--data-type", "discrete"],
            "config_content": "data_type = discrete\n",
            "expected_contains": ["Mk", "nst=1"]
        },
        {
            "name": "Default GTR",
            "cli_args": [],
            "config_content": "",
            "expected_contains": ["GTR"]
        }
    ]
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\nTest {i}: {test_case['name']}")
        print("-" * 60)
        
        # Test CLI version
        cli_model_line, cli_code = run_pandecay_command(test_case['cli_args'])
        cli_model = extract_model_from_line(cli_model_line)
        
        # Test config version
        config_model_line, config_code = run_pandecay_command([], test_case['config_content'])
        config_model = extract_model_from_line(config_model_line)
        
        print(f"  CLI Model Display:    {cli_model}")
        print(f"  Config Model Display: {config_model}")
        
        # Check if models match
        if cli_model and config_model:
            models_match = cli_model == config_model
            print(f"  Models Match:         {'✅ YES' if models_match else '❌ NO'}")
            
            # Check if expected content is present
            cli_has_expected = cli_model and all(expected in cli_model for expected in test_case['expected_contains'])
            config_has_expected = config_model and all(expected in config_model for expected in test_case['expected_contains'])
            
            print(f"  CLI Has Expected:     {'✅ YES' if cli_has_expected else '❌ NO'}")
            print(f"  Config Has Expected:  {'✅ YES' if config_has_expected else '❌ NO'}")
            
            if models_match and cli_has_expected and config_has_expected:
                print(f"  Overall Result:       ✅ PASS")
            else:
                print(f"  Overall Result:       ❌ FAIL")
                if not models_match:
                    print(f"    Issue: CLI and config produce different model displays")
                if not cli_has_expected:
                    print(f"    Issue: CLI missing expected content: {test_case['expected_contains']}")
                if not config_has_expected:
                    print(f"    Issue: Config missing expected content: {test_case['expected_contains']}")
        else:
            print(f"  Overall Result:       ❌ FAIL - Could not extract model information")
            if cli_code != 0:
                print(f"    CLI Error (code {cli_code}): {cli_model_line}")
            if config_code != 0:
                print(f"    Config Error (code {config_code}): {config_model_line}")

def test_config_file_generation():
    """Test config file generation and parameter inclusion"""
    
    print("\n\n" + "=" * 80)
    print("CONFIG FILE GENERATION TEST")
    print("=" * 80)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
        config_file = f.name
    
    try:
        # Test config generation command
        cmd = ["python3", "-m", "src.main", "--generate-config", config_file]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0:
            print("✅ Config generation command succeeded")
            
            # Check file contents
            with open(config_file, 'r') as f:
                content = f.read()
            
            # Check for newly added parameters
            new_params = [
                "nst =", "base_freq =", "gamma_shape =", "prop_invar =",
                "rates =", "parsmodel =", "starting_tree =", "paup_block =", 
                "temp =", "bayes_model =", "mpi_processors =", 
                "mrbayes_parse_timeout =", "test_branches =", 
                "constraint_file =", "output_style ="
            ]
            
            missing_params = [param for param in new_params if param not in content]
            
            if not missing_params:
                print("✅ All required parameters present in generated config")
                print(f"  Config file size: {len(content)} characters")
            else:
                print(f"❌ Missing parameters in config: {missing_params}")
                
        else:
            print(f"❌ Config generation failed (code {result.returncode})")
            print(f"  Error: {result.stderr}")
            
    except Exception as e:
        print(f"❌ Exception during config generation test: {e}")
    finally:
        if os.path.exists(config_file):
            os.unlink(config_file)

def test_parameter_precedence():
    """Test that CLI parameters override config file parameters"""
    
    print("\n\n" + "=" * 80)
    print("PARAMETER PRECEDENCE TEST (CLI > Config)")
    print("=" * 80)
    
    # Create config with one model setting
    config_content = """
model = HKY
nst = 2
base_freq = equal
"""
    
    # CLI args that should override config
    cli_args = ["--model", "GTR", "--nst", "1", "--base-freq", "empirical"]
    
    alignment_file = create_test_alignment()
    
    try:
        # Create config file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False) as f:
            f.write(f"alignment = {alignment_file}\n")
            f.write(config_content)
            config_file = f.name
        
        # Run with CLI overrides
        cmd = ["python3", "-m", "src.main"] + cli_args + ["--config", config_file]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        # Extract model display
        output_lines = result.stdout.split('\n') + result.stderr.split('\n')
        model_line = None
        for line in output_lines:
            if "Model:" in line and "│" in line:
                model_line = line.strip()
                break
        
        model_display = extract_model_from_line(model_line)
        print(f"Model Display with CLI Overrides: {model_display}")
        
        # Check that CLI values took precedence
        if model_display:
            has_jc = "JC" in model_display  # nst=1 should make it JC
            has_nst1 = "nst=1" in model_display
            has_empirical = "basefreq=empirical" in model_display
            
            if has_jc and has_nst1 and has_empirical:
                print("✅ CLI parameters correctly override config file")
            else:
                print("❌ CLI parameters did not override config file")
                print(f"  Expected: JC with nst=1, basefreq=empirical")
                print(f"  Got: {model_display}")
        else:
            print(f"❌ Could not extract model display (code {result.returncode})")
            
    except Exception as e:
        print(f"❌ Exception during precedence test: {e}")
    finally:
        if os.path.exists(alignment_file):
            os.unlink(alignment_file)
        if 'config_file' in locals() and os.path.exists(config_file):
            os.unlink(config_file)

def main():
    """Run all CLI parameter consistency tests"""
    
    print("PANDECAY CLI PARAMETER CONSISTENCY TEST SUITE")
    print("Testing actual CLI behavior and config file integration")
    print("=" * 80)
    
    # Run test suites
    test_model_parameter_consistency()
    test_config_file_generation()
    test_parameter_precedence()
    
    print("\n\n" + "=" * 80)
    print("CLI PARAMETER CONSISTENCY TESTS COMPLETED")
    print("=" * 80)
    print("\nThese tests validate the fix for the original --nst parameter issue")
    print("and ensure complete consistency between CLI and config file parameters.")

if __name__ == "__main__":
    main()