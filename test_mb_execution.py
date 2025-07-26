#!/usr/bin/env python3
"""
Test MrBayes execution to understand the issue.
"""

import subprocess
import tempfile
import shutil
from pathlib import Path

def test_mrbayes_execution():
    """Test MrBayes execution like the monolithic system."""
    
    # Create temporary directory
    temp_dir = Path("test_mb_temp")
    temp_dir.mkdir(exist_ok=True)
    
    try:
        # Copy test data
        shutil.copy("test_simple.nex", temp_dir / "test_simple.nex")
        
        # Create simple MrBayes script
        script_content = """#NEXUS

BEGIN MRBAYES;
    execute test_simple.nex;
    set autoclose=yes nowarn=yes;
    lset nst=1 rates=equal;
    mcmc ngen=100 printfreq=100 samplefreq=100 nchains=1 savebrlens=yes;
    sump;
    sumt;
END;
"""
        
        script_file = temp_dir / "test.nex"
        script_file.write_text(script_content)
        
        print(f"Created script: {script_file}")
        print(f"Working directory: {temp_dir}")
        print(f"Files in directory: {list(temp_dir.glob('*'))}")
        
        # Run MrBayes like monolithic system
        cmd = ["mb", script_file.name]
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {temp_dir}")
        
        result = subprocess.run(
            cmd,
            cwd=str(temp_dir),
            capture_output=True,
            text=True,
            timeout=60
        )
        
        print(f"Return code: {result.returncode}")
        print(f"Stdout length: {len(result.stdout)}")
        print(f"Stderr length: {len(result.stderr)}")
        
        if result.stdout:
            print("=== STDOUT ===")
            print(result.stdout[-1000:])  # Last 1000 chars
        
        if result.stderr:
            print("=== STDERR ===")
            print(result.stderr[-500:])  # Last 500 chars
        
        # Check output files
        print("=== OUTPUT FILES ===")
        output_files = list(temp_dir.glob("*"))
        for f in output_files:
            print(f"  {f.name}: {f.stat().st_size} bytes")
        
        # Look for MrBayes output patterns
        run_files = list(temp_dir.glob("*.run*.p")) + list(temp_dir.glob("*.run*.t"))
        con_files = list(temp_dir.glob("*.con.tre"))
        lstat_files = list(temp_dir.glob("*.lstat"))
        
        print(f"Run files: {len(run_files)}")
        print(f"Consensus files: {len(con_files)}")
        print(f"Lstat files: {len(lstat_files)}")
        
        if con_files:
            print("SUCCESS: MrBayes generated consensus tree!")
            return True
        else:
            print("FAILED: No consensus tree generated")
            return False
            
    except Exception as e:
        print(f"ERROR: {e}")
        return False
    finally:
        # Clean up
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    test_mrbayes_execution()