#!/usr/bin/env python3
"""
Debug script to test AU test execution and parsing.
"""

import tempfile
import logging
from pathlib import Path
from src.utils.tree_manager import TreeManager

# Setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def debug_au_test():
    """Test AU test execution to see what's being generated."""
    
    # Create temp directory for debugging  
    temp_dir = Path(tempfile.mkdtemp(prefix="debug_au_"))
    print(f"Debug temp directory: {temp_dir}")
    
    # Setup tree manager
    tree_manager = TreeManager(temp_dir, debug=True)
    
    # Copy alignment to temp directory
    alignment_file = temp_dir / "alignment.nex"
    alignment_file.write_text(Path("test_simple.nex").read_text())
    
    try:
        print("Building ML tree...")
        
        # Build ML tree first
        ml_result = tree_manager.build_ml_tree(
            alignment_file=alignment_file,
            model_settings={},
            timeout=300
        )
        
        ml_tree = ml_result['tree_file']
        print(f"ML tree: {ml_tree}")
        
        # Build a constraint tree
        print("Building constraint tree...")
        constraint_taxa = ['Homo_sapiens', 'Pan_troglodytes', 'Gorilla_gorilla']
        
        constraint_result = tree_manager.test_constraint(
            alignment_file=alignment_file,
            constraint_taxa=constraint_taxa,
            constraint_id="test_constraint",
            model_settings={},
            timeout=300
        )
        
        constraint_tree = constraint_result.get('tree_file')
        print(f"Constraint tree: {constraint_tree}")
        
        if constraint_tree and constraint_tree.exists():
            print("Running AU test...")
            
            # Run AU test
            au_results = tree_manager.run_au_test(
                alignment_file=alignment_file,
                tree_files=[ml_tree, constraint_tree],
                timeout=600
            )
            
            print(f"AU results structure: {au_results}")
            print(f"AU results type: {type(au_results)}")
            if isinstance(au_results, dict):
                for key, value in au_results.items():
                    print(f"  {key}: {value} (type: {type(value)})")
        else:
            print("No constraint tree generated - skipping AU test")
        
        # Check what files exist in temp directory
        print(f"\nFiles in temp directory {temp_dir}:")
        for file_path in temp_dir.rglob("*"):
            if file_path.is_file():
                print(f"  {file_path.relative_to(temp_dir)} ({file_path.stat().st_size} bytes)")
                
                # If it's an AU test log, show some content
                if 'au_test' in file_path.name.lower():
                    content = file_path.read_text()
                    print(f"    Content preview: {content[:500]}...")
        
    except Exception as e:
        print(f"Error during debugging: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        print(f"\nDebug temp directory preserved: {temp_dir}")

if __name__ == "__main__":
    debug_au_test()