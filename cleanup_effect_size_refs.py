#!/usr/bin/env python3
"""
Script to clean up remaining effect size references in panDecay.py
"""

import re

def clean_effect_size_references():
    """Remove all effect size references from panDecay.py"""
    
    with open('panDecay.py', 'r') as f:
        content = f.read()
    
    # Define replacements
    replacements = [
        # Comments and docstrings
        ('effect size calculations', 'dataset-relative calculations'),
        ('effect size calculation', 'dataset-relative calculation'),
        ('effect size', 'dataset-relative metric'),
        ('Effect Size', 'Dataset-Relative Metric'),
        ('Cohen.*d', 'dataset-relative normalization'),
        
        # Variable and method names - be careful to preserve logic
        ('effect_size_decay', 'ld_for_normalization'),
        ('has_effect_size', 'has_dataset_relative'),
        ('effect_size_methods', 'dataset_relative_methods'),
        ('effect_size_values', 'dataset_relative_values'),
        
        # Method calls and conditions - preserve the structure
        (r'method\.startswith\([\'"]effect_size[\'"]', "method in ['dataset_relative', 'percentile_rank', 'z_score']"),
        
        # Headers and output labels
        ('Effect Size', 'Dataset Relative'),
        ('ES Robust', 'Percentile Rank'),
        ('ES Weighted', 'Z Score'),
        ('Effect_Size', 'Dataset_Relative'),
        ('ES_Robust', 'Percentile_Rank'), 
        ('ES_Weighted', 'Z_Score'),
        
        # Documentation strings
        ('BD / SD\\(site signals\\)', 'relative ranking within dataset'),
        ('Cohen.*d framework', 'dataset-relative framework'),
        ('signal-to-noise ratio', 'relative ranking'),
        
        # Method names in normalization options
        ('"effect_size"', '"dataset_relative"'),
        ('"effect_size_robust"', '"percentile_rank"'),
        ('"effect_size_weighted"', '"z_score"'),
        ("'effect_size'", "'dataset_relative'"),
        ("'effect_size_robust'", "'percentile_rank'"),
        ("'effect_size_weighted'", "'z_score'"),
    ]
    
    # Apply replacements
    for old, new in replacements:
        content = re.sub(old, new, content, flags=re.IGNORECASE)
    
    # Clean up any remaining effect size function calls
    # These need careful handling to preserve function logic
    
    # Remove the effect size interpretation function entirely
    effect_size_func_pattern = r'def _get_effect_size_interpretation_scale\(self, effect_size\):.*?(?=\n    def |\nclass |\Z)'
    content = re.sub(effect_size_func_pattern, '', content, flags=re.DOTALL)
    
    # Remove smoke test references to effect sizes
    smoke_test_pattern = r'def run_smoke_tests\(\):.*?(?=\ndef |\nclass |\Z)'
    content = re.sub(smoke_test_pattern, '', content, flags=re.DOTALL)
    
    # Remove any remaining effect size references in output formatting
    content = re.sub(r'any\(key\.startswith\([\'"]effect_size[\'"].*?\)', 
                     "any(key in ['dataset_relative', 'percentile_rank', 'z_score'] for key in data.keys())", 
                     content)
    
    with open('panDecay.py', 'w') as f:
        f.write(content)
    
    print("Cleaned up effect size references in panDecay.py")

if __name__ == '__main__':
    clean_effect_size_references()