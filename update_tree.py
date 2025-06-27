#!/usr/bin/env python3
"""
Simple script to update the annotated tree with Bayesian results.
"""

# Read the current tree
with open("annotated_tree_combined.nwk", "r") as f:
    tree_content = f.read().strip()

# Read the Bayesian results from ml_decay_indices.txt
bayes_data = {}
with open("ml_decay_indices.txt", "r") as f:
    lines = f.readlines()
    in_data = False
    for line in lines:
        line = line.strip()
        if line.startswith("Clade_ID"):
            in_data = True
            continue
        if in_data and line and not line.startswith("-"):
            parts = line.split('\t')
            if len(parts) >= 9:
                clade_id = parts[0]
                bayes_decay = parts[6]
                bayes_factor = parts[7]
                bayes_data[clade_id] = (bayes_decay, bayes_factor)

# Update the tree with Bayesian results
for clade_id, (bd, bf) in bayes_data.items():
    # Find the clade annotation in the tree
    old_pattern = f"'{clade_id}|AU:"
    if old_pattern in tree_content:
        # Find the closing quote
        start = tree_content.find(old_pattern)
        end = tree_content.find("'", start + 1)
        old_annotation = tree_content[start:end+1]
        
        # Add Bayesian data before the closing quote
        new_annotation = old_annotation[:-1] + f"|BD:{bd}|BF:{bf}'"
        tree_content = tree_content.replace(old_annotation, new_annotation)

# Write the updated tree
with open("annotated_tree_combined_with_bayes.nwk", "w") as f:
    f.write(tree_content)
    f.write("\n")

print("Created annotated_tree_combined_with_bayes.nwk with Bayesian results included")

# Also create separate trees for Bayes decay and Bayes factor only
# Create Bayes decay tree
bd_tree = tree_content
for clade_id, (bd, bf) in bayes_data.items():
    old_pattern = f"'{clade_id}|"
    if old_pattern in bd_tree:
        # Find the full annotation
        start = bd_tree.find(old_pattern)
        end = bd_tree.find("'", start + 1)
        old_annotation = bd_tree[start:end+1]
        # Replace with just Bayes decay
        new_annotation = f"'{clade_id}|BD:{bd}'"
        bd_tree = bd_tree.replace(old_annotation, new_annotation)

with open("annotated_tree_bayes_decay.nwk", "w") as f:
    f.write(bd_tree)
    f.write("\n")

print("Created annotated_tree_bayes_decay.nwk with only Bayes decay values")

# Create Bayes factor tree
bf_tree = tree_content
for clade_id, (bd, bf) in bayes_data.items():
    old_pattern = f"'{clade_id}|"
    if old_pattern in bf_tree:
        # Find the full annotation
        start = bf_tree.find(old_pattern)
        end = bf_tree.find("'", start + 1)
        old_annotation = bf_tree[start:end+1]
        # Replace with just Bayes factor
        new_annotation = f"'{clade_id}|BF:{bf}'"
        bf_tree = bf_tree.replace(old_annotation, new_annotation)

with open("annotated_tree_bayes_factor.nwk", "w") as f:
    f.write(bf_tree)
    f.write("\n")

print("Created annotated_tree_bayes_factor.nwk with only Bayes factor values")