import logging
from ete3 import Tree

# Setup basic configuration for logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Load the tree
try:
    tree = Tree("/Users/yma18/Documents/GitHub/2024_Rhizo_micro_silva/dada2_nochim.tree")  # Ensure this is the correct path to your tree file
    logging.info("Tree loaded successfully.")
except Exception as e:
    logging.error(f"Failed to load tree: {e}")

# Load the list of IDs to keep
ids_to_keep = set()
try:
    with open("/Users/yma18/Documents/GitHub/2024_Rhizo_micro_silva/ASV_ID.txt") as file:  # Ensure this is the correct path to your ID list file
        ids_to_keep = {line.strip() for line in file}
    logging.info(f"IDs to keep loaded: {len(ids_to_keep)} entries.")
except Exception as e:
    logging.error(f"Failed to load IDs: {e}")

# Check if the tree contains any of the IDs to keep
leaf_names = set(tree.get_leaf_names())  # This gets all leaf names as a set of strings
if not leaf_names.intersection(ids_to_keep):  # Checks for intersection between two sets
    logging.warning("None of the IDs to keep are in the tree.")
else:
    logging.info(f"Found {len(leaf_names.intersection(ids_to_keep))} matching IDs in the tree.")

# Prune the tree to keep only the desired IDs
try:
    tree.prune(ids_to_keep, preserve_branch_length=True)
    logging.info("Tree pruning successful.")
except Exception as e:
    logging.error(f"Error during tree pruning: {e}")

# Save the pruned tree
try:
    tree.write(outfile="/Users/yma18/Documents/GitHub/2024_Rhizo_micro_silva/pruned_tree.txt")
    logging.info("Pruned tree saved successfully.")
except Exception as e:
    logging.error(f"Failed to save pruned tree: {e}")

