import pandas as pd

pd.set_option('display.max_rows', 150)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1000)

import numpy as np
import pandas as pd

def calculate_log2_fold_change(treatment_counts, control_counts):
    """
    Calculate the observed log2 fold change for each guide in the treatment and control conditions.

    Parameters:
    - treatment_counts: A DataFrame where rows represent guides and columns represent replicates for the treatment condition.
    - control_counts: A DataFrame where rows represent guides and columns represent replicates for the control condition.

    Returns:
    - log2_fold_changes: A Series containing the log2 fold changes for each guide.
    """

    # Pseudocount to avoid log(0)
    pseudocount = 32

    # Number of replicates
    R_T = treatment_counts.shape[1]  # Number of treatment replicates
    R_C = control_counts.shape[1]    # Number of control replicates

    # Log-transform and median-normalize treatment and control counts
    Tg = np.log2(treatment_counts + pseudocount) - np.median(np.log2(treatment_counts + pseudocount), axis=0)
    Cg = np.log2(control_counts + pseudocount) - np.median(np.log2(control_counts + pseudocount), axis=0)

    # Calculate the mean of normalized replicates for treatment and control
    mean_Tg = Tg.mean(axis=1)  # Mean across treatment replicates
    mean_Cg = Cg.mean(axis=1)  # Mean across control replicates

    # Compute the observed log2 fold change
    log2_fold_changes = mean_Tg - mean_Cg

    return log2_fold_changes

# Example usage
# Sample data (replace these with your actual counts)
treatment_counts = pd.DataFrame({
    'replicate_1': [100, 200, 150],
    'replicate_2': [110, 210, 160],
    'replicate_3': [105, 205, 155],
})

control_counts = pd.DataFrame({
    'replicate_1': [95, 190, 140],
    'replicate_2': [100, 185, 145],
    'replicate_3': [97, 188, 142],
})

# Calculate log2 fold changes
log2_fold_changes = calculate_log2_fold_change(treatment_counts, control_counts)

print("Log2 Fold Changes:")
print(log2_fold_changes)
