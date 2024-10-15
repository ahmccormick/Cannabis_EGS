########################################
#Step 1
#########################################
import allel
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

# Step 1: Read VCF file
vcf_file = "/Users/annamccormick/R/cannabis_GEAV/Inputs/Cannabis_sativa_PRJNA734114_filtered.vcf.gz"
callset = allel.read_vcf(vcf_file, fields=['samples', 'calldata/GT'])  # Include only relevant fields for speed

# Extract the genotype data
gt = allel.GenotypeArray(callset['calldata/GT'])

# Step 2: Convert genotypes to genind equivalent (as a numpy array)
# Assuming you want to work with diploid data (like ploidy = 2 in R)
genotypes = gt.to_n_alt(fill=-1)

# Subset the data by individuals (subset_ids)
subset_ids = ["SRR14708197", "SRR14708198", "SRR14708199", "SRR14708200", "SRR14708201",
              "SRR14708202", "SRR14708203", "SRR14708204", "SRR14708205", "SRR14708206",
              "SRR14708207", "SRR14708208", "SRR14708209", "SRR14708210", "SRR14708211",
              "SRR14708212", "SRR14708213", "SRR14708217", "SRR14708222", "SRR14708228",
              "SRR14708229", "SRR14708230", "SRR14708231", "SRR14708232", "SRR14708233",
              "SRR14708234", "SRR14708235", "SRR14708236", "SRR14708237", "SRR14708240",
              "SRR14708243", "SRR14708244", "SRR14708246", "SRR14708248", "SRR14708251",
              "SRR14708254", "SRR14708258", "SRR14708259", "SRR14708272", "SRR14708273",
              "SRR14708275", "SRR14708276", "SRR14708277", "SRR14708278"]

# Find the indices of the samples that match the subset_ids
subset_idx = [i for i, sample in enumerate(callset['samples']) if sample in subset_ids]

# Print the names of the samples being subsetted
print("Subsetted Sample IDs:")
for i in subset_idx:
    print(callset['samples'][i])

# Subset the genotype data to match the selected samples
genotypes_subset = genotypes[:, subset_idx]

# Step 3: Compute genetic distance matrix
dist_matrix = pairwise_distances(genotypes_subset.T, metric='hamming')

# Step 4: Sampling core collection from the subset
# Use only the indices from the subset for sampling
core_indices = np.random.choice(np.arange(len(subset_idx)), size=25, replace=False)

# Map the core indices back to the SRR sample names
core_sample = [subset_ids[i] for i in core_indices]

# Print the core sample names
print("\nCore Sample IDs:")
for sample in core_sample:
    print(sample)

# Save the core sample (SRR names) to CSV
core_sample_data = pd.DataFrame(core_sample, columns=['Sample_ID'])
core_sample_data.to_csv("/Users/annamccormick/R/cannabis_GEAV/Outputs/cannabis_core_n25_2.csv", index=False)