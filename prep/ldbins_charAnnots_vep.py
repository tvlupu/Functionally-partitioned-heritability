import pandas as pd # type: ignore
import sys

# Parameters
chrom=sys.argv[1]
bin_num = sys.argv[2]
curr_annot = sys.argv[3]
ld_bin = sys.argv[4]
maf_low = sys.argv[5]
maf_high = sys.argv[6]

# Load and prepare data
print(f'LD score file: chr{chrom}_score.ld')
ldscores = pd.read_csv(f'chr{chrom}_score.ld', sep=' ', header=None)
ldscores.columns = ['SNP', 'chr', 'bp', 'MAF', 'mean_rsq', 'var_num', 'max_rsq', 'ldscore']

# Filter to chromosome and sort descending by ldscore
ldscores = ldscores.loc[ldscores['chr'] == int(chrom)].sort_values(by='ldscore', ascending=False)

# Add variant annotations
annotsfile = pd.read_csv(f'cat_vep_chr{chrom}.txt', sep=' ', header=None)
annotsfile.columns = ['VAR', 'ANNOT']

# Remove duplicate vars, keeping the first one (highest ldscore for ldscores, cod>rna>prox>intr>reg>other for rest)
annotsfile = annotsfile.drop_duplicates(subset=['SNP'], keep='first')
ldscores_subset = ldscores.drop_duplicates(subset=['SNP'], keep='first')

ldscores_subset = pd.merge(ldscores_subset, annotsfile, on='SNP', how='left')

# Replace missing annotations with "other"
ldscores_subset['ANNOT'] = ldscores_subset['ANNOT'].fillna("other")
print('Num variants in chromosome:', len(ldscores_subset))

# +
# Subset to relevant annotation
temp_dt = ldscores_subset[(ldscores_subset['ANNOT'] == curr_annot) & (ldscores_subset['MAF'] > float(maf_low)) \
                         & (ldscores_subset['MAF'] <= float(maf_high))]

print(f"LD bin: {ld_bin}")
print(f"MAF range: {maf_low} to {maf_high}")
print(f"Variant annotation: {curr_annot}")
print('Num variantss after filtering on MAF and annotation:', len(temp_dt))


# LD binning based on high or low LD
output_file = f"chr{chrom}_{bin_num}_vep.varlist"
q2 = temp_dt['ldscore'].quantile(1/2)

if ld_bin == 'lowbin':
    varlist = temp_dt[temp_dt['ldscore'] < q2]['SNP']
else:  # high LD bin
    varlist = temp_dt[temp_dt['ldscore'] >= q2]['SNP']

# Write variant list
print('Num variants in bin ', {bin_num}, ': ',  len(varlist))
varlist.to_csv(output_file, index=False, header=False, sep='\t', quoting=3)
print(f"Variant list written to: {output_file}")
