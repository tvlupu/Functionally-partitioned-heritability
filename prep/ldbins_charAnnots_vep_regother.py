import pandas as pd # type: ignore
import sys

# parameters
chrom=sys.argv[1]
bin_num = sys.argv[2]
curr_annot = sys.argv[3]
ld_bin = sys.argv[4]
maf_low = sys.argv[5]
maf_high = sys.argv[6]
other_bin = sys.argv[7]

# load data
print(f'LD score file: chr{chrom}_score.ld')
ldscores = pd.read_csv(f'chr{chrom}_score.ld', sep=' ', header=None)
ldscores.columns = ['SNP', 'chr', 'bp', 'MAF', 'mean_rsq', 'var_num', 'max_rsq', 'ldscore']

# filter to chromosome and sort descending by ldscore
# drop duplicate variants
ldscores = ldscores.loc[ldscores['chr'] == int(chrom)].sort_values(by='ldscore', ascending=False)
ldscores = ldscores.drop_duplicates(subset=['SNP'], keep='first')

# drop duplicate variants, keeping the first one (highest ldscore for ldscores, cod>rna>prox>intr>reg>other for rest)
annotsfile_reg = pd.read_csv(f'chr{chrom}_{curr_annot}', sep=' ', header=None)
annotsfile_reg.columns = ['SNP', 'ANNOT']
annotsfile_reg = annotsfile_reg.drop_duplicates(subset=['SNP'], keep='first')

# merge ldscore file with annotsfile_reg, keeping all variants from annotsfile_reg only
ldscores_subset = pd.merge(ldscores, annotsfile_reg, on='SNP', how='right')
print(f'Num variants {curr_annot} in chromosome:', len(ldscores_subset))

# +
# subset to relevant annotation
temp_dt = ldscores_subset[(ldscores_subset['ANNOT'] == curr_annot) & (ldscores_subset['MAF'] > float(maf_low)) & (ldscores_subset['MAF'] <= float(maf_high))]

print(f"LD bin: {ld_bin}")
print(f"MAF range: {maf_low} to {maf_high}")
print(f"variants annotation: {curr_annot}")
print('Num variants after filtering on MAF and annotation:', len(temp_dt))


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
# -

# add variant annotations
annotsfile_nonreg = pd.read_csv(f'chr{chrom}_vep', header=None)
annotsfile_nonreg.columns = ['SNP']
annotsfile_nonreg['ANNOT'] = "excl"

# +
#ld score subset for everything that isn't reg or previous annotations
annotsfile_all = pd.concat([annotsfile_nonreg, annotsfile_reg], ignore_index=True)
annotsfile_all = annotsfile_all.drop_duplicates(subset=['SNP'], keep='first')

ldscores_subset = pd.merge(ldscores, annotsfile_all, on='SNP', how='left')

# replace missing annotations with "other"
ldscores_subset['ANNOT'] = ldscores_subset['ANNOT'].fillna("other")
print('Num variants in chromosome:', len(ldscores_subset))

# +
# subset to "other" annotation
temp_dt = ldscores_subset[(ldscores_subset['ANNOT'] == "other") & (ldscores_subset['MAF'] > float(maf_low)) \
                         & (ldscores_subset['MAF'] <= float(maf_high))]

print(f"LD bin: {ld_bin}")
print(f"MAF range: {maf_low} to {maf_high}")
print(f"variants annotation: other")
print('Num variants after filtering on MAF and annotation:', len(temp_dt))


# LD binning based on high or low LD
output_file = f"chr{chrom}_{other_bin}_vep.varlist"
q2 = temp_dt['ldscore'].quantile(1/2)

if ld_bin == 'lowbin':
    varlist = temp_dt[temp_dt['ldscore'] < q2]['SNP']
else:  # high LD bin
    varlist = temp_dt[temp_dt['ldscore'] >= q2]['SNP']

# Write variant list
print('Num variants in bin ', {other_bin}, ': ',  len(varlist))
varlist.to_csv(output_file, index=False, header=False, sep='\t', quoting=3)
print(f"Variant list written to: {output_file}")
