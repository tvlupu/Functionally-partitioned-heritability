# convert function
# input file prefix is 'gcta', output file prefix is 'conv'

import numpy as np
import struct
import pandas as pd
import sys


def convert_grm_gcta(in_prefix, out_prefix, size=4):

# 1. FILE NAMES
    gcta_bin = in_prefix + ".grm.bin"
    gcta_N = in_prefix + ".grm.N.bin"
    gcta_id = in_prefix + ".grm.id"
    print(f"Input files: {gcta_bin}, {gcta_id}, {gcta_N}")
    out_bin = out_prefix + ".grm.bin"
    out_iid = out_prefix + ".grm.iid"


# 2. READ + PREP DATA
    # Read IDs
    with open(gcta_id) as f:
        id_list = [line.split()[1] for line in f]  # take IID (2nd col)
    n = len(id_list)
 
    # Number of elements in lower triangle (including diagonal)
    num_elements = n * (n + 1) // 2
 
    # Read GCTA GRM binary file
    with open(gcta_bin, "rb") as f:
        grm = np.fromfile(f, dtype=np.float32, count=num_elements)
    
    # Read GCTA N binary file
    with open(gcta_N, "rb") as f:
        Nfile = np.fromfile(f, dtype=np.float32, count=num_elements)
    nm = max(Nfile)
    grm *= Nfile
    
    # Initialize GRM matrix
    grm_mat = np.zeros((n, n), dtype=np.float32)
 
    # Fill lower triangle from GCTA
    tril_indices = np.tril_indices(n)
    grm_mat[tril_indices] = grm
 
    # Make transposed GRM (i.e. upper triangle) because Python writes things weird later
    grm_mat = grm_mat.T 
        # for symmetrical GRM:
        #+ grm_mat - np.diag(np.diag(grm_mat))

# 3. WRITE .IID FILE (assume 2 columns: FID, IID)
    ids = pd.read_csv(gcta_id, delim_whitespace=True, header=None)
    ids.columns = ['FID', 'IID']
    id_list = ids['IID'].astype(str)
    id_list.to_csv(out_iid, sep="\t", index=False, header=False)
 
    print(f"Completed writing the iid file: {out_iid}")


# 4. WRITE .BIN FILE
    # Get lower triangle including diagonal
    grm_values = grm_mat.astype(np.float32)
 
    with open(out_bin, "wb") as f:
        # Write number of individuals (as 4-byte int)
        f.write(np.array([n], dtype=np.int32).tobytes())
 
        # Write number of markers (as 4-byte float)
        f.write(np.array([nm], dtype=np.float32).tobytes())
 
        # Write GRM values (as 4-byte float)
        f.write(grm_values[np.triu_indices(n)].tobytes()) # remove the 0s in the matrix while writing string
 
    print(f"Completed writing the bin file: {out_bin}")
    return grm_mat

#### convert
for arg in sys.argv:
    print(arg)
convert_grm_gcta(sys.argv[1], sys.argv[2])
