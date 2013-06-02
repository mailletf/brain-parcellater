
import scipy.sparse as S
import numpy as N

def sparse_matrix_from_dotfile(filename):
    # Load all lines
    lines = [x.strip().split(" ") for x in open(filename)]

    # Find max values for target and seed voxels to determine
    # how big the dense matrix would be
    num_target_voxels = N.max( [ int(x[2]) for x in lines ])
    num_seed_voxels = N.max( [ int(x[0]) for x in lines ])
    #max_value = N.max( [ int(x[-1]) for x in lines ])

    print " > Creating sparse matrix of dimensions %d x %d" % \
        (num_seed_voxels, num_target_voxels)

    mat = S.dok_matrix((num_seed_voxels,num_target_voxels), int)
    for l in lines:
        try:
            sv, tv, v = int(l[0]), int(l[2]), int(l[4])
            mat[sv-1,tv-1] = v
        except Exception as e:
            print e
            print sv, tv, v
            return

    return mat


#def cross_correlate_matrix(mat):



