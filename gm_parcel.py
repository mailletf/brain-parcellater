
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

    #num_target_voxels = 100
    #num_seed_voxels = 25

    mat = [S.dok_matrix((1,num_target_voxels), int) for x in xrange(num_seed_voxels)]
    for l in lines:
        try:
            sv, tv, v = int(l[0]), int(l[2]), int(l[4])
            #if tv > num_target_voxels: break
            #if sv >= num_seed_voxels: continue
            mat[sv-1][0,tv-1] = v
        except Exception as e:
            print e
            print sv, tv, v
            return

    return mat


def cross_correlate_matrix(mat):
    # Cache the norm of all vectors
    print " > Computing norms..."
    norms = [N.sqrt(mat[i].multiply(mat[i]).sum(1)) for i in xrange(len(mat))]

    # Calculate cross correlation
    cc_mat = N.zeros((len(mat), len(mat)))

    print " > Computing cross-correlation matrix..."
    import time
    for i in xrange(len(mat)):
        print "%d/%d at %s" % (i, len(mat), time.ctime())
        for j in xrange(len(mat)):
            if cc_mat[i,j] != 0: continue
            cc_mat[i,j] = mat[i].dot(mat[j].T)[0,0] / (norms[i] * norms[j])
            cc_mat[j,i] = cc_mat[i,j]

    return cc_mat




