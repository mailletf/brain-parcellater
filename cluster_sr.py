import os, math, time
from optparse import OptionParser
import pylab
import numpy as N

import conn_ccorr_mat as CC
import io_util as io


def do_spec_reorder(B):
    C = B + 1

    Q = -C
    for i in xrange(Q.shape[0]):
        Q[i,i] = -(N.sum(Q[i,:]) + Q[i,i])
    
    t = N.zeros(Q.shape)
    for i in xrange(len(Q)):
        for j in xrange(len(Q)):
            t[i,j] = 1.0 / N.sqrt( N.sum(C[:,j]) )

    D = t * Q * t

    # calcule eigenvalues and get the vector for the
    # 2nd smallest value
    eigenvalues, eigenvectors = N.linalg.eig(D)
    sorted_values = N.argsort(eigenvalues)
    v = eigenvectors[sorted_values[1]]
    v2 = t.dot(v)
    
    permutation_vector = N.argsort(v2)
    print permutation_vector
    print B[permutation_vector]
    
    pylab.imshow(B)
    pylab.show()
    pylab.imshow(B[permutation_vector])
    pylab.show()


def run(options, args):
    cc_mat = CC.get_ccmat(options, args)
    do_spec_reorder(cc_mat)


if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("", "--conn-mat-file", dest="conn_mat_filename",
                              help="connectivity matrix filename")
    parser.add_option("", "--voxel-coords-file", dest="voxel_coords_filename",
                              help="voxel coordinates filename")
    parser.add_option("", "--ascii-out-file", dest="ascii_out_filename",
                              help="ascii matrix output filename")
    parser.add_option("-s", "--samples", dest="samples", default=0, type="int",
                              help="number of samples (ie number of streamlines started for each voxel)")
    parser.add_option("-t", "--threshold", dest="threshold", default=0, type="float",
                              help="Consider tracts connected if P(connected) > THRESHOLD %, where P(connected)=matrix_value/samples")

    parser.add_option("", "--max-seed-voxels", dest="max_seed_voxels", default=0, type="int",
                              help="Limit the number of seed voxels used to the first N")
    parser.add_option("", "--max-target-voxels", dest="max_target_voxels", default=0, type="int",
                              help="Limit the number of target voxels used to the first N")

    parser.add_option("", "--display-corr-matrix", action="store_true", dest="display_corr_matrix")

    (options, args) = parser.parse_args()

    if not options.conn_mat_filename:
        raise Exception("Must specify filename. Run with --help for instructions.")

    run(options, args)

