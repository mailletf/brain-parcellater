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
    
    parser.add_option_group(CC.get_option_parser_group(parser))

    (options, args) = parser.parse_args()

    if not options.conn_mat_filename:
        raise Exception("Must specify filename. Run with --help for instructions.")

    run(options, args)

