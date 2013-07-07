import os, math, time
from optparse import OptionParser
import pylab
import numpy as N
import random


from scipy.sparse.linalg import lobpcg
from pyamg import smoothed_aggregation_solver
#from helper import trimesh, graph_laplacian

import conn_ccorr_mat as CC
import io_util as io

def get_fake_data():
    m = N.zeros((50,50))
    m[40:,:10] = 1
    m[:10,40:] = 1
    
    indexes = range(len(m))
    random.shuffle(indexes)

    return (m, indexes)

def do_spec_reorder(W):
    D = N.zeros(W.shape)
    for i in xrange(len(W)):
        D[i,i] = N.sum(W[i,:])

    E = N.zeros(W.shape)
    for i in xrange(len(D)):
        E[i,i] = 1.0/N.sqrt(D[i,i])

    w2 = E.dot((D-W)).dot(E)


    # calcule eigenvalues and get the vector for the
    # 2nd smallest value
    eigenvalues, eigenvectors = N.linalg.eigh(w2)
    sorted_values = N.argsort(eigenvalues)
    v = eigenvectors[:, sorted_values[1]]

    v2 = E.dot(v)
    permutation_vector = N.argsort(v2)
    
    fig = pylab.figure()
    f1 = fig.add_subplot(1,2,1)
    #f1.set_title(text="Cross-correlation matrix")
    f1.imshow(W, interpolation="nearest")

    f2 = fig.add_subplot(1,2,2)
    #f2.set_title(text="Spectral reordering")
    Wr1 = W[permutation_vector]
    Wr1 = Wr1[:,permutation_vector]
    f2.imshow(Wr1, interpolation="nearest")
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

