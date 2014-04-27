# conn_ccorr_mat.py
# Copyright (C) 2013 Francois Maillet, Martha Shiell
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#

import os, math, time
from optparse import OptionGroup
from sklearn.cluster import KMeans
from scipy.spatial.distance import euclidean
import scipy.sparse as S
import pylab
import numpy as N


def corr_matrix_cache_filename(filename, samples=0, threshold=0, max_seed_voxels=0, max_target_voxels=0, binarize=True):
    return "%s_%d_%0.4f_%d_%d_%d.npy" % (filename, samples, threshold, max_seed_voxels, max_target_voxels, binarize)


def sparse_matrix_from_dotfile(filename, samples, threshold, max_seed_voxels=0, max_target_voxels=0, binarize=True):
    has_max_sv = max_seed_voxels > 0
    has_max_tv = max_target_voxels > 0

    # Load all lines
    lines = [x.strip().split(" ") for x in open(filename)]

    # Find max values for target and seed voxels to determine
    # how big the dense matrix would be
    all_seeds = [ int(x[0]) for x in lines ]
    num_seed_voxels = N.max(all_seeds)
    all_seeds = set(all_seeds)
    num_target_voxels = N.max( [ int(x[2]) for x in lines ])


    # Find threshold value
    if threshold >= 1:
        raise Exception("Threshold must be 0<t<1")
    if samples==0:
        raise Exception("Samples must be >0")


    print " > Creating sparse matrix of dimensions %d x %d" % \
        (num_seed_voxels, num_target_voxels)

    if has_max_sv:
        num_seed_voxels = N.min((num_seed_voxels, max_seed_voxels))
    if has_max_tv:
        num_target_voxels = N.min((num_target_voxels, max_target_voxels))

    mat = [S.dok_matrix((1,num_target_voxels), int) for x in xrange(num_seed_voxels)]
    num_total = 0
    num_skipped_threshold = 0
    for l in lines:
        try:
            num_total += 1
            sv, tv, v = int(l[0]), int(l[2]), float(l[4])
            if has_max_tv and tv > num_target_voxels: break
            if has_max_sv and sv >= num_seed_voxels: continue

            if sv in all_seeds:
                all_seeds.remove(sv)

            if v<=0:
	       print "  Warning. Target voxel with value of %0.2f!    Seed:%d   Target:%d" % (v, sv, tv)

            # Make sure P(connected) is above threshold
            if v/samples < threshold:
                num_skipped_threshold += 1
                continue

            mat[sv-1][0,tv-1] = 1 if binarize else v

        except Exception as e:
            print e
            print sv, tv, v
            return

    if len(all_seeds)>0:
        print "  Warning the following seed voxels were absent from the sparse matrix file:"
        print all_seeds

    print "    Skipped %d/%d due to thresholding" % (num_skipped_threshold, num_total)
    return mat


def cross_correlate_matrix_binary(mat):
    # Create set-representation of matrix
    mat2 = [ set([y[1] for y in x.keys()]) for x in mat ]

    # Cache the norm of all vectors
    print " > Computing norms..."
    # because the value of each cell is always going to be 1, we can simplify
    norms = [N.sqrt(len(mat2[i])) for i in xrange(len(mat2))]

    # Calculate cross correlation
    cc_mat = N.zeros((len(mat), len(mat)))

    num_norm_zero = 0
    for i in xrange(len(mat)):
        print "%d/%d at %s" % (i, len(mat), time.ctime())
        num_nan = 0
        for j in xrange(len(mat)):
            if cc_mat[i,j] != 0: continue
            
            norm_mult = norms[i] * norms[j]
            if norm_mult == 0:
                num_norm_zero += 1
                continue

            cc_mat[i,j] = len(mat2[i].intersection(mat2[j])) / norm_mult
            if math.isnan(cc_mat[i,j]):
                cc_mat[i,j] = 0
                num_nan += 1
            cc_mat[j,i] = cc_mat[i,j]

        if num_nan>0:
            print "WARNING!! %d NaN values in the last line!" % num_nan

    print "Skipped %d cells due to 0 norm" % num_norm_zero
    return cc_mat



def cross_correlate_matrix_nonbinary(mat):
    # Cache the norm of all vectors
    print " > Computing norms..."
    # because the value of each cell is always going to be 1, we can simplify
    norms = [ N.sqrt(mat[i].multiply(mat[i]).sum(1)) for i in xrange(len(mat)) ]
    # Calculate cross correlation
    cc_mat = N.zeros((len(mat), len(mat)))

    print " > Computing cross-correlation matrix..."
    for i in xrange(len(mat)):
        print "%d/%d at %s" % (i, len(mat), time.ctime())
        num_nan = 0
        for j in xrange(len(mat)):
            if cc_mat[i,j] != 0: continue
            cc_mat[i,j] = mat[i].dot(mat[j].T)[0,0] / (norms[i] * norms[j])
            if math.isnan(cc_mat[i,j]):
                cc_mat[i,j] = 0
                num_nan += 1
            cc_mat[j,i] = cc_mat[i,j]

        if num_nan>0:
            print "WARNING!! %d NaN values in the last line!" % num_nan

    return cc_mat


def get_ccmat(options, args):
    cache_filename = corr_matrix_cache_filename(options.conn_mat_filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels, options.binarize)
    
    if not os.path.exists(cache_filename):
        mat = sparse_matrix_from_dotfile(options.conn_mat_filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels, options.binarize)
        if options.binarize:
            cc_mat = cross_correlate_matrix_binary(mat)
        else:
            cc_mat = cross_correlate_matrix_nonbinary(mat)

        # Saving cached version
        print " > Caching cross-correlation matrix to: %s" % cache_filename
        N.save(cache_filename, cc_mat)

    else:
        print " > Loading from cache: %s" % cache_filename
        cc_mat = N.load(cache_filename)

    if options.euclidian_const_scale != 0:
        print "     Max val in cross-corr matrix: %0.4f" % N.max(cc_mat)
        print "    Mean val in cross-corr matrix: %0.4f" % N.mean(cc_mat)
        print " > Applying euclidan distance matrix with scale parameter: %0.2f" % options.euclidian_const_scale
        dist_mat = get_seed_voxel_euclidian_dist_mat(options, args)
        dist_mat *= options.euclidian_const_scale
        print "     Max distance (post scale): %0.4f" % N.max(dist_mat)
        print "     Avg distance (post scale): %0.4f" % N.mean(dist_mat)
        cc_mat = cc_mat + dist_mat

    return cc_mat


def get_seed_voxel_euclidian_dist_mat(options, args):
    print " > Computing euclidian distance matrix between seed voxels..."
    # Load coords for all seed voxels
    seed_coords = []
    for line in open(options.voxel_coords_filename):
        pline = [x for x in line.strip().split(" ") if len(x)>0]
        seed_coords.append(N.asarray([int(pline[0]), int(pline[1]), int(pline[2])]))

    dist_mat = N.zeros((len(seed_coords), len(seed_coords))) - 1
    for i in xrange(len(dist_mat)):
        for j in xrange(len(dist_mat)):
            if dist_mat[i,j] != -1: continue
            dist = euclidean(seed_coords[i], seed_coords[j])
            dist_mat[i,j] = dist
            dist_mat[j,i] = dist

    dist_mat /= N.max(dist_mat)

    print "     Max distance: %0.4f" % N.max(dist_mat)
    print "     Avg distance: %0.4f" % N.mean(dist_mat)
    print "    Done euclidian distance mat."
    return dist_mat


def get_option_parser_group(parser):

    group = OptionGroup(parser, "Connection & correlation matrices options")

    group.add_option("", "--conn-mat-file", dest="conn_mat_filename",
                              help="connectivity matrix filename")
    group.add_option("", "--voxel-coords-file", dest="voxel_coords_filename",
                              help="voxel coordinates filename")
    group.add_option("-s", "--samples", dest="samples", default=0, type="int",
                              help="number of samples (ie number of streamlines started for each voxel)")
    group.add_option("-t", "--threshold", dest="threshold", default=0, type="float",
                              help="Consider tracts connected if P(connected) > THRESHOLD %, where P(connected)=matrix_value/samples")

    group.add_option("", "--no-binarize", dest="binarize", default=True, action="store_false",
                              help="Do not binarise input sparse matrix")
    group.add_option("", "--euclidian-const-scale", dest="euclidian_const_scale", default=0, type="float",
                              help="Add distance constraint to cross corr matrix by adding scaled euclidian matrix before kmeans. This parameter is multiplied to euclidian dist matrix before adding it. 0=off.")
    group.add_option("", "--max-seed-voxels", dest="max_seed_voxels", default=0, type="int",
                              help="Limit the number of seed voxels used to the first N")
    group.add_option("", "--max-target-voxels", dest="max_target_voxels", default=0, type="int",
                              help="Limit the number of target voxels used to the first N")
    group.add_option("", "--display-corr-matrix", action="store_true", dest="display_corr_matrix",
                              help="Display the correlation matrix")
    return group


