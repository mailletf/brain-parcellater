#!/usr/bin/env python
#
# gm_parcel.py
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
from optparse import OptionParser
from sklearn.cluster import KMeans
import scipy.sparse as S
import pylab
import numpy as N


def corr_matrix_cache_filename(filename, samples=0, threshold=0, max_seed_voxels=0, max_target_voxels=0):
    return "%s_%d_%0.4f_%d_%d.npy" % (filename, samples, threshold, max_seed_voxels, max_target_voxels)


def sparse_matrix_from_dotfile(filename, samples, threshold, max_seed_voxels=0, max_target_voxels=0):
    has_max_sv = max_seed_voxels > 0
    has_max_tv = max_target_voxels > 0

    # Load all lines
    lines = [x.strip().split(" ") for x in open(filename)]

    # Find max values for target and seed voxels to determine
    # how big the dense matrix would be
    num_target_voxels = N.max( [ int(x[2]) for x in lines ])
    num_seed_voxels = N.max( [ int(x[0]) for x in lines ])

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

            # Make sure P(connected) is above threshold
            if v/samples < threshold:
                num_skipped_threshold += 1
                continue

            # Set all values to one
            mat[sv-1][0,tv-1] = 1

        except Exception as e:
            print e
            print sv, tv, v
            return

    print "    Skipped %d/%d due to thresholding" % (num_skipped_threshold, num_total)
    return mat


def cross_correlate_matrix(mat):
    # Cache the norm of all vectors
    print " > Computing norms..."
    # because the value of each cell is always going to be 1, we can simplify
    # sqrt(mat[i].multiply(mat[i]).sum(1)) = sqrt(len(mat[i].nonzero()))
    norms = [N.sqrt(len(mat[i].nonzero())) for i in xrange(len(mat))]

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


def do_kmeans(cc_mat, num_clusters, keep_top=1):
    kmeans = KMeans(init='k-means++', n_clusters=num_clusters, n_init=25, precompute_distances=True, max_iter=2500)
    fitted = kmeans.fit(cc_mat)
    dist_from_clusters = fitted.transform(cc_mat)
    memberships = fitted.predict(cc_mat)

    # For each cluster, keep only top KEEP_TOP elements
    for cid in xrange(num_clusters):
        num_thresholded = 0
        c_members = N.where(memberships==cid)[0]
        
        c_distances = N.sort( [dist_from_clusters[i][cid] for i in c_members] )
        num_to_keep = N.floor(keep_top * len(c_members))
        threshold_value = c_distances[-num_to_keep]

        for i in c_members:
            if dist_from_clusters[i][cid] < threshold_value:
                num_thresholded += 1
                memberships[i] = -1

        print "   For cluster %d, %d/%d values were removed due to thresholding" % (cid, num_thresholded, len(c_members))

    return memberships


def write_to_ascii(voxel_coords_filename, clusters, output_ascii_filename, dimx, dimy, dimz):
    # load up coord file
    coords = {}
    for line in open(voxel_coords_filename):
        pline = [x for x in line.strip().split(" ") if len(x)>0]
        if int(pline[3])!=0: raise Exception("Only supports t=0")

        # check to see if we have thresholded out this voxel when doing
        # k-means. in that case, its value will be -1
        cluster_value = clusters[int(pline[4])-1]
        if cluster_value == -1: continue

        coords["%s-%s-%s" % (pline[0], pline[1], pline[2])] = str(cluster_value+1)


    # write output file
    writer = open(output_ascii_filename, "w")
    for z in xrange(dimz):
        for y in xrange(dimy):
            to_write = [coords["%s-%s-%s" % (x, y, z)] if "%s-%s-%s" % (x, y, z) in coords else "0" for x in xrange(dimx)]
            writer.write(" ".join(to_write)+"\n")
    writer.close()



def run(options, args):
    cache_filename = corr_matrix_cache_filename(options.conn_mat_filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels)
    
    if not os.path.exists(cache_filename):
        mat = sparse_matrix_from_dotfile(options.conn_mat_filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels)
        cc_mat = cross_correlate_matrix(mat)

        # Saving cached version
        print " > Caching cross-correlation matrix to: %s" % cache_filename
        N.save(cache_filename, cc_mat)

    else:
        cc_mat = N.load(cache_filename)

    if options.display_corr_matrix:
        pylab.imshow(cc_mat)
        pylab.show()
    
    else:
        print " > Running k-means..."
        if options.kmeans_keep_only <= 0 or options.kmeans_keep_only>1:
            raise Exception("kmeans-keep-only needs to be 0<x<=1")
        clusters = do_kmeans(cc_mat, options.num_clusters, options.kmeans_keep_only)
        print clusters
        print " > Writing ascii matrix file to: %s" % options.ascii_out_filename
        dimx = 91
        dimy = 109
        dimz = 91
        write_to_ascii(options.voxel_coords_filename, clusters, options.ascii_out_filename, dimx, dimy, dimz)



if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("", "--conn-mat-file", dest="conn_mat_filename",
                              help="connectivity matrix filename")
    parser.add_option("", "--voxel-coords-file", dest="voxel_coords_filename",
                              help="voxel coordinates filename")
    parser.add_option("", "--ascii-out-file", dest="ascii_out_filename",
                              help="ascii matrix output filename")
    parser.add_option("-c", "--clusters", dest="num_clusters", default=2, type="int",
                              help="number of k-means clusters")
    parser.add_option("", "--kmeans-keep-only", dest="kmeans_keep_only", default=1.0, type="float",
                              help="after having done k-means, only keep closest KMEAN-KEEP-ONLY % of cluster members")
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



