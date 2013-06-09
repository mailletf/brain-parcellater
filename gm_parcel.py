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
import pylab
import numpy as N

import conn_ccorr_mat as CC
import io_util as io


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



def run(options, args):
    cc_mat = CC.get_ccmat(options, args)

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
        io.write_to_ascii(options.voxel_coords_filename, clusters, options.ascii_out_filename, dimx, dimy, dimz)


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

    parser.add_option("", "--no-binarize", dest="binarize", default=True, action="store_false",
                              help="Do not binarise input sparse matrix")
    parser.add_option("", "--max-seed-voxels", dest="max_seed_voxels", default=0, type="int",
                              help="Limit the number of seed voxels used to the first N")
    parser.add_option("", "--max-target-voxels", dest="max_target_voxels", default=0, type="int",
                              help="Limit the number of target voxels used to the first N")

    parser.add_option("", "--display-corr-matrix", action="store_true", dest="display_corr_matrix")

    (options, args) = parser.parse_args()

    if not options.conn_mat_filename:
        raise Exception("Must specify filename. Run with --help for instructions.")

    run(options, args)



