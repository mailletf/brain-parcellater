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
from optparse import OptionParser, OptionGroup
from sklearn.cluster import KMeans, spectral_clustering
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


def show_kmeans_clusters(cc_mat, clusters, num_clusters):
    #  itemindex=numpy.where(array==item)
    f, axarr = pylab.subplots(num_clusters, 1, sharex=True)
    for sp_idx, f in enumerate(axarr):
        #f = fig.add_subplot((num_clusters, 1, sb_idx+1), sharex=True)
        #f.set_title(text="Cross-correlation matrix")
        f.imshow(cc_mat[N.where(clusters==sp_idx)[0]], interpolation="nearest")

    pylab.show()


def run(options, args):
    cc_mat = CC.get_ccmat(options, args)

    if options.display_corr_matrix:
        pylab.imshow(cc_mat)
        pylab.show()
        return

    if options.clustering_algo == "kmeans":
        print " > Running k-means..."
        if options.kmeans_keep_only <= 0 or options.kmeans_keep_only>1:
            raise Exception("kmeans-keep-only needs to be 0<x<=1")
        clusters = do_kmeans(cc_mat, options.num_clusters, options.kmeans_keep_only)

    elif options.clustering_algo == "spectral":
        clusters = spectral_clustering(cc_mat, n_clusters=options.num_clusters,
            assign_labels="discretize")
            #eigen_solver='arpack', assign_labels="discretize")

    print clusters
    print " > Writing output file in '%s' format to: %s" % (options.out_format, options.out_filename)
    if options.out_format == "vol":
        io.write_to_ascii(options.voxel_coords_filename, clusters, options.out_filename,
            options.dimx, options.dimy, options.dimz)
    elif options.out_format == "surf":
        io.write_to_surface_format(clusters, options.out_filename, options.label_filename)

    if options.show_clusters:
        show_kmeans_clusters(cc_mat, clusters, options.num_clusters)


if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("", "--out-file", dest="out_filename",
                              help="output filename")
    parser.add_option("", "--label-file", dest="label_filename", default="",
                              help="label filename to be reordered when outputing results in surface format")
    parser.add_option("-o", "--out-format", dest="out_format", type="string",
                              default="vol", help="output file format: vol, surf")
    parser.add_option("-c", "--clusters", dest="num_clusters", default=2, type="int",
                              help="number of k-means clusters")

    parser.add_option("-a", "--algo", dest="clustering_algo", default="kmeans", type="choice",
                              help="clustering algorithm to use", choices=["kmeans", "spectral"])
    parser.add_option("", "--show-clusters", dest="show_clusters", default=False, action="store_true",
                              help="show cross-correlation matrix rows grouped by clusters")


    kmeans_group = OptionGroup(parser, "K-mean options")
    kmeans_group.add_option("", "--kmeans-keep-only", dest="kmeans_keep_only", default=1.0, type="float",
                              help="after having done k-means, only keep closest KMEAN-KEEP-ONLY % of cluster members")
    parser.add_option_group(kmeans_group)

    parser.add_option_group(CC.get_option_parser_group(parser))
    
    group = OptionGroup(parser, "Dimension of the coordinate space")
    group.add_option("-x", "--xdim", dest="dimx", default=91, type="int",
                              help="x size. default for 2mm")
    group.add_option("-y", "--ydim", dest="dimy", default=109, type="int",
                              help="y size. default for 2mm")
    group.add_option("-z", "--zdim", dest="dimz", default=91, type="int",
                              help="z size. default for 2mm")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if not options.out_format in ["vol", "surf"]:
        raise Exception("Invalid output type '%s'" % options.out_format)

    if not options.conn_mat_filename:
        raise Exception("Must specify filename. Run with --help for instructions.")

    run(options, args)



