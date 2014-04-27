# io_util.py
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

import os

def write_to_ascii(voxel_coords_filename, clusters, output_ascii_filename, dimx, dimy, dimz):
    # load up coord file
    coords = {}
    for line in open(voxel_coords_filename):
        pline = [x for x in line.strip().split(" ") if len(x)>0]
        if len(pline) != 5:
            print pline
            raise Exception("Expected split line to be len 5. Got %d. Line: %s" % (len(pline), line))
        if int(pline[3])!=0: raise Exception("Only supports t=0")

        # check to see if we have thresholded out this voxel when doing
        # k-means. in that case, its value will be -1
        cluster_idx = int(pline[4])-1
        if cluster_idx > len(clusters)-4:
            raise Exception("Cluster index out of bounds when writing out ascii file! Should you be running with surface output?")
        cluster_value = clusters[cluster_idx]
        if cluster_value == -1: continue

        coords["%s-%s-%s" % (pline[0], pline[1], pline[2])] = str(cluster_value+1)


    # write output file
    writer = open(output_ascii_filename, "w")
    for z in xrange(dimz):
        for y in xrange(dimy):
            to_write = [coords["%s-%s-%s" % (x, y, z)] if "%s-%s-%s" % (x, y, z) in coords else "0" for x in xrange(dimx)]
            writer.write(" ".join(to_write)+"\n")
    writer.close()

def write_to_surface_format(clusters, out_filename_prefix, label_filename):
    if not os.path.exists(label_filename):
        raise Exception("label_filename '%s' does not exist!" % label_filename)    

    label_lines = open(label_filename).readlines()

    assert(len(label_lines)-2 == len(clusters))

    # open writers
    writers = []
    for c_idx in xrange(max(clusters)+1):
        toutput_fn = "%s_c%d" % (out_filename_prefix, c_idx)
        print "  Opening writer: %s" % toutput_fn
        writers.append(open(toutput_fn, "w"))

    # writer headers
    for w in writers:
        for i in xrange(2):
            w.write(label_lines[i])

    # for each seed voxel, write in correct output file based on its cluster value
    for seed_voxel, line_to_write in enumerate(label_lines[2:]):
        cluster_val = clusters[seed_voxel]
        writers[cluster_val].write(line_to_write)


