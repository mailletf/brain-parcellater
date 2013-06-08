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

