
import os, math
from optparse import OptionParser
import scipy.sparse as S
import numpy as N


def corr_matrix_cache_filename(filename, samples=0, threshold=0, max_seed_voxels=0, max_target_voxels=0):
    return "%s_%d_%0.4f_%d_%d.npy" % (filename, samples, threshold, max_seed_voxels, max_target_voxels)


def sparse_matrix_from_dotfile(filename, samples, threshold, max_seed_voxels=0, max_target_voxels=0):
    """
    threshold
    """
    
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
    norms = [N.sqrt(mat[i].multiply(mat[i]).sum(1)) for i in xrange(len(mat))]
    if any((math.isnan(x) for x in norms)): raise Exception("No norm can be NaN!")

    # Calculate cross correlation
    cc_mat = N.zeros((len(mat), len(mat)))

    print " > Computing cross-correlation matrix..."
    import time
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



def run(options, args):
    cache_filename = corr_matrix_cache_filename(options.filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels)
    
    if not os.path.exists(cache_filename):
        mat = sparse_matrix_from_dotfile(options.filename, options.samples, options.threshold, options.max_seed_voxels, options.max_target_voxels)
        cc_mat = cross_correlate_matrix(mat)

        # Saving cached version
        print " > Caching cross-correlation matrix to: %s" % cache_filename
        N.save(cache_filename, cc_mat)

    print "wouhou!"



if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                              help="connectivity matrix FILE")
    parser.add_option("-s", "--samples", dest="samples", default=0, type="int",
                              help="number of samples (ie number of streamlines started for each voxel)")
    parser.add_option("-t", "--threshold", dest="threshold", default=0, type="float",
                              help="Consider tracts connected if P(connected) > THRESHOLD %, where P(connected)=matrix_value/samples")
    
    parser.add_option("", "--max-seed-voxels", dest="max_seed_voxels", default=0, type="int",
                              help="Limit the number of seed voxels used to the first N")
    parser.add_option("", "--max-target-voxels", dest="max_target_voxels", default=0, type="int",
                              help="Limit the number of target voxels used to the first N")

    (options, args) = parser.parse_args()

    if not options.filename:
        raise Exception("Must specify filename. Run with --help for instructions.")

    run(options, args)



