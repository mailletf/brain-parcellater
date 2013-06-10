
import numpy as N
import scipy.spatial.distance as D
import scipy.sparse as S
import scipy.stats as SS
import conn_ccorr_mat as CCM

import unittest

size = 4

class CrossCorrelationTest(unittest.TestCase):
    def setUp(self):
        def computeCrossCorrMat(m):
            cm = N.zeros(m.shape)
            for i in xrange(size):
                for j in xrange(size):
                    #cm[i,j] = N.correlate(m[i], m[j])
                    cm[i,j] = 1 - D.cosine(m[i], m[j])
                    #cm[i,j] = SS.pearsonr(m[i], m[j])[0]
            return cm

        self.mat = N.random.random((size,size))
        self.cc_mat = computeCrossCorrMat(self.mat)

        self.bin_mat = N.round(self.mat)
        self.cc_bin_mat = computeCrossCorrMat(self.bin_mat)

    def test_nonBinaryTest(self):
        test_mat = [S.dok_matrix((1,size), int) for x in xrange(size)]
        for i in xrange(size):
            for j in xrange(size):
                test_mat[i][0,j] = self.mat[i,j]

        cc_test_mat = CCM.cross_correlate_matrix_nonbinary(test_mat)
        N.testing.assert_almost_equal(cc_test_mat, self.cc_mat)


    def test_binaryTest(self):
        test_mat = [S.dok_matrix((1,size), int) for x in xrange(size)]
        for i in xrange(size):
            for j in xrange(size):
                if self.mat[i,j]>0:
                    test_mat[i][0,j] = self.bin_mat[i,j]

        cc_test_mat = CCM.cross_correlate_matrix_binary(test_mat)
        N.testing.assert_almost_equal(cc_test_mat, self.cc_bin_mat)

        

if __name__ == '__main__':
    unittest.main()

