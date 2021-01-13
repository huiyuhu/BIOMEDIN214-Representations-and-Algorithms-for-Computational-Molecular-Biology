import unittest

from align import *

TEST_INPUT_FILE = "test_example.input"
TEST_OUTPUT_FILE = "test_example.output"


class TestAlignmentClasses(unittest.TestCase):

    def test_match_matrix(self):
        """
        Tests match matrix object
        """
        match_matrix = MatchMatrix()
        match_matrix.set_score("A", "C", 5)
        self.assertEqual(match_matrix.get_score("A", "C"), 5)

    def test_score_matrix_score(self):
        """
        Tests score matrix object score set + get methods
        """
        ### FILL IN ###
        # this should be very similar to test match matrix
        score_matrix = ScoreMatrix('M', 5, 5)
        score_matrix.set_score( 1, 1, 1.1)
        self.assertEqual(score_matrix.get_score(1, 1), 1.1)
        return

    def test_score_matrix_pointers(self):
        """
        Tests score matrix object pointer set + get methods
        """
        score_matrix = ScoreMatrix('M', 5, 5)
        score_matrix.set_pointers(1, 1, 'M', 2, 2)
        self.assertEqual(score_matrix.get_pointers(1, 1)[0], ('M', (2, 2)))
        return

    def test_param_loading(self):
        """
        Tests AlignmentParameters "load_params_from_file()" function
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)
        self.assertEqual(align_params.seq_a, "AATGC")
        self.assertEqual(align_params.seq_b, "AGGC")
        self.assertTrue(align_params.global_alignment)
        self.assertEqual(align_params.dx, 0.1)
        self.assertEqual(align_params.ex, 0.5)
        self.assertEqual(align_params.dy, 0.6)
        self.assertEqual(align_params.ey, 0.3)
        self.assertEqual(align_params.alphabet_a, "ATGC")
        self.assertEqual(align_params.alphabet_b, "ATGCX")
        self.assertEqual(align_params.len_alphabet_a, 4)
        self.assertEqual(align_params.len_alphabet_b, 5)

        # test that the match match is set up correctly
        #  if this fails, make sure you are loading the asymmetric matrix properly!
        match_mat = align_params.match_matrix
        self.assertEqual(match_mat.get_score("A", "X"), 0.3)
        self.assertEqual(match_mat.get_score("C", "G"), -0.3)
        self.assertEqual(match_mat.get_score("G", "C"), 0)

    def test_update_ix(self):
        """
        Test AlignmentAlgorithm's update Ix
        """

        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5

        # create an alignment object
        align = Align(TEST_INPUT_FILE, TEST_OUTPUT_FILE)
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.ix_matrix.set_score(2, 2, 2.5)

        # run the method!
        align.update_ix(3, 2)

        score = align.ix_matrix.get_score(3, 2)
        pointers = align.ix_matrix.get_pointers(3, 2)
        self.assertEqual(pointers, [('M', (2, 2)), ('Ix', (2, 2))])
        self.assertEqual(score, 2)

        return
        ### FILL IN for checking pointers! ###
        # note - in this example, it should point to M -AND- Ix
        # check by hand!

    def test_update_m(self):
        """
        Test AlignmentAlgorithm's update M
        """
        ### FILL IN ###
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5

        # create an alignment object
        align = Align(TEST_INPUT_FILE, TEST_OUTPUT_FILE)
        align_params.seq_a="ATGC"
        align_params.seq_b="ATGC"
        align.align_params = align_params
        matchM = MatchMatrix()
        matchM.set_score('A', 'A', 1)
        matchM.set_score('T', 'T', 1)
        matchM.set_score('G', 'G', 1)
        matchM.set_score('C', 'C', 1)
        align_params.match_matrix = matchM

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.ix_matrix.set_score(2, 2, 2.5)
        align.iy_matrix.set_score(2, 2, 2)

        # run the method!
        align.update_m(3, 3)

        score = align.m_matrix.get_score(3, 3)
        pointers = align.m_matrix.get_pointers(3, 3)
        self.assertEqual(pointers, [('M', (2, 2))])
        self.assertEqual(score, 4)
        return

    def test_update_iy(self):
        """
        Test AlignmentAlgorithm's update Iy
        """
        ### FILL IN ###
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dx = 1
        align_params.ex = 0.5

        # create an alignment object
        align = Align(TEST_INPUT_FILE, TEST_OUTPUT_FILE)
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.iy_matrix.set_score(2, 2, 2.5)

        # run the method!
        align.update_iy(2, 3)

        score = align.iy_matrix.get_score(2, 3)
        pointers = align.iy_matrix.get_pointers(2, 3)
        self.assertEqual(pointers, [('M', (2, 2)), ('Iy', (2, 2))])
        self.assertEqual(score, 2)

        return

    def test_traceback_start(self):
        """
        Tests that the traceback finds the correct start
        Should test local and global alignment!
        """
        align_params = AlignmentParameters()
        align_params.seq_a = "ACGT"
        align_params.seq_b = "ACGT"
        align_params.global_alignment = True

        align = Align(TEST_INPUT_FILE, TEST_OUTPUT_FILE)
        align.align_params = align_params
        align.m_matrix.set_score(4, 4, 1.0)
        align.ix_matrix.set_score(3, 4, 1.0)
        start = align.find_traceback_start()

        ## Test max score
        self.assertEqual(start[0], 1.0)
        ## Test max locs
        self.assertEqual(start[1], [('M', 4, 4), ('Ix', 3, 4)])
        return


if __name__ == '__main__':
    unittest.main()
