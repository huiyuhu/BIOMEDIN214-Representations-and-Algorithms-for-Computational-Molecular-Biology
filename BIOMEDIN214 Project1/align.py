import sys
from typing import Set, Tuple  # NOTE: You may need to "pip install typing" locally if this import gives you errors
import numpy as np


def fuzzy_equals(a: float, b: float) -> bool:
    """
    Checks if two floating point numbers are equivalent.
    :param a: first number
    :param b: second number
    :return: True if a and be are within epsilon = 10^-6
    """
    epsilon = 10 ** (-6)
    return abs(a - b) < epsilon


def removeDuplicates(lst):
    return list(set([i for i in lst]))


#### ------- CLASSES ------- ####
class MatchMatrix(object):
    """
    A class representation of the match matrix, S. Stores the scores of matches between characters.
    """

    def __init__(self):
        self.matchM = {}  # dictionary is used for store score in match matrix

        pass

    def set_score(self, a: str, b: str, score: float) -> None:
        """
        Updates or adds a score for a specified match

        :param a: the character from sequence A
        :param b: the character from sequence B
        :param score: the score to set the match M(a,b)
        """

        self.matchM[(a, b)] = score

        pass

    def get_score(self, a: str, b: str) -> float:
        """
        Returns the score for a particular match, where a is the
        character from sequence A and b is from sequence B.

        :param a: the character from sequence A
        :param b: the character from sequence B
        :return: the score of that match, M(a,b)
        """

        return self.matchM[(a, b)]

        pass


class ScoreEntry(object):
    # For ScoreMatrix Class

    def __init__(self, row, col, name=None):
        self.row = row
        self.col = col
        self.name = name
        self.score = 0.0
        self.pointers = []


class ScoreMatrix(object):
    # score of best alignment ends
    """
    A class representation of the score matrices (M, Ix, Iy), which will be dynamically updated.
    The score matrix consists of a 2-D array of ScoreEntries that are updated during alignment
    and used to output the maximum alignment.
    """

    def __init__(self, name: str, nrow: int, ncol: int) -> None:
        """
        Initialize ScoreMatrix class.

        :param name: identifier for the score matrix, should be in {Ix, Iy, M}
        :param nrow: number of rows for ScoreMatrix
        :param ncol: number of columns for ScoreMatrix
        """

        self.nrow = nrow
        self.ncol = ncol
        self.name = name

        # Create a score matrix, each element is an instance of ScoreEntry.
        self.score_matrix = [[ScoreEntry(i, j, name) for j in range(ncol)] for i in range(nrow)]

    def get_score(self, row: int, col: int) -> float:
        """
        Return the current score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :return: current score at position (row, col)
        """

        return self.score_matrix[row][col].score

        pass

    def set_score(self, row: int, col: int, score: float) -> None:
        """
        Set the score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :param score: score to set at position (row, col)
        """

        self.score_matrix[row][col].score = score

        pass

    def get_pointers(self, row: int, col: int):
        """
        Return the indices of the entries being pointed to. Remember, these are essential for the final traceback.

        :param row: row index
        :param col: column index
        :return: a set of indices (represented as tuples) corresponding to other entries being pointed to for traceback
        ex: {(1,0), (1,1)}
        """
        return self.score_matrix[row][col].pointers

        pass

    def set_pointers(self, row: int, col: int, pointer_name: str, pointer_row: int, pointer_col: int) -> None:
        """
        Add pointers to each entry in the score matrix.

        :param row: row index
        :param col: column index
        :param pointers: set of pointers to add to your entry
        """

        self.score_matrix[row][col].pointers.append((pointer_name, (pointer_row, pointer_col)))

        pass

    def print_scores(self) -> str:
        """
        Returns a nicely formatted string containing the scores in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        :return: a string representation of the scores in the score matrix
        """
        # Print scores line by line
        for i in range(len(self.score_matrix)):
            for j in range(len(self.score_matrix[i])):
                print(self.score_matrix[i][j].score, end=' ')
            print()

        pass

    def print_pointers(self) -> str:
        """
        Returns a nicely formatted string containing the pointers in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!
        """
        # Print pointers line by line
        for i in range(len(self.score_matrix)):
            for j in range(len(self.score_matrix[i])):
                print('[', self.get_pointers(i, j), end=' ]')
            print()

        pass


class AlignmentParameters(object):
    """
    A class to hold the alignment parameters.
    """

    def __init__(self) -> None:
        """
        Initialize AlignmentParameters object with default parameters.
        """

        # The definition for all of these class variables are documented in the annotated version of the input file
        # on the P1 page on Canvas.
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.len_alphabet_a = 0
        self.alphabet_a = ""
        self.len_alphabet_b = 0
        self.alphabet_b = ""
        self.match_matrix = MatchMatrix()


    def load_params_from_file(self, input_file: str) -> None:
        """
        Read the parameters from an input file and update the alignment parameters accordingly

        :param input_file: path to the alignment input file (whose structure is defined on the project page)

        """

        with open(input_file) as f:
            # sequence of letters indicating sequence A,B
            self.seq_a = f.readline().strip()
            self.seq_b = f.readline().strip()

            # An indication of whether local (1) or global (0) alignment is sought.
            if int(f.readline().strip()) == 0:
                self.global_alignment = True

            # A set of gap penalties for introducing gaps into A or B.
            gap_pen = f.readline().strip().split(' ')
            self.dx, self.ex, self.dy, self.ey = gap_pen
            self.dx = float(self.dx)
            self.ex = float(self.ex)
            self.dy = float(self.dy)
            self.ey = float(self.ey)

            self.len_alphabet_a = int(f.readline().strip())
            self.alphabet_a = f.readline().strip()
            self.len_alphabet_b = int(f.readline().strip())
            self.alphabet_b = f.readline().strip()

            # Match Matrix: score between an element in A and one in B
            raw = []  # a list for data input
            lines = f.readlines()

            # place scores in file in the match matrix
            for line in lines:
                raw.append(line.split())
            for i in raw:
                if len(i) != 5:
                    continue
                a = i[2]  # match matrix row
                b = i[3]  # match matrix col
                score = float(i[4])  # match matrix score (s)
                self.match_matrix.set_score(a, b, score)  # set score in match matrix

        pass


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by calling "align()"

    NOTE ON IMPLEMENTATION: You might find it helpful for code structure/organization to create additional functions
    that handle subtasks during the alignment. This is totally acceptable. However, you MUST implement the following
    functions with the expected behavior as outlined in the docstrings, which we will check in order to give
    at least partial credit for your P1 implementations:

    - populate_score_matrices
    - update_m, update_ix, and update_iy
    - find_traceback_start

    NOTE 2: Don't forget about that fuzzy_equals function at the top.

    """

    def __init__(self, input_file: str, output_file: str) -> None:
        """
        Initialize Align object.

        :param input_file: alignment input file path
        :param output_file: file path to write the output alignment
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()
        len_a = 0
        len_b = 0

        ### FILL IN ###
        # NOTE: be careful about how you initialize these!

        # Get params of the input file
        if input_file != "":
            self.align_params.load_params_from_file(self.input_file)

            len_a = len(self.align_params.seq_a) + 1
            len_b = len(self.align_params.seq_b) + 1

        # initiate the score matrix
        self.m_matrix = ScoreMatrix('M', len_a, len_b)
        self.ix_matrix = ScoreMatrix('Ix', len_a, len_b)
        self.iy_matrix = ScoreMatrix('Iy', len_a, len_b)
        self.match_dictionary = {}  # traceback optimization - a dictionary to save the result

    def align(self):
        """
        Main method for running the alignment.

        Note: there is no static typing on the method, as you can choose to return arbitrary
        intermediates/output if it's helpful for debugging. The essential minimal functionality that this
        method must have is to write the resulting alignments to the output file in the format specified in the
        project page
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        self.traceback()
        self.write_output()

    def populate_score_matrices(self) -> None:
        """
        Populate the score matrices based on the data in align_params. Should call update(i,j) for each entry
        in the score matrices.
        """
        len_a = len(self.align_params.seq_a)  #length of seq a
        len_b = len(self.align_params.seq_b)  #length of seq b

        for i in range(1, len_a + 1):
            for j in range(1, len_b + 1):
                self.update(i, j)  # call update() for each entry
        pass

    def update(self, row: int, col: int) -> None:
        """
        Update all matrices at a given row and column index.

        :param row: index of row to update
        :param col: index of column to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row: int, col: int) -> None:
        """
        Update matrix M.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###
        seq_a = self.align_params.seq_a
        seq_b = self.align_params.seq_b

        ## Match Matrix score S

        s = self.align_params.match_matrix

        ## Score Calculation

        m = self.m_matrix.get_score(row - 1, col - 1) + s.get_score(seq_a[row - 1], seq_b[col - 1])
        ix = self.ix_matrix.get_score(row - 1, col - 1) + s.get_score(seq_a[row - 1], seq_b[col - 1])
        iy = self.iy_matrix.get_score(row - 1, col - 1) + s.get_score(seq_a[row - 1], seq_b[col - 1])

        ## Find Max Score among M, Ix, Iy
        max_score = max(m, ix, iy)

        ## No neg score for local alignment - set to zero if neg
        if not self.align_params.global_alignment:
            self.m_matrix.set_score(row, col, max(max_score, 0))
        else:
            self.m_matrix.set_score(row, col, max_score)

        ## Pointer
        if fuzzy_equals(m, max_score):
            self.m_matrix.set_pointers(row, col, 'M', row - 1, col - 1)
        if fuzzy_equals(ix, max_score):
            self.m_matrix.set_pointers(row, col, 'Ix', row - 1, col - 1)
        if fuzzy_equals(iy, max_score):
            self.m_matrix.set_pointers(row, col, 'Iy', row - 1, col - 1)

        pass

    def update_ix(self, row: int, col: int) -> None:
        """
        Update matrix Ix.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###

        ## Score Calculation
        m = self.m_matrix.get_score(row - 1, col) - self.align_params.dy
        ix = self.ix_matrix.get_score(row - 1, col) - self.align_params.ey

        ## Find Max
        max_score = max(m, ix)

        ## No neg score for local alignment
        if not self.align_params.global_alignment:
            self.ix_matrix.set_score(row, col, max(max_score, 0))
        else:
            self.ix_matrix.set_score(row, col, max_score)

        ## Pointer
        if fuzzy_equals(m, max_score):
            self.ix_matrix.set_pointers(row, col, 'M', row - 1, col)
        if fuzzy_equals(ix, max_score):
            self.ix_matrix.set_pointers(row, col, 'Ix', row - 1, col)

        pass

    def update_iy(self, row: int, col: int) -> None:
        """
        Update matrix Iy.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###

        ## Score Calculation
        m = self.m_matrix.get_score(row, col - 1) - self.align_params.dx
        iy = self.iy_matrix.get_score(row, col - 1) - self.align_params.ex

        ## Find Max
        max_score = max(m, iy)

        ## No neg score for local alignment
        if not self.align_params.global_alignment:
            self.iy_matrix.set_score(row, col, max(max_score, 0))
        else:
            self.iy_matrix.set_score(row, col, max_score)

        ## Pointer
        if fuzzy_equals(m, max_score):
            self.iy_matrix.set_pointers(row, col, 'M', row, col - 1)
        if fuzzy_equals(iy, max_score):
            self.iy_matrix.set_pointers(row, col, 'Iy', row, col - 1)

        pass

    def find_traceback_start(self) -> Tuple[float, Set[Tuple[str, int, int]]]:
        """
        Find the location(s) to start the traceback and the corresponding best score.
        NOTE: Think carefully about how to set this up for local alignment.

        :return: The value of the best score and the location(s) to start the traceback to produce this score.
        The expected format is (max_val, max_loc), where max_val is the best score and max_loc is a set containing
        the positions to start the traceback that produce the best score, where each position is represented by a tuple
        [ex. (5.5, {('M',1,2), ('Ix',3,4)}) ].
        """
        ### FILL IN ###

        seq_a = self.align_params.seq_a
        seq_b = self.align_params.seq_b

        ## Global: traceback start from best score in the last row/column
        score_list_m = []
        score_list_ix = []
        score_list_iy = []

        if self.align_params.global_alignment:
            for i in range(1, len(seq_a) + 1):  # last column
                score_list_m = np.append(score_list_m, self.m_matrix.get_score(i, len(seq_b)))

            for j in range(1, len(seq_b) + 1):  # last row
                score_list_m = np.append(score_list_m, self.m_matrix.get_score(len(seq_a), j))

            for i in range(1, len(seq_a) + 1):  # last column
                score_list_ix = np.append(score_list_ix, self.ix_matrix.get_score(i, len(seq_b)))
            for j in range(1, len(seq_b) + 1):  # last row
                score_list_ix = np.append(score_list_ix, self.ix_matrix.get_score(len(seq_a), j))

            for i in range(1, len(seq_a) + 1):  # last column
                score_list_iy = np.append(score_list_iy, self.iy_matrix.get_score(i, len(seq_b)))
            for j in range(1, len(seq_b) + 1):  # last row
                score_list_iy = np.append(score_list_iy, self.iy_matrix.get_score(len(seq_a), j))

            list_mx = np.append(score_list_m, score_list_ix)
            list_mxy = np.append(list_mx, score_list_iy)
            max_val = max(list_mxy)

            max_loc = []
            for i in range(1, len(seq_a) + 1):  # last column
                if fuzzy_equals(self.m_matrix.get_score(i, len(seq_b)), max_val):
                    max_loc.append(('M', i, len(seq_b)))
                if fuzzy_equals(self.ix_matrix.get_score(i, len(seq_b)), max_val):
                    max_loc.append(('Ix', i, len(seq_b)))
                if fuzzy_equals(self.iy_matrix.get_score(i, len(seq_b)), max_val):
                    max_loc.append(('Iy', i, len(seq_b)))

            for j in range(1, len(seq_b) + 1):  # last row
                if fuzzy_equals(self.m_matrix.get_score(len(seq_a), j), max_val):
                    max_loc.append(('M', len(seq_a), j))
                if fuzzy_equals(self.ix_matrix.get_score(len(seq_a), j), max_val):
                    max_loc.append(('Ix', len(seq_a), j))
                if fuzzy_equals(self.iy_matrix.get_score(len(seq_a), j), max_val):
                    max_loc.append(('Iy', len(seq_a), j))

            max_loc = set(max_loc)

            return max_val, max_loc

        ## Local: traceback start from best score anywhare in the matrix
        else:
            max_value = 0.0
            for i in range(len(seq_a) + 1):
                for j in range(len(seq_b) + 1):
                    if self.m_matrix.get_score(i, j) > max_value:
                        max_value = self.m_matrix.get_score(i, j)
                    if self.ix_matrix.get_score(i, j) > max_value:
                        max_value = self.ix_matrix.get_score(i, j)
                    if self.iy_matrix.get_score(i, j) > max_value:
                        max_value = self.iy_matrix.get_score(i, j)
            max_loc = []
            for i in range(1, len(seq_a) + 1):
                for j in range(1, len(seq_b) + 1):
                    if fuzzy_equals(self.m_matrix.get_score(i, j), max_value):
                        max_loc.append(('M', i, j))
                    if fuzzy_equals(self.ix_matrix.get_score(i, j), max_value):
                        max_loc.append(('Ix', i, j))
                    if fuzzy_equals(self.iy_matrix.get_score(i, j), max_value):
                        max_loc.append(('Iy', i, j))

            max_val = max_value
            max_loc = set(max_loc)

            return max_val, max_loc

        pass

    def traceback_helper(self, loc):
        """
        Traceback Recursion: Find the traceback path between M, Ix, Iy
        Input: start point Tuple (matrix, (row, col))
        Output: A list of trackback path
        """

        results = []  # a list to save the result

        key = loc[0] + str(loc[1][0]) + str(loc[1][1])

        if key in self.match_dictionary:
            return self.match_dictionary[key]  # get the result if it is already in the dictionary

        # get the pointers and scores in the score matrices
        if loc[0] == 'M':
            pointers = self.m_matrix.get_pointers(loc[1][0], loc[1][1])
            score = self.m_matrix.get_score(loc[1][0], loc[1][1])
        elif loc[0] == 'Ix':
            pointers = self.ix_matrix.get_pointers(loc[1][0], loc[1][1])
            score = self.ix_matrix.get_score(loc[1][0], loc[1][1])
        elif loc[0] == 'Iy':
            pointers = self.iy_matrix.get_pointers(loc[1][0], loc[1][1])
            score = self.iy_matrix.get_score(loc[1][0], loc[1][1])

        # Local alignment - stop when score reach zero.
        if self.align_params.global_alignment:
            if loc[1][0] == 0 or loc[1][1] == 0:
                self.match_dictionary[key] = [[]]
                return [[]]
            if len(pointers) == 0:
                self.match_dictionary[key] = [[loc]]
                return [[loc]]
        else:
            if score == 0:
                self.match_dictionary[key] = [[]]
                return [[]]

        for pointer in pointers:
            results = results + self.traceback_helper(pointer)  # Recursion

        new_result = []

        for r in results:
            r = [loc] + r
            new_result.append(r)
        self.match_dictionary[key] = new_result
        return new_result

        pass

    def print_pathlist(self, path_list):
        """
        Print the traceback path list
        Input: a list of path, each path is a list of Tuple (matrix, (row, col))
        Output: A list of trackback path
        """
        # Path list have multiple path, print each pointers in each path
        for path in path_list:
            path_str = ''
            for pointer in path:
                path_str = path_str + pointer[0] + '(' + str(pointer[1][0]) + ',' + str(pointer[1][1]) + ')-->'
            print(path_str[0:-3])

    def traceback(self):  ### FILL IN additional arguments ###
        """
        Perform a traceback.

        NOTE: There is no static typing on the method, as you have freedom to choose how you'd like to implement
        this method, including which arguments to include and what to return.

        HINT: It is extremely helpful for debugging to include a way to print the traceback path.
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """
        path_list = []
        max_locs = self.find_traceback_start()

        # Use traceback_helper to get all the path
        for max_loc in max_locs[1]:
            # Convert the loc to a tuple that traceback_helper could use
            start = (max_loc[0], (max_loc[1], max_loc[2]))
            paths = self.traceback_helper(start)
            path_list = path_list + paths

        seq_a = self.align_params.seq_a
        seq_b = self.align_params.seq_b

        alignment_pairs = []

        # Get the alignment based on path list
        for path in path_list:
            alignment_a = ''
            alignment_b = ''
            for pointer in path:
                if pointer[0] == 'M':  # A[i] - B[j]
                    alignment_a = seq_a[pointer[1][0] - 1] + alignment_a
                    alignment_b = seq_b[pointer[1][1] - 1] + alignment_b
                if pointer[0] == 'Ix':  # A[i]:Gap
                    alignment_a = seq_a[pointer[1][0] - 1] + alignment_a
                    alignment_b = '_' + alignment_b
                if pointer[0] == 'Iy':  # Gap:B[j]
                    alignment_a = '_' + alignment_a
                    alignment_b = seq_b[pointer[1][1] - 1] + alignment_b

            alignment_pair = (alignment_a, alignment_b)  ## Combine seqA and seqB to a tuple
            alignment_pairs.append(alignment_pair)  ## save all the pairs

        return removeDuplicates(alignment_pairs)  ## remove duplicated alignments

        pass

    def write_output(self) -> None:
        """
        Write the output of an alignment to the output file.
        Output includes score and alignment result.
        """
        max_val = self.find_traceback_start()[0]
        alignment_pairs = self.traceback()

        ## Write Score
        out_file = open(self.output_file, "w")
        out_file.write(str(round(max_val, 1)))
        out_file.write('\n')

        ## Write Alignment
        for pair in alignment_pairs:
            seq_a = pair[0]
            seq_b = pair[1]
            out_file.write('\n')
            out_file.write(seq_a + '\n')
            out_file.write(seq_b + '\n')

        out_file.close()

        pass


# DO NOT MODIFY THE CODE BELOW!
def main():
    """
    Run the align function from command line, passing the input and output paths as arguments.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) != 3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__ == "__main__":
    main()

