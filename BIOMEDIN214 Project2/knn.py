import sys
from math import sqrt
import pandas as pd


# K-Nearest Neighbors classification


class KNN:
    def __init__(self):
        self.expression = {}
        self.sample = []
        self.sample_dict = {}
        self.assignments = {}
        self.assignment_list = []

    def load_data(self, express_file, sample_file):
        """
        read_file: read the expression data file - tab-delimited
        Inputs:
            Expression file: a row for each gene, a column for each sample
            Sample file: a row for each sample,
                        second column is the sample label, with a 0 - control and a 1 - patient.
        :return: none
        """
        f_e = open(express_file)
        f_s = open(sample_file)
        expression_df = pd.read_csv(f_e, sep='\t').drop(['SYMBOL'], axis=1)
        for column in expression_df.columns:
            self.expression[column] = expression_df[column].tolist()
        sample_df = pd.read_csv(f_s, sep='\t', header=None, names=['sample', 'label'])
        self.sample = list(zip(list(sample_df['sample']), list(sample_df['label'])))
        self.sample_dict = sample_df.set_index('sample')['label'].to_dict()

        f_e.close()

        f_s.close()

        pass

    def euclidean_distance(self, train, test):
        """
        Euclidean distance, the smaller the value, the more similar two records will be.
        :return: distance between two vectors
        """
        distance = 0.0
        for i in range(0, len(train)):
            distance += (train[i] - test[i]) ** 2
        return sqrt(distance)

        pass

    def get_assignments(self, k, fn):
        """
        Input:
        :param k: number of neighbors
        :param fn: between [0, 1]
        :return: should return the class assignments for all samples
        """

        # Leave one out
        for i in range(len(self.sample)):
            test_sample_name = self.sample[i][0]
            test = self.expression[test_sample_name]
            distance_dic = {}

            for j in range(len(self.sample)):
                if j == i:
                    continue
                train_sample_name = self.sample[j][0]
                distance = self.euclidean_distance(self.expression[train_sample_name], test)
                distance_dic[train_sample_name] = distance

            sorted_distance = [(k, v) for k, v in sorted(distance_dic.items(), key=lambda item: item[1])]
            positives = 0
            for j in range(k):
                sample_name = sorted_distance[j][0]
                positives += self.sample_dict[sample_name]
            predict_value = 0
            if 1.0 * positives / k > fn:
                predict_value = 1

            self.assignments[test_sample_name] = predict_value
            self.assignment_list.append(predict_value)
        # print(self.assignments)
        return self.assignment_list

        pass

    def calc_metrics(self, k, fn):
        """
        Input:
        :param k: number of neighbors
        :param fn: between [0, 1]
        :return: a list of float values [sensitivity,specificity] of a KNN classifier
        """
        truth = self.sample_dict
        _ = self.get_assignments(k, fn)
        pred = self.assignments
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        for sample_name in truth:
            if truth[sample_name] == 1 and pred[sample_name] == 1:
                TP += 1
            elif truth[sample_name] == 0 and pred[sample_name] == 1:
                FP += 1
            elif truth[sample_name] == 1 and pred[sample_name] == 0:
                FN += 1
            elif truth[sample_name] == 0 and pred[sample_name] == 0:
                TN += 1

        sensitivity = TP / (TP + FN)
        specificity = TN / (TN + FP)
        # accuracy = (TP + TN)/(TP+FN+TN+FP)

        return sensitivity, specificity

        pass


def main():
    # Check that arguments are in present in the command line as expected

    if len(sys.argv) != 3:
        print("Please specify an expression file and an sample file as args.")
        return

    express_file = sys.argv[1]
    sample_file = sys.argv[2]
    """
    # input variables
    express_file = 'GSE25628_filtered_expression.txt'
    sample_file = 'GSE25628_samples.txt'
    """

    k = 5
    fn = 0.5
    # create an align object and run
    knn = KNN()
    knn.load_data(express_file, sample_file)
    knn.get_assignments(k, fn)
    knn.calc_metrics(k, fn)


if __name__ == "__main__":
    main()
