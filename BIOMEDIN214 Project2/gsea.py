import sys
from math import sqrt
import pandas as pd
from numpy import cumsum
from random import shuffle


class GSEA:
    def __init__(self):
        self.expression = pd.DataFrame()
        self.sample = pd.DataFrame()
        self.kegg = {}

    def load_data(self, express_file, sample_file, kegg_file):
        """
        Input:
        :param express_file: expression
        :param sample_file: sample
        :param kegg_file: kegg
        :return: read the data and save to data frame/dict, return none.
        """
        f_e = open(express_file)
        f_s = open(sample_file)
        f_k = open(kegg_file)

        self.expression = pd.read_csv(f_e, sep='\t', index_col='SYMBOL')
        self.sample = pd.read_csv(f_s, sep='\t', header=None, names=['sample', 'label'])
        for line in f_k:
            geneset = line.strip().split()
            self.kegg[geneset[0]] = geneset[2:]

        f_e.close()
        f_s.close()
        f_k.close()

        pass

    def get_gene_rank_order(self):
        """
        calculate fc (take average), get ranked gene list
        :return: ranked gene list
        """
        ctrl_name = self.sample[(self.sample['label'] == 0)]['sample'].to_list()
        pt_name = self.sample[(self.sample['label'] == 1)]['sample'].to_list()

        # the values have already been log normalized
        ctrl_fc = self.expression[ctrl_name].mean(1)
        pt_fc = self.expression[pt_name].mean(1)

        # Calculate
        rank_fc = (pt_fc - ctrl_fc).sort_values(ascending=False)
        ranked_gene = rank_fc.index.to_list()

        return ranked_gene

        pass

    def score_parameters(self, geneset):
        """
        calculate points for hit and penalty for miss, called by random walk
        :param geneset: a string, name from kegg
        :return: point for hit and penalty for miss
        """
        total_set = self.expression.index.to_list()
        gene_set = self.kegg[geneset]
        filtered_gene_set = list(set(gene_set) - (set(gene_set) - set(total_set)))

        # calculate parameters
        # number of total genes and number of genes in the gene set
        # Attention!! filter out genes in the gene set that are not in the expression data
        n_total = len(total_set)
        n_gene_set = len(filtered_gene_set)

        # points you will get for hit and miss
        p_hit = sqrt((n_total - n_gene_set) / float(n_gene_set))
        p_miss = -sqrt(n_gene_set / float(n_total - n_gene_set))

        return p_hit, p_miss

        pass

    def random_walk(self, geneset, ranked_gene):
        """
        ES helper, called by get_enrichment_score
        :param geneset: a string, from kegg
        :param ranked_gene: a list of ranked gene
        :return: a list of points
        """
        gene_set = self.kegg[geneset]
        p_hit, p_miss = self.score_parameters(geneset)
        # ranked_gene = self.get_gene_rank_order()
        points = []

        for gene in ranked_gene:
            if gene in gene_set:
                points.append(p_hit)
            else:
                points.append(p_miss)

        return points

        pass

    def get_enrichment_score(self, geneset):
        """
        Calculate enrichment score (ES), should be a float correct to two decimal places for a given gene set
        :param geneset: the string for the gene set name corresponding to the gene set
                        (such as KEGG_CITRATE_CYCLE_TCA_CYCLE)
        :return: enrichment score, keep 2 decimal point
        """
        ranked_gene = self.get_gene_rank_order()
        points = self.random_walk(geneset, ranked_gene)
        scores = cumsum(points)
        enrichment_score = round(max(scores), 2)

        return enrichment_score

        pass

    def permuting_labels(self):
        """
        permuting sample labels - patients/controls
        :return: two list, one is patients, one is controls
        """
        # ctrl_count = len(self.sample[(self.sample['label'] == 0)]['sample'].to_list())
        # pt_count = len(self.sample[(self.sample['label'] == 1)]['sample'].to_list())
        label_list = self.sample['label'].to_list()
        sample_names = self.sample['sample'].to_list()
        permuting_sample_df = pd.DataFrame(columns=label_list, index=list(range(100)))
        # permuting_labels_dict = {}
        for i in range(100):
            shuffle(sample_names)
            permuting_sample_df.loc[i, :] = sample_names
            # permuting_labels_dict[i] = label_list
        return permuting_sample_df['0'], permuting_sample_df['1']

        pass

    def perm_enrichment_score(self, geneset):
        """
        iterate 100 times, to calculate ES for random labeled samples
        :param geneset: string
        :return: a list of ES
        """
        permuting_ctrl, permuting_pt = self.permuting_labels()
        perm_enrichment_scores = []
        for i in range(100):
            permuting_ctrl_name = permuting_ctrl.loc[i, :].to_list()
            permuting_pt_name = permuting_pt.loc[i, :].to_list()
            permuting_ctrl_fc = self.expression[permuting_ctrl_name].mean(1)
            permuting_pt_fc = self.expression[permuting_pt_name].mean(1)

            # Rank the genes
            perm_rank_fc = (permuting_pt_fc - permuting_ctrl_fc).sort_values(ascending=False)
            perm_ranks = perm_rank_fc.index.to_list()

            # Enrichment score
            perm_points = self.random_walk(geneset, perm_ranks)
            perm_scores = cumsum(perm_points)

            # Perm_ES is a list size = 100
            perm_enrichment_scores.append(round(max(perm_scores), 2))
            # print('round: ' + str(i))

        return perm_enrichment_scores

        pass

    def calculate_p_val(self, geneset):
        """
        calculate a p-value for gene set by counting the number of times in the permuted iterations
        :param geneset: a string
        :return: p-value for the geneset
        """
        perm_ES_list = self.perm_enrichment_score(geneset)
        actual_ES = self.get_enrichment_score(geneset)

        greater_count = 0
        for ES in perm_ES_list:
            if ES >= actual_ES:
                greater_count += 1

        p_val = greater_count / len(perm_ES_list)

        return p_val

    def get_sig_sets(self, p):
        """
        Correct p-values for number of tests (equal to the number of gene sets) via Bonferroni
        and then find the significant gene sets
        :param p: threshold
        :return: significant gene sets
        """
        kegg_list = list(self.kegg.keys())
        correct_threshold = p / len(self.kegg)
        sig_sets = []

        i = 0
        for geneset in kegg_list:
            i += 1
            p_val = self.calculate_p_val(geneset)
            if p_val < correct_threshold:
                sig_sets.append(geneset)
            # print("Compute set " + geneset + " (" + str(i) + '/' + str(len(kegg_list)) + ")")

        return sig_sets


def main():
    # Check that arguments are in present in the command line as expected

    if len(sys.argv) != 4:
        print("Please specify an expression file, an sample file and a geneset file as args.")
        return

    express_file = sys.argv[1]
    sample_file = sys.argv[2]
    kegg_file = sys.argv[3]

    """
    # input variables
    express_file = 'GSE25628_filtered_expression.txt'
    sample_file = 'GSE25628_samples.txt'
    kegg_file = 'c2.cp.kegg.v6.2.symbols.filtered.gmt'
    """
    p = 0.05

    gsea = GSEA()
    gsea.load_data(express_file, sample_file, kegg_file)
    gsea.get_gene_rank_order()
    # sig_genesets = gsea.get_sig_sets(p)
    gsea.get_sig_sets(p)
    # print(sig_genesets)


if __name__ == "__main__":
    main()
