from __future__ import print_function
import random
import re
from Bio import pairwise2
from multiprocessing import Pool
import math
from itertools import repeat
from Bio.Seq import Seq
import numpy as np
import csv


def calculate_alignment(sequences, i):
    print("seq " + str(i))
    seq1 = sequences[i]
    scores = list(repeat(0, len(sequences)))
    max_penalty = (len(sequences) * len(sequences))/2 - len(sequences)
    for j, seq2 in enumerate(sequences):
        if j > i:
            score = pairwise2.align.globalxx(seq1, seq2)[0][2]
            score_revcomp = pairwise2.align.globalxx(str(Seq(seq1).reverse_complement()), seq2)[0][2]
            max_score = max(score, score_revcomp)
            if max_score < 7:
                scores[j] = 7-max_score
            else:
                scores[j] = -1 * max_penalty
    return scores


def f(m, v):
    s = 0
    for i in xrange(len(v)):
        for j in xrange(i + 1, len(v)):
            s = s + m[v[i]][v[j]]
    return s


def evaluate(cur_score, new_score, t):
    if new_score >= cur_score:
        return True
    elif t > 0 and math.exp((-(cur_score - new_score)) / t) > random.uniform(0, 1):
        return True
    else:
        return False



def simulated_annealing(m):
    v = random.sample(xrange(len(m)), 100)
    cur_score = f(m, v)
    scores = []
    cooling_iterations = 300
    t = 10000.0
    decr = t / cooling_iterations
    no_change = 0
    i = 0

    while no_change < 500:

        if t > decr:
            t -= decr
        else:
            t = 0

        p = random.uniform(0, 1)
        scores.append(cur_score)
        it_score = cur_score

        j = random.randint(0, len(v) - 1)
        if p <= 0.1:
            v_new = list(v)
            v_new.pop(j)
            if len(v) > 2:
                new_score = f(m, v_new)
                if evaluate(cur_score, new_score, t):
                    v = v_new
                    cur_score = new_score

        else:
            indices = [x for x in xrange(len(m)) if x not in v]
            random.shuffle(indices)
            for k in indices:
                v_new = list(v)

                if p > 0.4:
                    v_new[j] = k
                else:
                    v_new.append(k)

                new_score = f(m, v_new)
                if evaluate(cur_score, new_score, t):
                    v = v_new
                    cur_score = new_score

        if cur_score > it_score:
            print(no_change, it_score, cur_score, t)
            no_change = 0

        else:
            no_change += 1

        i += 1

    outfile = open('scores_new.txt', 'w')
    for score in scores:
        outfile.write(str(score) + '\n')
    outfile.close()

    return v


class SequenceFilters():
    @staticmethod
    def gc_content(sequences, num_gc):
        return [act_seq for act_seq in sequences if
                len(re.findall("c|g", act_seq, flags=re.IGNORECASE)) == num_gc]

    @staticmethod
    def repeats(sequences, num_repeats, length_repeats):
        new_sequences = []

        for seq in sequences:
            has_repeat = False
            for l in range(length_repeats):
                m = re.search("(([gatc]{" + str(l + 1) + "})\\2{" + str(num_repeats) + "})", seq, flags=re.IGNORECASE)
                if m is not None:
                    has_repeat = True
            if has_repeat is not True:
                new_sequences.append(seq)

        return new_sequences


    @staticmethod
    def similarity(sequences):
        sequences = random.sample(sequences, 1000)
        score_matrix = []

        pool = Pool(processes=5)
        results = []
        for i, seq1 in enumerate(sequences):
            results.append(pool.apply_async(calculate_alignment, args=(sequences, i,)))

        pool.close()
        pool.join()

        for res in results:
            score_matrix.append(res.get())

        for i in xrange(len(score_matrix)):
            for j in xrange(i + 1, len(score_matrix)):
                score_matrix[j][i] = score_matrix[i][j]

        # with open("score_matrix_new.txt") as tsv:
        #     reader = csv.reader(tsv, delimiter="\t")
        #     reader.next()
        #     for line in reader:
        #         del line[0]
        #         del line[-1]
        #         score_matrix.append(map(int, line))
        #

        best = open('score_matrix_new.txt', 'w')
        for i in xrange(0, len(sequences)):
            best.write(sequences[i] + '\t')
        best.write('\n')
        for i in xrange(0, len(sequences)):
            best.write(sequences[i] + '\t')
            for j in xrange(0, len(sequences)):
                best.write(str(score_matrix[i][j]) + '\t')
            best.write('\n')
        best.close()


        # G = nx.gnm_random_graph(100, 800)
        # score_matrix = nx.adjacency_matrix(G)
        v = simulated_annealing(score_matrix)

        best = open('best_new.txt', 'w')
        for i in v:
            best.write(sequences[i] + '\t')
        best.write('\n')
        for i in v:
            best.write(sequences[i] + '\t')
            for j in v:
                ind_i = sorted([i, j])
                best.write(str(score_matrix[ind_i[0]][ind_i[1]]) + '\t')
            best.write('\n')
        best.close()

        # G = nx.Graph()
        # for i, scores in enumerate(score_matrix):
        # for j, score in enumerate(scores):
        #         if score < 6.0:
        #             G.add_edge(i, j + i + 1)
        #
        # # cliques = clique.max_clique(G)
        # # print([sequences[i] for i in cliques])
        # p = MCP(G)
        # r = p.solve('glpk')
        # print(r.ff, r.solution)

