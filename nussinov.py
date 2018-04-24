#!/usr/bin/env python 

'''Created by Tobias Landenberger on 23.04.18 1:06 PM'''
import argparse


# Initialize NxN DP Matrix with zeroes
def initializeMatrix(n):
    dp = [[0 for i in range(n)] for j in range(n)]
    return dp


def delta(i, j, seq, minLoopLength):
    if i >= j-minLoopLength:
        score = 0
    else:
        if (seq[i] == 'G' and seq[j] == 'C') or (seq[i] == 'C' and seq[j] == 'G'):
            score = int(args.score_GC)
        elif (seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A'):
            score = int(args.score_AU)
        elif (seq[i] == 'G' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'G'):
            score = int(args.score_GU)
        else:
            score = 0
    return score


def fourth_case(i, j, dp):
    tmp = 0
    for k in range(i + 1, j - 1):
        tmp = max(tmp, (dp[i][k] + dp[k + 1][j]))
    return tmp


def gamma(i, j, dp, seq, minLoopLength):
    return (max(dp[i + 1][j],
                dp[i][j - 1],
                dp[i + 1][j - 1] + delta(i, j, seq, minLoopLength),
                fourth_case(i, j, dp)))


def printMatrix(dp):
    for i in range(len(dp)):
        print(dp[i])


def traceback(i, j, dp, seq, struc, minLoopLength):
    if i < j:
        if dp[i][j] == dp[i+1][j]:
            traceback(i+1, j, dp, seq, struc, minLoopLength)
        elif dp[i][j] == dp[i][j-1]:
            traceback(i, j-1, dp, seq, struc, minLoopLength)
        elif dp[i][j] == dp[i+1][j-1] + delta(i, j, seq, minLoopLength):
            struc.append((i+1, j+1))
            traceback(i+1, j-1, dp, seq, struc, minLoopLength)
        else:
            for k in range(i+1, j-1):
                if dp[i][j] == dp[i][k] + dp[k+1][j]:
                    traceback(i, k, dp, seq, struc, minLoopLength)
                    traceback(k+1, j, dp, seq, struc, minLoopLength)
                    break


def createBpseq(seq, bp, maxScore):
    out = ""
    out += "Filename: " + str(args.input) + "\n"
    out += "Min-Loop: " + str(args.min_loop_length) + "\n"
    out += "GC: " + str(args.score_GC) + "\n"
    out += "AU: " + str(args.score_AU) + "\n"
    out += "GU: " + str(args.score_GU) + "\n"
    out += "Score: " + str(maxScore) + "\n"
    out += "\n"

    for i in range(0, len(seq)):
        out += str(i+1) + " " + seq[i]
        tpl = [(x, y) for x, y in bp if x == i + 1]
        if tpl:
            out += " " + str(tpl[0][1]) + "\n"
        else:
            out += " 0" + "\n"


    return out


def createDotBracket(seq, bp, maxScore):
    out = "\n"
    out += seq + "\n"

    for i in range(0, len(seq)):
        tpl = [(x, y) for x, y in bp if x == i + 1]
        tpl2 = [(x, y) for x, y in bp if y == i + 1]
        if tpl:
            out += "("
        elif tpl2:
            out += ")"
        else:
            out += "."

    return out


def nussinov(seq, minLoopLength):
    n = len(seq)
    dp = initializeMatrix(n)
    bp = []

    # Fill the matrix
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            dp[i][j] = gamma(i, j, dp, seq, minLoopLength)

    # Traceback
    traceback(i, j, dp, seq, bp, minLoopLength)

    maxScore = dp[0][n - 1]

    #printMatrix(dp)

    print(createBpseq(seq, bp, maxScore))

    print(createDotBracket(seq, bp, maxScore))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Predicting rna secondary structure with Nussinov algorithm""")
    parser.add_argument('-i', "--input", help='Input file in fasta format')
    parser.add_argument("--min-loop-length", default=3, help='Set minimum loop length')
    parser.add_argument("--score-GC", default=1, help='Set GC-score for scoring function')
    parser.add_argument("--score-AU", default=1, help='Set AU-score for scoring function')
    parser.add_argument("--score-GU", default=1, help='Set Gu-score for scoring function')

    args = parser.parse_args()

    seq = ""
    # Read Fasta file
    with open("/Users/Tobias/Documents/Studium/SS_18/SSBI/Assignments/test.fasta") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                continue
            else:
                seq += line.strip()

    nussinov(seq, int(args.min_loop_length))
