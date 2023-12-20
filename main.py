import numpy as np

# Redefine the scoring matrix for the alignments, including the gap scoring
scoring_matrix = {
    'A': {'A': 1, 'G': -0.8, 'T': -0.2, 'C': -2.3, '-': -0.6},
    'G': {'A': -0.8, 'G': 1, 'T': -1.1, 'C': -0.7, '-': -1.5},
    'T': {'A': -0.2, 'G': -1.1, 'T': 1, 'C': -0.5, '-': -0.9},
    'C': {'A': -2.3, 'G': -0.7, 'T': -0.5, 'C': 1, '-': -1},
    '-': {'A': -0.6, 'G': -1.5, 'T': -0.9, 'C': -1, '-': 0}  # Set gap-gap score to 0
}

# Redefine the sequences to be aligned
seq1 = "TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC"
seq2 = "AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"

# Initialize the matrix for dynamic programming
dp_matrix = np.zeros((len(seq1)+1, len(seq2)+1))

# Initialize the first row and column of the dp_matrix according to gap penalties
for i in range(1, len(seq1)+1):
    dp_matrix[i][0] = dp_matrix[i-1][0] + scoring_matrix[seq1[i-1]]['-']
for j in range(1, len(seq2)+1):
    dp_matrix[0][j] = dp_matrix[0][j-1] + scoring_matrix['-'][seq2[j-1]]

# Fill the dynamic programming matrix and compute the score
for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
        match = dp_matrix[i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
        delete = dp_matrix[i-1][j] + scoring_matrix[seq1[i-1]]['-']
        insert = dp_matrix[i][j-1] + scoring_matrix['-'][seq2[j-1]]
        dp_matrix[i][j] = max(match, delete, insert)

# Function to trace back the optimal alignment
def traceback(seq1, seq2, dp_matrix):
    align1, align2 = '', ''
    i, j = len(seq1), len(seq2)
    while i > 0 and j > 0:
        score_current = dp_matrix[i][j]
        score_diagonal = dp_matrix[i-1][j-1]
        score_up = dp_matrix[i][j-1]
        score_left = dp_matrix[i-1][j]
        if score_current == score_diagonal + scoring_matrix[seq1[i-1]][seq2[j-1]]:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + scoring_matrix[seq1[i-1]]['-']:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        else:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1
    return align1[::-1], align2[::-1]

# Get the optimal alignment
alignment1, alignment2 = traceback(seq1, seq2, dp_matrix)

# The highest score is at the bottom-right corner of the matrix
highest_score = dp_matrix[-1][-1]

# Print the results
highest_score, alignment1, alignment2

print("Highest Score:", highest_score)
print(alignment1)
print(alignment2)