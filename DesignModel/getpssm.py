import re
import pandas as pd
import math

file = open('2000train.fasta', 'r')
res = open('2000pssm.txt', 'w')
content = file.read()
seq = re.findall(r'\n[A-Z]+', content)
pssm = [[0] * 61, [0] * 61, [0] * 61, [0] * 61]
for i in range(0, 4501):
    # print(seq[i])
    for j in range(1, 62):
        if seq[i][j] == 'A':
            pssm[0][j - 1] = pssm[0][j - 1] + 1
        if seq[i][j] == 'T':
            pssm[1][j - 1] = pssm[1][j - 1] + 1
        if seq[i][j] == 'C':
            pssm[2][j - 1] = pssm[2][j - 1] + 1
        if seq[i][j] == 'G':
            pssm[3][j - 1] = pssm[3][j - 1] + 1
print(pssm)
for i in range(0, 4):
    for j in range(0, 61):
        pssm[i][j] = round(2 * math.log(pssm[i][j] / 4501 * 4, 2), 2)
for i in range(0, len(pssm)):
    res.write(str(pssm[i]))
    res.write("\n")
