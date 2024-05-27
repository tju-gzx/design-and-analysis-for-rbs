

# # 从 fasta 序列文件中读取所有序列
# with open('rbs_seqs.fasta') as f:
#     content = f.read()
#
# # 提取所有序列的核苷酸序列
# seqs = re.findall(r'\n([A-Z]+)', content)
#
# # 计算频率矩阵
# freq = [[0] * 4 for i in range(30)]
# for seq in seqs:
#     for i in range(len(seq)):
#         if seq[i] == 'A':
#             freq[i][0] += 1
#         elif seq[i] == 'T':
#             freq[i][1] += 1
#         elif seq[i] == 'C':
#             freq[i][2] += 1
#         elif seq[i] == 'G':
#             freq[i][3] += 1
#
# # 计算背景概率值
# total_count = sum(item for row in freq for item in row)
# bg_freq = [sum(freq[i]) / total_count for i in range(len(freq))]
#
# # 将频率矩阵转化为百分比矩阵
# total_seqs = len(seqs)
# for i in range(len(freq)):
#     total_count = sum(freq[i])
#     for j in range(len(freq[i])):
#         freq[i][j] /= total_count
#         freq[i][j] *= 100
#
# # 计算 PSSM 矩阵
# pssm = [[0] * 4 for i in range(30)]
# for i in range(len(freq)):
#     for j in range(len(freq[i])):
#         pssm[i][j] = round(2 * math.log2(freq[i][j] / bg_freq[i]), 2)
#
# # 将 PSSM 矩阵写入文件
# with open('pssm.txt', 'w') as f:
#     for row in pssm:
#         f.write('\t'.join(str(val) for val in row))
#         f.write('\n')
#
# print(pssm)
import re
import pandas as pd
import math

file = open('RBSclu.fasta', 'r')
res = open('pssm_cluster_2023.txt', 'w')
content = file.read()
seq = re.findall(r'\n[A-Z]+', content)
pssm = [[0.0] * 4 for i in range(30)]
for i in range(len(seq)):
    for j in range(1, len(seq[i])):
        if seq[i][j] == 'A':
            pssm[j - 1][0] += 1
        elif seq[i][j] == 'T':
            pssm[j - 1][1] += 1
        elif seq[i][j] == 'C':
            pssm[j - 1][2] += 1
        elif seq[i][j] == 'G':
            pssm[j - 1][3] += 1

# Apply pseudocounts
for i in range(30):
    for j in range(4):
        pssm[i][j] += 1

# Normalize counts and compute log-odds scores
for i in range(30):
    total_count = sum(pssm[i])
    for j in range(4):
        pssm[i][j] = round(2 * math.log2(pssm[i][j] / len(seq) * 4), 2)

# Write PSSM to file
for i in range(30):
    res.write('\t'.join(map(str, pssm[i])) + '\n')

# Close files
file.close()
res.close()






