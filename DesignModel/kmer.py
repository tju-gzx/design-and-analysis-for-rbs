from collections import defaultdict
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
def count_kmers(seq, k):
    """计算序列中所有k-mer的出现次数"""
    kmer_counts = defaultdict(int)
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

def kmer_frequencies(seq_file, ks):
    """计算指定序列文件中指定k值的k-mer频率"""
    kmer_freqs = {}
    for k in ks:
        kmer_freqs[k] = {}
    with open(seq_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq).upper()
            total_kmers = len(seq) - k + 1
            for k in ks:
                kmer_counts = count_kmers(seq, k)
                for kmer, count in kmer_counts.items():
                    kmer_freqs[k][kmer] = kmer_freqs[k].get(kmer, 0) + count / total_kmers
    return kmer_freqs

def js_divergence(p, q):
    """计算两个概率分布的JS散度"""
    m_keys = set(p.keys()).union(set(q.keys()))
    m = np.array([0.5 * (p.get(k, 0) + q.get(k, 0)) for k in m_keys])
    kl_pm = np.sum([p.get(k, 0) * np.log2(p.get(k, 0) / m[i]) for i, k in enumerate(m_keys)])
    kl_qm = np.sum([q.get(k, 0) * np.log2(q.get(k, 0) / m[i]) for i, k in enumerate(m_keys)])
    return 0.5 * (kl_pm + kl_qm)

# 读入两个序列文件
seq_file1 = 'GAselect2.fasta'
seq_file2 = 'rbs_seqs.fasta'

# 计算指定K值下的K-mer频率
ks = [2,4,6,8]
kmer_freqs1 = kmer_frequencies(seq_file1, ks)
kmer_freqs2 = kmer_frequencies(seq_file2, ks)

# # 计算每个K值下的JS散度，并保存序列1和序列2的K-mer频率
# js_divergences = []
# freqs1 = []
# freqs2 = []
# for k in ks:
#     p = kmer_freqs1[k]
#     q = kmer_freqs2[k]
#     jsd = js_divergence(p, q)
#     js_divergences.append(jsd)
#     for kmer in p.keys():
#         freqs1.append(p[kmer])
#         freqs2.append(q.get(kmer, 0))
#     print(f'K={k}, JS Divergence={jsd}')

# # 绘制散点图
# import matplotlib.pyplot as plt
# plt.scatter(freqs1, freqs2)
# plt.xlabel('GA k-mer frequency')
# plt.ylabel('Origin k-mer frequency')
# plt.title('JS Divergence Scatter Plot')
# plt.show()
# # 绘制子图
# fig, axs = plt.subplots(2, 2, figsize=(10, 10))
# for i, k in enumerate(ks):
#     row = i // 2
#     col = i % 2
#     p = kmer_freqs1[k]
#     q = kmer_freqs2[k]
#     jsd = js_divergence(p, q)
#     print(f'K={k}, JS Divergence={jsd}')
#     freqs1 = []
#     freqs2 = []
#     for kmer in p.keys():
#         freqs1.append(p[kmer])
#         freqs2.append(q.get(kmer, 0))
#     axs[row, col].scatter(freqs1, freqs2)
#     axs[row, col].set_xlabel('GA k-mer frequency')
#     axs[row, col].set_ylabel('Origin k-mer frequency')
#     axs[row, col].set_title(f'JS Divergence Scatter Plot (K={k})')
# plt.show()
# 计算每个K值下的JS散度，并保存序列1和序列2的K-mer频率
js_divergences = []
freqs1 = []
freqs2 = []
epsilon = 0.0001  # 设置一个极小值
for k in ks:
    p = kmer_freqs1[k]
    q = kmer_freqs2[k]
    # 计算两个序列中所有K-mer的并集
    kmers = set(p.keys()) | set(q.keys())
    # 归一化K-mer频率
    p_sum = sum(p.get(kmer, epsilon) for kmer in kmers)
    q_sum = sum(q.get(kmer, epsilon) for kmer in kmers)
    p = {kmer: (p.get(kmer, 0) + epsilon) / p_sum for kmer in kmers}
    q = {kmer: (q.get(kmer, 0) + epsilon) / q_sum for kmer in kmers}
    jsd = js_divergence(p, q)
    js_divergences.append(jsd)
    for kmer in kmers:
        freqs1.append(p.get(kmer, 0))
        freqs2.append(q.get(kmer, 0))
    print(f'K={k}, JS Divergence={jsd}')