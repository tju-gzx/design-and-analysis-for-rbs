from Bio import SeqIO
import pandas as pd

# # 读取Fasta文件并将序列保存到列表中
# sequences = []
# with open("born_rbs.fasta", "r") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         sequences.append(str(record.seq))
#
# # 将列表转换为DataFrame并保存为CSV文件
# df = pd.DataFrame(sequences, columns=["Sequence"])
# df.to_csv("born_rbs.csv", index=False)

import pandas as pd

# 读取文本文件并将其保存为DataFrame
df = pd.read_csv("2000GAborn.txt", sep="\t", header=None)

# 将DataFrame-保存为CSV文件
df.to_csv("2000GA.csv", index=False)