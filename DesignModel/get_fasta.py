
import pandas as pd
import pandas as pd
#
# # 读取 Excel 文件中的数据
# df = pd.read_excel('2000训练集.xlsx')
#
# # 将数据写入文本文件
# with open('2000train.txt', 'w') as f:
#     for index, row in df.iterrows():
#         line = '\t'.join([str(x) for x in row]) + '\n'
#         f.write(line)

# 读取txt文件，假设每行为一个序列
with open("rbs-select2.txt", "r") as f:
    sequences = f.readlines()

# 将每个序列转换为fasta格式
fasta_sequences = []
for i, seq in enumerate(sequences):
    title = f">Sequence{i+1}"
    fasta_seq = f"{title}\n"
    while len(seq) > 80:
        fasta_seq += f"{seq[:80]}\n"
        seq = seq[80:]
    fasta_seq += seq.strip()
    fasta_sequences.append(fasta_seq)

# 将fasta格式的序列写入文件
with open("GAselect2.fasta", "w") as f:
    f.write("\n".join(fasta_sequences))



