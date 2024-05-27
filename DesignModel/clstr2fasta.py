import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def write_sequences_to_fasta(first_sequences, output_file):
    with open(output_file, 'w') as f_out:
        for seq_id, seq in first_sequences.items():
            seq_record = SeqRecord(Seq(seq), id=seq_id)
            SeqIO.write(seq_record, f_out, 'fasta')

# 读取clstr文件并提取每个Cluster组中第一行的序列名
def extract_first_sequence_names(clstr_file):
    first_sequence_names = []
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                cluster_id = line.strip().split()[1]
            else:
                seq_info = line.strip().split('\t')[-1]
                if seq_info.endswith('*'):
                    seq_id = re.search(r'>\S+', seq_info).group()[1:].split()[0].rstrip('...')
                    if seq_id not in first_sequence_names:
                        first_sequence_names.append(seq_id)
    return first_sequence_names

# 从fasta文件中提取每个Cluster组中第一行的序列的DNA序列
def extract_first_sequences(fasta_file, first_sequence_names):
    first_sequences = {}
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        if seq_record.id in first_sequence_names:
            first_sequences[seq_record.id] = str(seq_record.seq)
    return first_sequences

# 测试代码
fasta_file = 'rbs_seqs.fasta'
clstr_file = 'rbs.clstr'

# 提取每个Cluster组中第一行的序列名
first_sequence_names = extract_first_sequence_names(clstr_file)

# 从fasta文件中提取每个Cluster组中第一行的序列的DNA序列
first_sequences = extract_first_sequences(fasta_file, first_sequence_names)


# 将提取的序列写入新的fasta文件中
write_sequences_to_fasta(first_sequences, 'RBSclu.fasta')

# 打印每个Cluster组中第一行的序列的DNA序列
for seq_id, seq in first_sequences.items():
    print(f'Sequence ID: {seq_id}')
    print(f'Sequence: {seq}\n')