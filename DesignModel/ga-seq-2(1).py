import random
import numpy as np

# 随机生成num个长30bp的序列
def ori_popular(num):
    popular = [[] for _ in range(num)]
    for i in range(num):
        for j in range(0, 30):
            x = random.randint(0, 3)
            # 0123分别代表ATCG
            popular[i].append(x)
    return popular


def select2_2(popular_gene, select_score):
    # 30个clusters的pssm
    #  pssm_clusters=[[20.13, 19.5, 17.79, 18.43], [20.27, 19.58, 17.71, 18.09], [20.24, 19.56, 17.7, 18.18], [20.29, 19.5, 17.72, 18.17], [20.34, 19.49, 17.94, 17.86], [20.31, 19.64, 17.71, 17.89], [20.27, 19.59, 17.79, 17.99], [20.27, 19.63, 17.87, 17.85], [20.16, 19.75, 17.71, 18.02], [20.46, 19.42, 17.4, 18.17], [20.49, 19.35, 17.37, 18.24], [20.47, 19.36, 17.24, 18.36], [20.62, 18.82, 17.14, 18.79], [20.86, 18.68, 16.57, 18.74], [20.93, 18.29, 15.96, 19.17], [20.86, 17.4, 15.56, 19.89], [20.45, 16.67, 15.09, 20.68], [19.76, 15.95, 14.24, 21.41], [19.46, 15.86, 13.77, 21.61], [19.57, 16.59, 14.42, 21.39], [19.04, 17.85, 14.85, 21.3], [19.2, 18.8, 15.93, 20.69], [19.84, 19.02, 16.2, 20.01], [20.32, 19.11, 16.52, 19.28], [20.27, 19.59, 17.3, 18.39], [19.92, 19.79, 17.92, 18.23], [20.36, 18.77, 17.97, 18.76], [20.55, 19.18, 17.7, 18.08], [19.76, 19.42, 18.99, 17.94], [22.25, 17.35, 8.81, 16.8]]
    # for i in range(0,10):
    #   pssm_i = read_cluster_to_pssm(i)
    #   pssm_clusters.append(pssm_i)
    # 一共有30个cluster，这里改成10试试
    pssm_clusters = np.loadtxt('pssm_cluster_2023.txt')
    pssm_clusters = pssm_clusters.transpose()
    new_generation = popular_gene
    for i in range(10):
        new_generation = select2(new_generation, select_score, pssm_clusters[i])
    return new_generation

# 与select相比后面改成了小于
def select2(popular_gene, select_score, pssm): # 基因与目标个体的距离，要与每个聚类平均距离最大
    fitness = []
    for i in range(len(popular_gene)):
        value = 0
        for j in range(0, len(popular_gene[0])):
            if popular_gene[i][j] == 0:
                value = value + pssm[j]
            if popular_gene[i][j] == 1:
                value = value + pssm[j]
            if popular_gene[i][j] == 2:
                value = value + pssm[j]
            if popular_gene[i][j] == 3:
                value = value + pssm[j]
        fitness.append(value)
    new_popular = []
    for i in range(0, len(fitness)):
        # print(fitness[i])
        if fitness[i] < select_score:      # 这里改成了小于
            new_popular.append(popular_gene[i])
    # print(fitness)
    return new_popular

# 适应度函数。
def select(popular_gene, select_score):
    fitness = []
    # pssm的顺序为ATCG
    # pssm[i][j] = round(2 * math.log(pssm[i][j] / 4501 * 4, 2), 2)
    # pssm = [[20.13, 19.5, 17.79, 18.43], [20.27, 19.58, 17.71, 18.09], [20.24, 19.56, 17.7, 18.18], [20.29, 19.5, 17.72, 18.17], [20.34, 19.49, 17.94, 17.86], [20.31, 19.64, 17.71, 17.89], [20.27, 19.59, 17.79, 17.99], [20.27, 19.63, 17.87, 17.85], [20.16, 19.75, 17.71, 18.02], [20.46, 19.42, 17.4, 18.17], [20.49, 19.35, 17.37, 18.24], [20.47, 19.36, 17.24, 18.36], [20.62, 18.82, 17.14, 18.79], [20.86, 18.68, 16.57, 18.74], [20.93, 18.29, 15.96, 19.17], [20.86, 17.4, 15.56, 19.89], [20.45, 16.67, 15.09, 20.68], [19.76, 15.95, 14.24, 21.41], [19.46, 15.86, 13.77, 21.61], [19.57, 16.59, 14.42, 21.39], [19.04, 17.85, 14.85, 21.3], [19.2, 18.8, 15.93, 20.69], [19.84, 19.02, 16.2, 20.01], [20.32, 19.11, 16.52, 19.28], [20.27, 19.59, 17.3, 18.39], [19.92, 19.79, 17.92, 18.23], [20.36, 18.77, 17.97, 18.76], [20.55, 19.18, 17.7, 18.08], [19.76, 19.42, 18.99, 17.94], [22.25, 17.35, 8.81, 16.8]]
    pssm = np.loadtxt('pssm2023.txt', delimiter='\t')
    pssm = np.transpose(pssm)
    for i in range(len(popular_gene)):
        value = 0
        for j in range(0, len(popular_gene[0])):
            if popular_gene[i][j] == 0:
                value = value + pssm[0][j]
            if popular_gene[i][j] == 1:
                value = value + pssm[1][j]
            if popular_gene[i][j] == 2:
                value = value + pssm[2][j]
            if popular_gene[i][j] == 3:
                value = value + pssm[3][j]
        fitness.append(value)
    new_popular = []
    for i in range(0, len(fitness)):
        if fitness[i] > select_score:
            new_popular.append(popular_gene[i])
    print(fitness)
    # 使用zip()配合sorted()结合索引和值
    items = list(zip(range(len(fitness)), fitness))
    items = sorted(items, key=lambda x: x[1], reverse=True)

    # 提取前3个最大值的索引
    top_3_indices = [item[0] for item in items[:3]]
    best_gene = [popular_gene[i] for i in top_3_indices]
    return new_popular, best_gene


def crossover(population, pc):
    # 两两进行交换产生子序列
    pop_len = len(population)
    for i in range(0, pop_len - 1, 1):
        cpc = random.random()
        # print(cpc)
        if cpc < pc:
            cpoint = random.randint(0, len(population[0]))
            # print(cpoint)
            temporary1 = []
            temporary2 = []
            temporary1.extend(population[i][0:cpoint])
            temporary1.extend(population[i + 1][cpoint:len(population[i])])
            temporary2.extend(population[i + 1][0:cpoint])
            temporary2.extend(population[i][cpoint:len(population[i])])
            population.append(temporary1)
            population.append(temporary2)


def mutation(population, pm):
    pop_len = len(population)
    for i in range(0, pop_len, 1):
        for j in range(0, len(population[0]), 1):
            cpm = random.random()
            if cpm < pm:
                x = random.randint(1, 3)
                population[i][j] = (population[i][j] + x) % 4
                # print(i, j)


def geneticAlgorithm(num, select_score, select_score2,target_score, pc, pm, generation):
    ori_popular_gene = ori_popular(num)
    new_popular_gene = ori_popular_gene
    rate = pow(target_score - select_score, 1 / generation)
    for i in range(0, generation):
        select_score = select_score * rate
        # 严格筛选条件
        new_popular_gene, last_best_gene = select(new_popular_gene, select_score)
        # 这里加了select2_2
        # new_popular_gene = select2_2(new_popular_gene, select_score2)
        crossover(new_popular_gene, pc)
        mutation(new_popular_gene, pm)
        # 保留上一代最优的3个序列
        new_popular_gene.extend(last_best_gene)
        new_popular_gene = np.unique(new_popular_gene, axis=0)
        print('Process generation', i)
        print('当前种群数量', len(new_popular_gene))
    # print(new_popular_gene)
    res = []
    for i in range(0, len(new_popular_gene)):
        gene_str = ''
        for j in range(0, len(new_popular_gene[0])):
            if new_popular_gene[i][j] == 0:
                gene_str = gene_str + 'A'
            if new_popular_gene[i][j] == 1:
                gene_str = gene_str + 'T'
            if new_popular_gene[i][j] == 2:
                gene_str = gene_str + 'C'
            if new_popular_gene[i][j] == 3:
                gene_str = gene_str + 'G'
        res.append(gene_str)
    # print(res)
    return res


# num表示随机生成的基因数目
# select_score用于fitness函数，若某一序列的fitness.value小于select_score则被淘汰
# pc表示发生某一条序列与下一条交叉并产生两条新序列的概率，pm表示每个碱基发生变异的概率
# 随机生成序列所需时间较长
# target_score为目标序列的
res = geneticAlgorithm(num=9000,
                       select_score=1,
                       select_score2=100,
                       target_score=25,
                       pc=0.1,
                       pm=0.01,
                       generation=50)
print('共得到序列个数：', len(res))
fo = open("ga_RBS_seq.txt", "w")
for i in range(0, len(res)):
    fo.write(res[i])
    fo.write("\n")
fo.close()