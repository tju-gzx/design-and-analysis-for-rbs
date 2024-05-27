# 打开文件
with open('2000trainpssm.txt', 'r') as f:
        # 读取文件中的所有数字并转换成float类型
        numbers = []
        for line in f.readlines():
            # 将点号替换成空格
            line = line.replace('/t', ' ')
            # 将每行字符串转换成浮点数
            nums = [float(x) for x in line.strip().split()]
            numbers.extend(nums)

# 计算所有数字之和
total = sum(numbers)

# 输出结果
print('所有数字之和为:', total)
