def read_matrix_dimensions(file_path):
    with open(file_path, 'r') as file:
        # 读取所有行到一个列表中
        lines = file.readlines()

        # 获取列数，即第一行的元素数量
        num_columns = len(lines[0].split())

        # 获取行数，即文件中的行数
        num_rows = len(lines)

        return num_rows, num_columns

# 使用函数并打印结果
file_path = '../data/calibration/matrix_output_m_128_k_768.txt'  # 替换为你的文件路径
rows, columns = read_matrix_dimensions(file_path)
print(f"Matrix dimensions: {rows} x {columns}")