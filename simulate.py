import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

# 读取h5ad文件
# adata = ad.read_h5ad('integrated_data.h5ad')

# # 设定参数
# chip_size_mm = 1  # 1mm = 1,000,000nm
# cell_density_range_mm = (50, 100)
# cell_diameter_range_um = (10, 20)  # diameter 10-20um = 10,000-20,000nm
# spot_spacing_nm = 500
# gene_per_cell_density_range = (500, 1000)


# 模拟基因表达
def sample_gene(cell_id, adata, num_genes, gene_expression=1, top_genes=2000):
    """
    :param cell_id: 细胞ID
    :param adata: Anndata数据结构的数据
    :param num_genes: 每个细胞需要抽取的基因数量
    :param top_genes: 筛选条件，从表达量前2000的基因中抽取
    :param gene_expression: 筛选条件，从表达量大于 gene_expression 的基因中抽取
    :return: 返回所有抽取的基因名称及其表达量
    """
    # 获取细胞中表达量大于1的基因的索引
    cell_index = np.where(adata.obs_names == cell_id)[0][0]
    cell_expression = adata.X[cell_index, :]
    genes_above_1 = np.where(cell_expression > gene_expression)[0]  # 筛选表达量大于1的基因的索引

    # 如果基因少于top_genes，则直接使用所有表达量大于1的基因
    if len(genes_above_1) < top_genes:
        top_gene_indices = genes_above_1
    else:
        # 对基因的表达量进行排序，并取前top_genes个
        sorted_indices = np.argsort(cell_expression[genes_above_1])[::-1]  # 输出从大到小的索引
        top_gene_indices = genes_above_1[sorted_indices[:top_genes]]

    # 从高表达基因中随机选取num_genes个基因
    sampled_genes_indices = np.random.choice(top_gene_indices, size=num_genes, replace=False)

    # 返回抽取到的基因的名称和对应的基因表达量
    sampled_genes = adata.var_names[sampled_genes_indices]
    sampled_expression = cell_expression[sampled_genes_indices]

    return sampled_genes, sampled_expression


def is_overlap(new_cell, cells):
    """
    检查新生成的细胞是否与已存在的细胞重叠。
    """
    for cell in cells:
        if abs(new_cell[0] - cell[0]) < (new_cell[2] / 2 + cell[2] / 2) and \
                abs(new_cell[1] - cell[1]) < (new_cell[2] / 2 + cell[2] / 2):
            return True
    return False


def generate_cell(cell_diameter_range_um, chip_size_mm, cell_position):
    # 确定细胞的大小
    while True:
        cell_size = np.random.randint(low=cell_diameter_range_um[0], high=cell_diameter_range_um[1] + 1) * 1e3
        margin = cell_size / 2

        # 生成细胞位置坐标
        cell_position_x = np.random.rand() * (chip_size_mm * 1e6 - cell_size) + margin
        cell_position_y = np.random.rand() * (chip_size_mm * 1e6 - cell_size) + margin
        new_cell = (cell_position_x, cell_position_y, cell_size)

        # 检查是否重叠
        if not is_overlap(new_cell, cell_position):
            cell_position.append(new_cell)
            return new_cell, cell_size


def draw_picture(spatial):
    """
    :param spatial: DataFrame数据
    :return: 将细胞位置通过散点图展示，并保存图像
    """
    # 绘制细胞位置
    plt.figure(figsize=(10, 10))
    # 创建颜色映射
    colormap = plt.cm.viridis
    norm = plt.Normalize(vmin=0, vmax=len(spatial['cell_type'].unique()) - 1)
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])

    # 创建颜色字典
    cell_type_to_color = {cell_type: colormap(norm(i)) for i, cell_type in
                          enumerate(sorted(spatial['cell_type'].unique()))}

    # 绘制散点图
    scatters = []
    for cell_type in sorted(spatial['cell_type'].unique()):
        subset = spatial[spatial['cell_type'] == cell_type]
        sizes = (subset['cell_size'] ** 2) / 1e6
        scatters.append(plt.scatter(subset['x_position'] / 1e6, subset['y_position'] / 1e6,
                                    s=sizes, alpha=0.5, color=cell_type_to_color[cell_type],
                                    marker='s', label=cell_type))

    # 添加图例
    plt.legend(handles=scatters, title='Cell Type')

    plt.title('Simulated Spatial Transcriptomics Chip')
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')
    plt.savefig('simulated_spatial_transcriptomics_chip.png')
    plt.show()


# 获取细胞包含的spot点及其坐标
def cell_contains_spots(cell_center_x, cell_center_y, cell_size_nm, spot_spacing_nm, chip_size_mm, spot_size_nm):
    """
    :param cell_center_x: 细胞的x坐标
    :param cell_center_y: 细胞的y坐标
    :param cell_size_nm: 细胞的大小，如范围10-20um
    :param spot_spacing_nm: 相邻两个spot点间隔距离，如500nm/715nm
    :param chip_size_mm: 芯片大小
    :param spot_size_nm: spot点大小
    :return: 获取细胞所包含的spot点及其坐标
    """
    # 转换单位
    chip_size_nm = chip_size_mm * 1e6

    # 边缘spot点与芯片边界的距离
    distance_spot_chip = (spot_spacing_nm - spot_size_nm) / 2

    # 计算细胞边界
    half_cell_size = cell_size_nm // 2
    cell_min_x = cell_center_x - half_cell_size
    cell_max_x = cell_center_x + half_cell_size
    cell_min_y = cell_center_y - half_cell_size
    cell_max_y = cell_center_y + half_cell_size

    # 初始化spot坐标列表
    spots_in_cell = []

    # 计算spot的起始索引（基于间隔和边界条件）
    start_x = max(0, int((cell_min_x - distance_spot_chip) // spot_spacing_nm))
    end_x = min(int((cell_max_x + distance_spot_chip) // spot_spacing_nm), int(chip_size_nm // spot_spacing_nm))
    start_y = max(0, int((cell_min_y - distance_spot_chip) // spot_spacing_nm))
    end_y = min(int((cell_max_y + distance_spot_chip) // spot_spacing_nm), int(chip_size_nm // spot_spacing_nm))

    # 遍历spot
    for i in range(start_x, end_x):
        for j in range(start_y, end_y):
            spot_x = i * spot_spacing_nm + spot_size_nm // 2 + distance_spot_chip  # spot中心坐标
            spot_y = j * spot_spacing_nm + spot_size_nm // 2 + distance_spot_chip

            # 检查spot是否在细胞内
            if cell_min_x <= spot_x <= cell_max_x and cell_min_y <= spot_y <= cell_max_y:
                spots_in_cell.append((spot_x, spot_y))

    return spots_in_cell


def generate_gem(spatial_df, gem_file_path):
    # 查看前几行数据
    print(spatial_df.head(5))

    # 初始化一个空的DataFrame来存储.gem格式的数据
    gem_data = []

    # 遍历 spatial_df 的每一行
    for index, row in spatial_df.iterrows():
        genes = row['genes']
        gene_expression = row['gene_expression']
        gene_spots = row['gene_spots']

        # 确保基因数量、表达量和spots的键一致
        assert len(genes) == len(gene_expression) == len(gene_spots)

        # 遍历每个基因
        for gene, expr, (x, y) in zip(genes, gene_expression, gene_spots.values()):
            gem_data.append({
                'geneID': gene,
                'x': x,
                'y': y,
                'Count': expr
            })

    # 将列表转换为DataFrame
    gem_df = pd.DataFrame(gem_data)

    # 将DataFrame保存到.gem文件
    gem_df.to_csv(gem_file_path, sep='\t', index=False, header=True)  # 使用制表符分隔，无索引，无表头


def simulate(h5ad_path, chip_size_mm, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm,
             gene_per_cell_density_range):
    """
    :param adata: Anndata数据结构文件
    :param chip_size_mm: 芯片尺寸，单位毫米
    :param cell_density_range_mm: 细胞数量取值范围，如（50-100）/平方毫米
    :param cell_diameter_range_um: 细胞的直径，单位微米,如（10-20）微米
    :param spot_spacing_nm: 相邻spot点中心的距离，单位纳米
    :param gene_per_cell_density_range: 每个细胞中抽取的基因数量范围，如（500-1000）个
    :return: 1.包含基因标签、基因坐标、基因表达量的.gem文件 2.芯片转录组图像
    """

    # 读取h5ad文件
    adata = ad.read_h5ad(h5ad_path)

    # 计算细胞总数量
    num_cell = chip_size_mm ** 2 * np.random.randint(low=cell_density_range_mm[0], high=cell_density_range_mm[1] + 1)

    # 统计出现次数top5的细胞类型并计算其权重
    top_cell_types = adata.obs['cell_type'].value_counts().head(5)
    top_cell_types_index = top_cell_types.index
    type_proportion = top_cell_types / top_cell_types.sum()

    # 生成细胞类型
    cell_types = np.random.choice(top_cell_types_index, size=num_cell, p=type_proportion)

    # 统计抽取的各细胞类型数量
    type_counts = pd.Series(cell_types).value_counts()

    # 为每种类型创建索引表，并从每种类型的细胞中随机选择
    cell_ids_types = []
    for cell_type, count in type_counts.items():
        type_mask = adata.obs['cell_type'] == cell_type
        type_ids = adata.obs.index[type_mask]
        select_ids = np.random.choice(type_ids, size=count, replace=False)
        cell_ids_types.extend([(id_, cell_type) for id_ in select_ids])

    # 存储结果
    spatial_data = []
    num_genes_per_cell = np.random.randint(low=gene_per_cell_density_range[0], high=gene_per_cell_density_range[1] + 1)

    # 初始化细胞位置
    cell_position = []
    # bug：cell_types里的细胞类型和cell_ids里的细胞标签没有对应（改成用细胞标签进行循环）
    for i, (cell_id, cell_type) in enumerate(cell_ids_types):

        # 获取当前细胞的位置
        pos, cell_size = generate_cell(cell_diameter_range_um, chip_size_mm, cell_position)
        # 确定细胞大小并从每个细胞中抽取基因
        sampled_genes, sampled_expression = sample_gene(cell_id, adata, num_genes_per_cell)

        # 获取细胞所包含的spot坐标
        spots = cell_contains_spots(pos[0], pos[1], cell_size, spot_spacing_nm, chip_size_mm, 220)

        # 使用索引来随机选择spots
        indices = np.arange(len(spots))  # 创建一个从0到len(spots)-1的索引数组
        gene_spots = {gene: spots[np.random.choice(indices, replace=False)] for gene in sampled_genes}

        # 存储数据
        spatial_data.append({
            'cell_id': cell_id,
            'x_position': pos[0],
            'y_position': pos[1],
            'cell_type': cell_type,
            'cell_size': cell_size,
            'genes': list(sampled_genes),
            'gene_expression': sampled_expression.tolist(),
            'gene_spots': gene_spots  # 存储每个基因对应的spot坐标
        })

    # 将结果转换为DataFrame
    spatial_df = pd.DataFrame(spatial_data)

    return spatial_df


# simulate(adata, chip_size_mm, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm,
#          gene_per_cell_density_range)