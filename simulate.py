import anndata as ad
import numpy as np
import pandas as pd
import matplotlib1.pyplo as plt

# a)基于1）中单细胞转录组数据的单细胞为模版进行抽样模拟；
# b)以华大Stereo-Seq时空转录组测序芯片500nm制式为模版；
# c)规格设置为1mm * 1mm尺寸芯片，不要求考虑模拟trackline，但要了解trackline排布规律及用途；
# d)模拟的细胞密度范围为50-100个细胞/平方毫米，要求随机确定一个数值；
# e)随机化确定细胞中心坐标，邻域直径参考人类细胞，要求在10-20微米（不要求特定细胞形状，可以为圆形或正方形）；
# f)在随机中心的直径领域内产生指定数目的随机坐标（X，Y）。
# g)模拟的细胞类型和数目要求从单细胞中top5细胞类型选取并且根据各类型比例加权获得数目，要求从各单细胞种类中随机抽取相应数目的单细胞ID；
# h)模拟的每个细胞基因数据需随机抽取从单细胞RNA中抽取500-1000个基因表达的基因标签，要求筛选条件：Count>1, Top 2000；
# i)给每个抽取的单细胞分配坐标中心，给每个细胞中的基因标签分配位置坐标；
# j)结果保存成华大空间转录组数据GEM格式（.gem）文件，（参考提示1）；
# k)并绘制芯片转录组图像，（参考提示2）。

# 操作步骤：设定参数——计算细胞总数量并选取——确定细胞坐标——计算芯片上spot数量并分配坐标——计算每个细胞中的基因数并选取——将spot随机分配给基因——随机分配基因表达量

# 读取h5ad文件
adata = ad.read_h5ad('result/integrated_data.h5ad')

# 设定参数
chip_size_mm = 1  # 1mm = 1,000,000nm
cell_density_range_mm = (50, 100)
cell_diameter_range_um = (10, 20)  # diameter 10-20um = 10,000-20,000nm
spot_spacing_nm = 500
gene_per_cell_density_range = (500, 1000)


# 第一个spot距离芯片边界的距离：（500 - 220）/ 2 = 140nm
# 每行每列包含的spot数量： 1,000,000 / 500 = 2,000个； 芯片上spot数量： 2,000 * 2,000 = 4,000,000个


# 模拟基因表达
def sample_gene(cell_id, adata, num_genes, top_genes=2000):
    # 获取细胞中表达量大于1的基因的索引
    cell_index = np.where(adata.obs_names == cell_id)[0][0]
    cell_expression = adata.X[cell_index, :]
    genes_above_1 = np.where(cell_expression > 1)[0]  # 筛选表达量大于1的基因的索引

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


# 获取细胞包含的spot点及其坐标
def cell_contains_spots(cell_center_x, cell_center_y, cell_size_nm, spot_spacing_nm, chip_size_mm, spot_size_nm):
    # 转换单位
    chip_size_nm = chip_size_mm * 1e6

    # 计算细胞边界
    half_cell_size = cell_size_nm // 2
    cell_min_x = cell_center_x - half_cell_size
    cell_max_x = cell_center_x + half_cell_size
    cell_min_y = cell_center_y - half_cell_size
    cell_max_y = cell_center_y + half_cell_size

    # 初始化spot坐标列表
    spots_in_cell = []

    # 计算spot的起始索引（基于间隔和边界条件）
    start_x = max(0, int((cell_min_x - 140) // spot_spacing_nm))
    end_x = min(int((cell_max_x + 140) // spot_spacing_nm), chip_size_nm // spot_spacing_nm)
    start_y = max(0, int((cell_min_y - 140) // spot_spacing_nm))
    end_y = min(int((cell_max_y + 140) // spot_spacing_nm), chip_size_nm // spot_spacing_nm)

    # 遍历spot
    for i in range(start_x, end_x):
        for j in range(start_y, end_y):
            spot_x = i * spot_spacing_nm + spot_size_nm // 2 + 140  # spot中心坐标
            spot_y = j * spot_spacing_nm + spot_size_nm // 2 + 140

            # 检查spot是否在细胞内
            if cell_min_x <= spot_x <= cell_max_x and cell_min_y <= spot_y <= cell_max_y:
                spots_in_cell.append((spot_x, spot_y))

    return spots_in_cell


# 检查新生成的细胞是否与存在的细胞重叠    不重叠返回True，重叠返回False
def is_overlap(new_cell, cells):
    for cell in cells:
        if abs(new_cell[0] - cell[0]) < cell_diameter_range_um[1] * 1e3 and abs(new_cell[1] - cell[1]) < cell_diameter_range_um[1] * 1e3:
            return True
        return False


def simulate(adata, chip_size_mm, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm,
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

    # 统计出现次数top5的细胞类型并计算其权重
    top_cell_types = adata.obs['cell_type'].value_counts().head(5)
    top_cell_types_index = top_cell_types.index
    type_proportion = top_cell_types / top_cell_types.sum()

    # 计算细胞总数量
    num_cell = chip_size_mm ** 2 * np.random.randint(low=cell_density_range_mm[0], high=cell_density_range_mm[1] + 1)

    # 生成细胞类型和位置
    cell_position = []
    margin = cell_diameter_range_um[1] / 2 * 1e3
    while len(cell_position) < num_cell:
        cell_position_x = np.random.rand() * (chip_size_mm * 1e6 - cell_diameter_range_um[1] * 1e3) + margin
        cell_position_y = np.random.rand() * (chip_size_mm * 1e6 - cell_diameter_range_um[1] * 1e3) + margin
        new_cell = (cell_position_x, cell_position_y)
        if not is_overlap(new_cell, cell_position):
            cell_position.append(new_cell)
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

    # bug：cell_types里的细胞类型和cell_ids里的细胞标签没有对应（改成用细胞标签进行循环）
    for i, (cell_id, cell_type) in enumerate(cell_ids_types):
        # 获取当前细胞的位置
        pos = cell_position[i]

        # 确定细胞大小并从每个细胞中抽取基因
        cell_size = np.random.randint(low=cell_diameter_range_um[0], high=cell_diameter_range_um[1] + 1) * 1e3
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
    gem_df.to_csv('simulated_spatial_transcriptomics.gem', sep='\t', index=False, header=True)  # 使用制表符分隔，无索引，无表头

    # 绘制细胞位置
    plt.figure(figsize=(10, 10))
    # 创建一个颜色映射
    colormap = plt.cm.viridis
    norm = plt.Normalize(vmin=0, vmax=len(spatial_df['cell_type'].unique()) - 1)
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])  # 初始化标量映射对象

    # 创建一个颜色字典
    cell_type_to_color = {cell_type: colormap(norm(i)) for i, cell_type in
                          enumerate(sorted(spatial_df['cell_type'].unique()))}

    # 绘制散点图
    scatters = []
    for cell_type in sorted(spatial_df['cell_type'].unique()):
        subset = spatial_df[spatial_df['cell_type'] == cell_type]
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


simulate(adata, chip_size_mm, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm,
         gene_per_cell_density_range)
