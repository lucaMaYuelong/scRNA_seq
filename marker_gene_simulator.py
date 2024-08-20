import numpy as np
import anndata as ad
import pandas as pd
# n种cell type
# N个cell
# 为每种cell type分配k个cell 分配方式有两种：1.随机分配(存在上限和下限)；2.指定数量分配
# m个gene
# 为每个cell分配h个gene 随机分配且符合正态分布
# 每种cell type有y个marker gene y也是参数 y的分配：1.可以是每种cell type都包含y个marker gene 2.自定义？
# marker gene: 在一种cell type里面表达量高，在剩余cell type里面表达量低；表达量高应该是表达量低的两倍及以上
# gene表达量生成：你这里可以选择生成负二项分布或者泊松分布

# 一般情况下，假设一个细胞中有1000个基因，那么其中非零表达的基因数量大约为总基因数量的10%-20%


# 模拟生成基因
def generate_gene_expression(num_gene, num_cell):
    # 非零表达基因数量为总基因数量的10%-30%
    bio_p_min = 0.09
    bio_p_max = 0.29

    # 为每个基因生成一个独立的表达概率
    p_value = np.random.uniform(bio_p_min, bio_p_max, size=num_gene)

    # 初始化基因表达矩阵，默认为0（对于未表达的基因）
    gene_expression = np.zeros((num_gene, num_cell), dtype=int)

    # 对于每个基因，生成一个长度为num_cell的二项分布数组，然后填充到gene_expression中
    for i in range(num_gene):
        successes = np.random.binomial(n=1, p=p_value[i], size=num_cell)
        gene_expression[i, :] = successes

    # 找到所有非零表达的元素的扁平化索引
    flat_indices = np.flatnonzero(gene_expression)

    # 对非零表达的基因，生成符合正态分布的表达量，并限制在1到250之间
    num_nonzero = len(flat_indices)
    random_values = np.random.normal(loc=125, scale=50, size=num_nonzero)
    random_values = np.clip(random_values, 1, 250)

    # 使用扁平化索引来更新gene_expression中的值
    gene_expression.ravel()[flat_indices] = random_values

    return gene_expression


# 模拟生成细胞及其类型 generate_way: 1.random(随机生成） 2.custom(自定义)
def generate_cell(num_cell, num_cell_type):
    cell_types = []
    for i in range(num_cell_type):
        cell_types.append(f"cell_type_{i + 1}")  # 使用 append 方法添加元素}"

    cells = []
    for i in range(num_cell):
        # 随机选择一个细胞类型
        cell_type = np.random.choice(cell_types)
        # 生成细胞ID
        cell_id = f"cell_{i + 1}"
        # 将细胞ID和类型作为一个元组添加到列表中
        cells.append((cell_id, cell_type))
    return cells


def save_cells_to_txt(cells, filename="cells.txt"):
    """
    将细胞数据保存到txt文件中。

    :param cells: 一个列表，每个元素是一个包含细胞ID和类型的元组
    :param filename: 保存文件的名称
    """
    with open(filename, "w") as file:
        for cell_id, cell_type in cells:
            # 写入每一行，格式为"cell_i cell_type_i"
            file.write(f"{cell_id} {cell_type}\n")


def select_marker_genes(gene_expression, cells, num_marker_gene):
    # cells 是一个包含 (cell_id, cell_type) 元组的列表
    cell_types = np.array([cell[1] for cell in cells])
    cell_indices_by_type = {}
    for i, cell_type in enumerate(np.unique(cell_types)):
        cell_indices_by_type[cell_type] = np.where(cell_types == cell_type)[0]

    # 对每种细胞类型，找到高表达的基因
    selected_genes = set()
    marker_genes = {}
    for cell_type, indices in cell_indices_by_type.items():
        # 计算每种基因在该细胞类型中的平均表达量
        gene_mean_expression = np.mean(gene_expression[:, indices], axis=1)

        # 过滤掉已经被选为marker的基因
        available_genes = np.arange(gene_mean_expression.size)
        available_genes = available_genes[~np.isin(available_genes, selected_genes)]

        # 如果没有足够的新基因，则调整选择的基因数量
        num_to_select = min(num_marker_gene, len(available_genes))
        if num_to_select == 0:
            continue

            # 选择剩余基因中表达量最高的作为候选
        top_genes_indices = np.argsort(gene_mean_expression[available_genes])[-num_marker_gene:]
        top_genes = available_genes[top_genes_indices]

        # 更新已选基因集
        selected_genes.update(top_genes)

        # 存储marker基因
        marker_genes[cell_type] = top_genes

    # 修改这些基因在对应细胞类型中的表达量
    for cell_type, genes in marker_genes.items():
        cell_indices = cell_indices_by_type[cell_type]
        for gene in genes:
            gene_expression[gene, cell_indices] = np.random.randint(1000, 1501, size=cell_indices.shape)

    return gene_expression


# 修改combine_cell_gene函数以包含marker gene的选择
def combine_cell_gene_with_markers(num_cell, num_gene, num_cell_type, num_marker_gene):
    cells = generate_cell(num_cell, num_cell_type)
    cell_ids = [cell[0] for cell in cells]
    cell_types = [cell[1] for cell in cells]

    gene_expression = generate_gene_expression(num_gene, num_cell)
    gene_expression = select_marker_genes(gene_expression, cells, num_marker_gene)

    adata = ad.AnnData(X=gene_expression.T)
    adata.obs_names = cell_ids
    adata.obs['cell_type'] = cell_types
    adata.var_names = ['Gene_' + str(i + 1) for i in range(num_gene)]

    return adata


# 结果生成txt文件
def make_gene_txt(genes, filename='genes.txt'):
    with open(filename, 'w') as file:
        for gene in genes:
            file.write(f"{gene[0]}\t{gene[1]}\n")


# genes = generate_gene(num_gene=1000)
# make_gene_txt(genes)

# cells = generate_cell(100, 5)
# save_cells_to_txt(cells)

adata = combine_cell_gene_with_markers(100, 1000, 5, 10)
print(adata)
