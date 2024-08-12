import pandas as pd
import anndata as ad

count_txt = 'source/scRNA_count.txt'
cell_type_txt = 'source/scRNA_celltype.txt'


def combine_data(count_txt, cell_type_txt):
    """
    :param count_txt: 单细胞基因表达量文件    第一行为基因标签，第一列为细胞ID，其余为基因表达量
    :param celltype_txt: 细胞类型文件     第一列为细胞ID，第二列为细胞类型
    :return: 输出h5ad文件
    content:
        将两个文件内容整合，统计每种细胞类型的细胞数目。打印细胞数据结果，将整合的的数据
        保存为h5ad格式文件。
    """
    # 读取基因表达量文件
    cell_counts = pd.read_csv(count_txt, sep='\t')

    # 读取细胞类型注释文件
    cell_types = pd.read_csv(cell_type_txt, sep='\t', header=None, names=['cell_id', 'cell_type'])

    # 源文件中细胞ID为列，基因标签为行
    adata = ad.AnnData(X=cell_counts.T)  # 转置counts DataFrame

    # 将细胞类型信息添加到adata.obs中
    adata.obs['cell_type'] = cell_types.set_index('cell_id')['cell_type'].reindex(adata.obs_names).values

    # 统计每种细胞类型的细胞数目
    cell_type_counts = adata.obs['cell_type'].value_counts()
    print(cell_type_counts)

    # 保存为h5ad格式文件
    adata.write_h5ad('integrated_data.h5ad')

combine_data(count_txt, cell_type_txt)
