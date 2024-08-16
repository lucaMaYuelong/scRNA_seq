import argparse
import combine_data as cd
import simulate as sl
import bin_operation as bo
import anndata as ad

# 使用argparse 通过命令行设置参数
def parse_args():
    parse = argparse.ArgumentParser(description='Simulation of combining cell data')
    parse.add_argument('-count', '--count_data', help='gene count file, in txt format')
    parse.add_argument('-type', '--celltype_data', help='cell type file, in txt format')
    parse.add_argument('-hfp', '--h5ad_file_path', default='integrated_data.h5ad', help='path of combine data in h5ad format')
    parse.add_argument('-gfp', '--gem_file_path', default='simulated_spatial_transcriptomics.gem', help='path of simulation '
                                                                                              'data in gem format')
    parse.add_argument('-bhfp', '--bin_h5ad_file_path', default='simulated_spatial_transcriptomics.h5ad', help='path of bin '
                                                                                                      'operation data '
                                                                                                      'in h5ad format')
    parse.add_argument('-cs', '--chip_size_mm', type=int, help='size of chip, unit is mm')
    parse.add_argument('-cnmin', '--cell_density_range_min_mm', type=int, help='minimum number of cells, unit is mm')
    parse.add_argument('-cnmax', '--cell_density_range_max_mm', type=int, help='maximum number of cells, unit is mm')
    parse.add_argument('-csmin', '--cell_diameter_range_min_um', type=int, help='minimum size of cells, unit is um')
    parse.add_argument('-csmax', '--cell_diameter_range_max_um', type=int, help='maximum size of cells, unit is um')
    parse.add_argument('-ss', '--spot_spacing_nm', type=int, help='spacing of adjacent spot, unit is nm')
    parse.add_argument('-gdmin', '--gene_density_min', type=int, help='minimum number of genes in each cells')
    parse.add_argument('-gdmax', '--gene_density_max', type=int, help='maximum number of genes in each cells')
    parse.add_argument('-bs', '--bin_size', type=int, help='size of bin')
    args = parse.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    count_data = args.count_data
    celltype_data = args.celltype_data
    h5ad_path = args.h5ad_file_path
    gem_path = args.gem_file_path
    bin_h5ad_path = args.bin_h5ad_file_path
    chip_size = args.chip_size_mm
    cell_density_range_mm = (args.cell_density_range_min_mm, args.cell_density_range_max_mm)
    cell_diameter_range_um = (args.cell_diameter_range_min_um, args.cell_diameter_range_max_um)
    spot_spacing_um = args.spot_spacing_nm
    gene_density = (args.gene_density_min, args.gene_density_max)
    bin_size = args.bin_size

    # combine data
    adata = cd.combine_data(count_data, celltype_data)

    # generate h5ad file
    cd.generate_h5ad_file(adata, h5ad_path)

    # simulate
    spatial_df = sl.simulate(h5ad_path, chip_size, cell_density_range_mm, cell_diameter_range_um, spot_spacing_um, gene_density)

    # generate gem file
    sl.generate_gem(spatial_df, gem_path)

    # generate image
    sl.draw_picture(spatial_df)

    # binning operation
    bin_adata = bo.bin_operation(gem_path, bin_size)

    # generate bin data in h5ad format
    bo.gengerate_bin_h5ad_file(bin_adata, bin_h5ad_path)

# -count 'scRNA_count.txt' -type 'scRNA_celltype.txt' -hfp 'integrated_data.h5ad' -gfp 'simulated_spatial_transcriptomics.gem' -bhfp 'binned_spatial_transcriptomics.h5ad' -cs 1 -cnmin 50 -cnmax 100 -csmin 10 -csmax 20 -ss 500 -gdmin 500 -gdmax 1000 -bs 10000