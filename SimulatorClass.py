import combine_data as cd
import simulate as sl
import bin_operation as bo


class Simulator:

    def __init__(self):
        print("调用构造方法")

    def combine_data(self, count_file, type_file):
        adata = cd.combine_data(count_file, type_file)
        return adata

    def generate_h5ad_file(self, adata, h5ad_path):
        cd.generate_h5ad_file(adata, h5ad_path)

    def simulate(self, h5ad_path, chip_size, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm, gene_density):
        spatial_df = sl.simulate(h5ad_path, chip_size, cell_density_range_mm, cell_diameter_range_um, spot_spacing_nm, gene_density)
        return spatial_df

    def generate_gem(self, spatial_df, gem_path):
        sl.generate_gem(spatial_df, gem_path)

    def draw_picture(self, spatial_df):
        sl.draw_picture(spatial_df)

    def bin_operation(self, gem_path, bin_size):
        bin_adata = bo.bin_operation(gem_path, bin_size)
        return bin_adata

    def gengerate_bin_h5ad_file(self, bin_adata, gengerate_bin_h5ad_file):
        bo.gengerate_bin_h5ad_file(bin_adata, gengerate_bin_h5ad_file)


if __name__ == '__main__':
    simulator = Simulator()
    adata = simulator.combine_data('source/scRNA_count.txt', 'source/scRNA_celltype.txt')
    simulator.generate_h5ad_file(adata, 'integrated_data.h5ad')
    spatial_df = simulator.simulate('integrated_data.h5ad', 1, (50, 100), (10, 20), 500, (500, 1000))
    simulator.generate_gem(spatial_df, 'simulated_spatial_transcriptomics.gem')
    simulator.draw_picture(spatial_df)
    bin_adata = simulator.bin_operation('simulated_spatial_transcriptomics.gem', 10000)
    simulator.gengerate_bin_h5ad_file(bin_adata, 'binned_spatial_transcriptomics.h5ad')
