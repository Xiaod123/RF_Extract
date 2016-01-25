# Simple test script
import matplotlib.pyplot as pl
import extraction as ex

#pl.ion() # interactive plotting on to avoid blocking

def create_plot(freq_mat, data_mat, length_vec, color_keys, plot_name):
	print("Creating plot: {0:s}".format(plot_name) )
	pl.figure(1)
	pl.clf()
	pl.hold(True)
	
	for idx, freq in enumerate(freq_mat):
		data = data_mat[idx]
		length = length_vec[idx]
		
		pl.semilogy(freq, data, color_keys[length])
		pl.semilogy(freq, -data, color_keys[length] + ":")
		
		pl.savefig(plot_name)

#pad_l_s2p_file= "test/500_5um_4.s2p"
#pad_2l_s2p_file= "test/1000_5um_4B.s2p"

base_dir = "test"
base_dir = base_dir + "/"

pad_l_s2p_file= "500_3um_1.csv"
pad_2l_s2p_file= "1000_3um_1.csv"

pad_l_s2p_file= base_dir + pad_l_s2p_file
pad_2l_s2p_file= base_dir + pad_2l_s2p_file

z0= 50
method_extract = "distributed"

freq_mat = []
R_mat = []
L_mat = []
G_mat = []
C_mat = []

color_keys = { 2000:'k', 1000:'b', 500:'c', 100:'g', 50:'r' }

(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_csv_name=base_dir + "*_3um_*.csv", output_tag="_3um")
create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_3um.pdf")

(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_csv_name=base_dir + "*_3um_*.csv", skip_deembed=True, output_tag="_3um_no_deembed")
create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_3um_no_deembed.pdf")

#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name=base_dir + "*_5um_*.s2p", output_tag="_5um")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_5um.pdf")
#
#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name=base_dir + "*.s2p", output_tag="_all")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_all.pdf")
#
#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name=base_dir + "*_3um_*.s2p", skip_deembed=True, output_tag="_3um_no_deembed")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_3um_no_deembed.pdf")
#
#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name=base_dir + "*_5um_*.s2p", skip_deembed=True, output_tag="_5um_no_deembed")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_5um_no_deembed.pdf")
#
#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name=base_dir + "*.s2p", skip_deembed=True, output_tag="_all_no_deembed")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_all_no_deembed.pdf")
#
#	
#	
#
