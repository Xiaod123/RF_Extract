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

pad_l_s2p_file= "test/500_5um_4.s2p"
pad_2l_s2p_file= "test/1000_5um_4B.s2p"

pad_l_s2p_file="C:\Users\William\Dropbox\Research\Groups\I3DS\Projects\Wire Measurements\Will_Xuchen\W100\500_3um_1b.s2p"
pad_2l_s2p_file="C:\Users\William\Dropbox\Research\Groups\I3DS\Projects\Wire Measurements\Will_Xuchen\W100\1000_3um_1b.s2p"

z0= 50
method_extract = "distributed"

freq_mat = []
R_mat = []
L_mat = []
G_mat = []
C_mat = []

color_keys = { 2000:'k', 1000:'b', 500:'c', 100:'g', 50:'r' }

(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name="test/*_3um_*.s2p")
create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_d500-5.pdf")
#create_plot(freq_mat, L_mat, length_vec, color_keys, "L_plot.pdf")
#create_plot(freq_mat, G_mat, length_vec, color_keys, "G_plot.pdf")
#create_plot(freq_mat, C_mat, length_vec, color_keys, "C_plot.pdf")

(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name="test/*_3um_*.s2p", skip_deembed=True)
create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_no_deembed_d500-5.pdf")

#pad_l_s2p_file= "test/50_3um_1.s2p"
#pad_2l_s2p_file= "test/100_3um_4.s2p"
#
#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name="test/*_3um_*.s2p")
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_d50-3.pdf")

#(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(pad_l_s2p_file, pad_2l_s2p_file, z0_probe=z0, method=method_extract, skip_plots=True, struct_s2p_name="test/*_3um_*.s2p", skip_deembed=True)
#create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_no_deembed_d50-3.pdf")


	
	
	
	

