import extraction as ex
import argparse
import os
import os.path
import matplotlib.pyplot as pl


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("pad_L_csv_file", help="Filename for L structure measurement to be used for pad extraction")
	parser.add_argument("pad_2L_csv_file", help="Filename for 2L structure measurement to be used for pad extraction")
	parser.add_argument("--struct_csv_name", default="*.csv", help="Filename for structure to convert. If this argument is presented, ONLY file names conforming to this naming scheme will be processed. Accepts globs (i.e. input *_foo.s2p to process all files ending with _foo.s2p). Default is *.s2p (all s2p files)")
	parser.add_argument("--skip_deembed", default=False, action='store_true', help="Use this flag to skip pad deembedding. You will still need to input the pad L/2L filenames, but they will not be used")	
	parser.add_argument("--z0_real", type=float, default=50, help="Real portion of probe impedance. Default is 50 Ohms")
	parser.add_argument("--z0_imag", type=float, default=0, help="Imaginary portion of probe impedance. Default is 0 Ohms (Default impedance is 50 + 0j)")
	parser.add_argument("--skip_plots", action="store_true", default=False, help="Skip plotting for faster data extraction")
	parser.add_argument("--method", default="distributed", choices=["distributed", "lumped"], help="Type of RLGC extraction to perform. distributed (default) -- treats structure as transmission line and extracts from S in DB/DEG form. lumped -- treats structure as lumped element.") 
	parser.add_argument("--tag", default="", help="Output file tag")
	args = parser.parse_args()
	
	z0_probe = complex(args.z0_real, args.z0_imag)

	color_keys = { 2000:'k', 1000:'b', 500:'c', 100:'g', 50:'r' }

	curdir = os.path.split(os.getcwd())
	outdir = "extract_" + curdir[1]

	(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(args.pad_L_csv_file, args.pad_2L_csv_file, z0_probe, args.method, skip_plots=True, struct_csv_name="*_3um_*.csv", skip_deembed=True, output_tag="_3um_no_deembed", output_dir=outdir + "_3um_no_deembed")
	create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_3um_no_deembed.pdf", output_dir = outdir + "_3um_no_deembed")

	(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(args.pad_L_csv_file, args.pad_2L_csv_file, z0_probe, args.method, skip_plots=True, struct_csv_name="*_5um_*.csv", skip_deembed=True, output_tag="_5um_no_deembed", output_dir=outdir + "_5um_no_deembed")
	create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_5um_no_deembed.pdf", output_dir = outdir + "_5um_no_deembed")

	(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(args.pad_L_csv_file, args.pad_2L_csv_file, z0_probe, args.method, skip_plots=True, struct_csv_name="*_3um_*.csv", skip_deembed=False, output_tag="_3um", output_dir=outdir + "_3um")
	create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_3um.pdf", output_dir = outdir + "_3um")

	(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = ex.extract_rlgc(args.pad_L_csv_file, args.pad_2L_csv_file, z0_probe, args.method, skip_plots=True, struct_csv_name="*_5um_*.csv", skip_deembed=False, output_tag="_5um", output_dir=outdir + "_5um")
	create_plot(freq_mat, R_mat, length_vec, color_keys, "R_plot_5um.pdf", output_dir = outdir + "_5um")
	


def create_plot(freq_mat, data_mat, length_vec, color_keys, plot_name, output_dir=""):
	print("Creating plot: {0:s}".format(plot_name) )
	pl.figure(1)
	pl.clf()
	pl.hold(True)
	
	for idx, freq in enumerate(freq_mat):
		data = data_mat[idx]
		length = length_vec[idx]
		
		num_pos = 0
		num_neg = 0
		for el in data:
			if el > 0:
				num_pos += 1
			if el < 0:
				num_neg +=1
		if num_pos > 0:
			pl.semilogy(freq/1e9, data, color_keys[length])
		if num_neg > 0:
			pl.semilogy(freq/1e9, -data, color_keys[length] + ":")
		pl.xlabel('Frequency (GHz)')
		pl.ylabel("PUL R ($\Omega$/m)")
		
	output_name = os.path.join(output_dir, plot_name)
	pl.savefig(output_name)


if (__name__ == "__main__"):
	main()
