import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as pl
import math
import glob
import argparse
import rf_support as rfs
import os
import os.path


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
	parser.add_argument("--output_dir", default="extract", help="Directory to store outputs")
	args = parser.parse_args()
	
	z0_probe = complex(args.z0_real, args.z0_imag)
	(freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec) = extract_rlgc(args.pad_L_csv_file, args.pad_2L_csv_file, z0_probe, args.method, args.skip_plots, args.struct_csv_name, args.skip_deembed, args.tag, args.output_dir)

	

def extract_rlgc(pad_L_csv_filename, pad_2L_csv_filename, z0_probe=complex(50.0,0), method="distributed", skip_plots=False, struct_csv_name="*.csv", skip_deembed=False, output_tag = "", output_dir="extract"):

	file_list = glob.glob(struct_csv_name)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	
	print("Pad Deembedding file (L):  {0:s}".format(pad_L_csv_filename) )
	print("Pad Deembedding file (2L): {0:s}".format(pad_2L_csv_filename) )
	# Get pad deembedding parameters
	if not skip_deembed:
		(freq, abcd_pad, abcd_pad_inv, Sri_pad, Sdb_pad, Sdeg_pad) = get_pad_abcd(pad_L_csv_filename, pad_2L_csv_filename, z0_probe)
	else:
		abcd_pad_inv = [] # dummy value needed for extract_rlcg_from_measurement call below
	
	print("Extracting RLGC using {0:s} method...".format(method))
	freq_mat = []
	R_mat = []
	L_mat = []
	G_mat = []
	C_mat = []
	length_vec = []
	width_vec = []
	name_vec = []
	for filename in file_list:
		# First strip out any stuff from the path
		nfilename_arr = filename.split("\\")
		nfilename = nfilename_arr[-1]
		
		# now process the actual filename
		filename_arr = nfilename.split("_")
		trace_length_um = int(filename_arr[0])
		length_m = trace_length_um * 1e-6
		
		trace_width_name = filename_arr[1]
		trace_width_um = int( trace_width_name[0:-2] ) # get rid of the "um" in the width section
		data_final_arr = "_".join(filename_arr[2:]).split(".")
		data_final_str = data_final_arr[0]
		print("\tL: {0:d}um \t W: {1:d}um \t Sample: {2:s}".format(trace_length_um, trace_width_um, data_final_str) )
		
		# Construct output filename for each input file
		# Requires input files to be named as follows
		# $LENGTH_$WIDTHum_$WHATEVER.s2p
		# where $LENGTH is the structure length in microns
		# $WIDTH is the structure width in microns
		# and $WHATEVER is whatever's left over,  typically a sample number and maybe some other info
		structure_string = "L{0:d}um_W{1:d}um_{2:s}".format(trace_length_um, trace_width_um, data_final_str)
		rlgc_filename = "rlgc_" + structure_string + output_tag +  ".csv"

		(freq_hz, S, Z, T, Sdb, Sdeg) = rfs.get_rf_params_from_vna_csv(filename)
		Sdb_dut = Sdb
		Sdeg_dut = Sdeg
		abcd_dut = T
		
		(freq, R, L, G, C) = extract_rlcg_from_measurement( freq_hz, length_m, abcd_pad_inv, abcd_dut, z0_probe, method, skip_deembed)
		freq_mat.append(freq)
		R_mat.append(R)
		L_mat.append(L)
		G_mat.append(G)
		C_mat.append(C)
		name_vec.append(structure_string)
		length_vec.append(trace_length_um)
		width_vec.append(trace_width_um)

		write_rlgc(freq, R, L, G, C, rlgc_filename, output_dir)
		
		if not skip_plots:
			plot_rlgc(freq, R, L, G, C, structure_string + output_tag, output_dir)
			plot_s_params(freq, Sdb_dut, Sdeg_dut, structure_string + output_tag, output_dir)
	
	write_data(freq_mat[0], R_mat, name_vec, "R" + output_tag + ".csv", output_dir)
	write_data(freq_mat[0], L_mat, name_vec, "L" + output_tag + ".csv", output_dir)
	write_data(freq_mat[0], C_mat, name_vec, "C" + output_tag + ".csv", output_dir)
	write_data(freq_mat[0], G_mat, name_vec, "G" + output_tag + ".csv", output_dir)

	freq_min = 1e10
	freq_max = 2e20

	path_split = os.path.split(os.getcwd())
	cur_folder = path_split[-1]
	header_tag = cur_folder + output_tag
	write_averaged_data_freq_range(freq_mat[0], freq_min, freq_max, R_mat, length_vec, cur_folder + "_R" + output_tag + "_avg" + ".csv", output_dir, header_tag)

	output_dir = "C:\\Users\\William\\Dropbox\\Research\\Groups\\I3DS\\Projects\\Wire Measurements\\Will_Xuchen\\R_avg"
	write_averaged_data_freq_range(freq_mat[0], freq_min, freq_max, R_mat, length_vec, cur_folder + "_R" + output_tag + "_avg" + ".csv", output_dir, header_tag)
			
	return (freq_mat, R_mat, L_mat, G_mat, C_mat, name_vec, length_vec, width_vec)
	


def deembed_pads_from_measurement(abcd_pad_inv, abcd_dut, z0_probe = 50):
	# (abcd_dut_deembedded, Sri_dut, Sdb_dut, Sdeg_dut) = deembed_pads_from_measurement(abcd_pad_inv, abcd_dut, z0_probe = 50)
	
	num_freqs = len(abcd_dut)
	abcd_dut_deembedded = np.zeros( (num_freqs, 2, 2), dtype=complex )
	Sri_dut_deembedded = np.zeros( (num_freqs, 2, 2), dtype=complex )
	
	for idx, Pinv in enumerate(abcd_pad_inv):
		abcd_dut_deembedded_f = np.dot( Pinv, np.dot( abcd_dut[idx], Pinv) )
		abcd_dut_deembedded[idx] = abcd_dut_deembedded_f
		
	Sri_dut_deembedded = rfs.abcd2s(abcd_dut_deembedded, z0_probe, z0_probe)
	(Sdb_dut_deembedded, Sdeg_dut_deembedded) = rfs.sri2sdb(Sri_dut_deembedded)
	#Sri_dut_deembedded[idx] = Sri_dut_deembedded_f
		
	#Sri_dut = abcd2s(abcd_dut, z0_probe, z0_probe)
	#(Sdb_dut, Sdeg_dut) = sri2sdb(Sri_dut)
	
	return (abcd_dut_deembedded, Sri_dut_deembedded, Sdb_dut_deembedded, Sdeg_dut_deembedded)
	
	
	
def extract_rlcg_from_measurement( freq, length_m, abcd_pad_inv, abcd_meas, z0_probe = 50, method="distributed", skip_deembed=False):
	# (freq, R, L, G, C) = extract_rlcg_from_measurement( freq, length_m, abcd_pad_inv, abcd_dut, z0_probe = 50, method="distributed")
	# if skip_deembed = True then abcd_pad_inv is not used -- just pass an empty array (or whatever)
	
	if not skip_deembed:
		(abcd_dut, Sri_dut, Sdb_dut, Sdeg_dut) = deembed_pads_from_measurement(abcd_pad_inv, abcd_meas, z0_probe)
	else:
		Sri_dut = rfs.abcd2s(abcd_meas, z0_probe, z0_probe)
		(Sdb_dut, Sdeg_dut) = rfs.sri2sdb(Sri_dut)
	
	if method == "distributed":		
			(freq, R, L, G, C, gamma, attenuation, losstan, Zc) = distributed_rlgc_from_sdb(length_m, freq, Sdb_dut, Sdeg_dut, z0_probe)
	elif method == "lumped":
#		net_dut = rf.Network( f=freq*1e-9, s=Sri_dut, z0=z0_probe)
#		(freq, R, L, G, C, Zdiff, Ycomm, net) = lumped_rlgc_from_Network(net_dut, z0_probe)
		print("ERROR: NOT IMPLEMENTED")
	else:
		print("ERROR: NOT IMPLEMENTED")
		
	return (freq, R, L, G, C)
	
	
def distributed_rlgc_from_sdb(length_m, freq, Sdb, Sdeg, z0_probe=complex(50,0)):
	# length_m:	(m)	Length of structure being measured
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.
	
	Sri = rfs.sdb2sri(Sdb, Sdeg)
	abcd = rfs.s2abcd(Sri, z0_probe, z0_probe)
	
	d_vec = np.zeros( (len(abcd)), dtype=complex)
	c_vec = np.zeros( (len(abcd)), dtype=complex)
	
	for idx, T in enumerate(abcd):	
		d_vec[idx] = T[1][1]
		c_vec[idx] = T[1][0] # C vector (what I think needs to be used for Zc extraction)

	gamma = 1/length_m * np.arccosh(d_vec)
	Zc = 1/c_vec * np.sinh(gamma * length_m)

	R = ( gamma * Zc).real
	L = 1/2/math.pi/freq * ((gamma * Zc).imag )
	G = (gamma/Zc).real
	C = 1/2/math.pi/freq * (gamma/Zc).imag
	losstan = (gamma/Zc).real / (gamma/Zc).imag

	attenuation = 20*np.log10( abs(np.exp(-gamma*length_m)) )

	return ( freq, R, L, G, C, gamma, attenuation, losstan, Zc )

def lumped_rlgc_from_Network(net, z0_probe=complex(50,0)):
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.

	freq = net.f

	Z = net.z
	Y = net.y

	Zdiff = np.zeros( (len(Z)), dtype=complex)
	Ycomm = np.zeros( (len(Z)), dtype=complex)
	for idx in range(len(Zdiff)):
		zz = Z[idx]
		yy = Y[idx]
		
		Zdiff[idx] = z0_probe*(zz[0,0] - zz[0,1] - zz[1,0] + zz[1,1])
		Ycomm[idx] = yy[0,0] + yy[0,1] + yy[1,0] + yy[1,1]

	R = Zdiff.real
	L = 1/2/math.pi/freq *(Zdiff.imag)
	G = Ycomm.real
	C = 1/2/math.pi/freq * (Ycomm.imag)

	return (freq, R, L, G, C, Zdiff, Ycomm, net)


def get_pad_abcd(pad_L_s2p_filename, pad_2L_s2p_filename, z0_probe=complex(50,0)):
	
	(freq_L, Sdb_L, Sdeg_L) = rfs.get_sdb_from_vna_csv(pad_L_s2p_filename)
	(freq_2L, Sdb_2L, Sdeg_2L) = rfs.get_sdb_from_vna_csv(pad_2L_s2p_filename)

	# ABCD matrices
	S_L = rfs.sdb2sri(Sdb_L, Sdeg_L) 
	S_2L = rfs.sdb2sri(Sdb_2L, Sdeg_2L) 
	
	abcd_L = rfs.s2abcd( S_L, z0_probe)
	abcd_2L = rfs.s2abcd( S_2L, z0_probe)
	
	abcd_pad = np.zeros( np.shape(abcd_L), dtype=complex) # ABCD matrix structure for pad
	abcd_pad_inv = np.zeros( np.shape(abcd_L), dtype=complex) # inverse abcd matrix structure for pad
	
	# iterating across each frequency point
	for idx, abcd_L_mat in enumerate(abcd_L):
		abcd_2L_mat = abcd_2L[idx]
		
		abcd_L_inv = la.inv(abcd_L_mat)
		abcd_P_squared = la.inv( np.dot( abcd_L_inv, np.dot( abcd_2L_mat, abcd_L_inv) ) ) # PP = ( ML^-1 * M2L * ML^-1 )^-1
		abcd_P = la.sqrtm(abcd_P_squared) # ABCD matrix of the pad (single pad) for this frequency
		abcd_P_inv = la.inv(abcd_P)
		
		abcd_pad[idx] = abcd_P
		abcd_pad_inv[idx] = abcd_P_inv
	
	Sri_pad = rfs.abcd2s(abcd_pad, z0_probe, z0_probe)
	(Sdb_pad, Sdeg_pad) = rfs.sri2sdb(Sri_pad)
	
	return (freq_L, abcd_pad, abcd_pad_inv, Sri_pad, Sdb_pad, Sdeg_pad)


def write_s_db_deg( sdb, sdeg, freq, filename):
	outfile = open(filename,'w')

	freq_vec = freq
	s_db_list = sdb
	s_deg_list = sdeg

	for idx in range(len(freq_vec)):
		freq = freq_vec[idx]
		s_db_mat = s_db_list[idx]
		S11_db = s_db_mat[0][0]
		S12_db = s_db_mat[0][1]
		S21_db = s_db_mat[1][0]
		S22_db = s_db_mat[1][1]

		s_deg_mat = s_deg_list[idx]
		S11_deg = s_deg_mat[0][0]
		S12_deg = s_deg_mat[0][1]
		S21_deg = s_deg_mat[1][0]
		S22_deg = s_deg_mat[1][1]

		outstr = "{0:.8g},{1:.8g},{2:.8g},{3:.8g},{4:.8g},{5:.8g},{6:.8g},{7:.8g},{8:.8g}\n".format( freq, S11_db, S11_deg, S12_db, S12_deg, S21_db, S21_deg, S22_db, S22_deg)

		outfile.write(outstr)


def write_rlgc(freq, R, L, G, C, filename, output_dir=""):
	filename = os.path.join(output_dir, filename)
	outfile = open(filename, 'w')
	
	for idx, f in enumerate(freq):
		r = R[idx]
		l = L[idx]
		g = G[idx]
		c = C[idx]
		
		outstr = "{0:.8g},{1:.8g},{2:.8g},{3:.8g},{4:.8g}\n".format(f, r, l, g, c)
		outfile.write(outstr)
		

def write_data( freq, data_mat, name_mat, filename, output_dir=""):

	filename = os.path.join(output_dir, filename)
	outfile = open(filename, 'w')
	
	outstr = ",".join(name_mat)
	outfile.write("{0:s},{1:s}\n".format("Freq (Hz)", outstr) )
	
	data_mat = np.array(data_mat)

	for idx, f in enumerate(freq):
		data_str_vec = [ "{0:.8g}".format(el) for el in data_mat[:,idx] ]
		data_str = ",".join(data_str_vec)
		outfile.write("{0:.8g},{1:s}\n".format(f, data_str) )

def write_averaged_data_freq_range(freq, freq_min, freq_max, data_mat, length_vec, filename, output_dir="", header_tag=""):

	filename = os.path.join(output_dir, filename)
	outfile = open(filename, 'w')

	data_avg_vec = np.zeros( (len(length_vec)), dtype=float)

	for idx, length in enumerate(length_vec):
		data_avg = np.average( data_mat[idx][ (freq >=freq_min) & (freq <=freq_max) ] )
		data_avg_vec[idx] = data_avg

	sort_inds = np.argsort(length_vec)
	outfile.write("Length (um),{0:s}\n".format(header_tag))

	for idx in sort_inds:
		outfile.write("{0:d},{1:.8g}\n".format(length_vec[idx], data_avg_vec[idx]) )

		

def plot_rlgc(freq, R, L, G, C, structure_string, output_dir=""):
	freq_ghz = freq/1e9
	
	pl.figure(1, figsize=(9,13) )
	pl.clf()
	ax1 = pl.subplot(4,1,1)
	pl.plot(freq_ghz, R, "b", linewidth=2)
	pl.xlabel("Frequency (GHz)")
	#pl.ylabel("Resistance (Ohms)")
	pl.ylabel("R ($\Omega$/m)")
	ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	
	ax2 = pl.subplot(4,1,2)
	pl.plot(freq_ghz, L, "b", linewidth=2)
	pl.xlabel("Frequency (GHz)")
	#pl.ylabel("Inductance (H)")
	pl.ylabel("L (H/m)")
	ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	
	ax3 = pl.subplot(4,1,3)
	pl.plot(freq_ghz, G, "b", linewidth=2)
	pl.xlabel("Frequency (GHz)")
	#pl.ylabel("Conductance (S)")
	pl.ylabel("G (S/m)")
	ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	
	ax4 = pl.subplot(4,1,4)
	pl.plot(freq_ghz, C, "b", linewidth=2)
	pl.xlabel("Frequency (GHz)")
	#pl.ylabel("Capacitance (F)")
	pl.ylabel("C (F/m)")
	ax4.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	
	filename = structure_string + "_RLGC.pdf"
	filename = os.path.join(output_dir, filename)
	pl.savefig()
	
	
	
def plot_s_params(freq, Sdb, Sdeg, structure_string, output_dir=""):
	freq_ghz = freq/1e9
	
	S11_db = np.zeros( (len(Sdb)) )
	S12_db = np.zeros( (len(Sdb)) )
	S21_db = np.zeros( (len(Sdb)) )
	S22_db = np.zeros( (len(Sdb)) )
	
	S11_deg = np.zeros( (len(Sdb)) )
	S12_deg = np.zeros( (len(Sdb)) )
	S21_deg = np.zeros( (len(Sdb)) )
	S22_deg = np.zeros( (len(Sdb)) )
	
	for idx in range(len(Sdb)):
		S11_db[idx] = Sdb[idx][0][0]
		S12_db[idx] = Sdb[idx][0][1]
		S21_db[idx] = Sdb[idx][1][0]
		S22_db[idx] = Sdb[idx][1][1]
		
		S11_deg[idx] = Sdeg[idx][0][0]
		S12_deg[idx] = Sdeg[idx][0][1]
		S21_deg[idx] = Sdeg[idx][1][0]
		S22_deg[idx] = Sdeg[idx][1][1]
		
	pl.figure(1, figsize=(8.5,11) )
	pl.clf()
	pl.subplot(2,1,1)
	pl.hold(True)
	pl.plot(freq_ghz, S11_db, 'b', linewidth=2, label="S11")
	pl.plot(freq_ghz, S22_db, 'g', linewidth=2, label="S22")
	pl.xlabel("Frequency (GHz)")
	pl.ylabel("S Parameters (DB)")
	pl.grid()
	pl.legend()
	
	pl.subplot(2,1,2)
	pl.hold(True)
	pl.plot(freq_ghz, S12_db, 'b', linewidth=2, label="S12")
	pl.plot(freq_ghz, S21_db, 'g', linewidth=2, label="S21")
	pl.xlabel("Frequency (GHz)")
	pl.ylabel("S Parameters (DB)")
	pl.grid()
	pl.legend()
	
	filename = structure_string + "_Sdb.pdf"
	filename = os.path.join(output_dir, filename)
	pl.savefig(filename)
	
	pl.figure(2)
	pl.clf()
	pl.subplot(2,1,1)
	pl.hold(True)
	pl.plot(freq_ghz, S11_deg, 'b', linewidth=2, label="S11")
	pl.plot(freq_ghz, S22_deg, 'g', linewidth=2, label="S22")
	pl.xlabel("Frequency (GHz)")
	pl.ylabel("S Parameter Phase (Degrees)")
	pl.grid()
	pl.legend()
	
	pl.subplot(2,1,2)
	pl.hold(True)
	pl.plot(freq_ghz, S12_deg, 'b', linewidth=2, label="S12")
	pl.plot(freq_ghz, S21_deg, 'g', linewidth=2, label="S21")
	pl.xlabel("Frequency (GHz)")
	pl.ylabel("S Parameter Phase (Degrees)")
	pl.grid()
	pl.legend()
	
	filename = structure_string + "_Sdeg.pdf"
	filename = os.path.join(output_dir, filename)
	pl.savefig(filename)
	
		

if (__name__ == "__main__"):
	main()
