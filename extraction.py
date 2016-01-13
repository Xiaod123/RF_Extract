import scipy as sp
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib as pl
import skrf as rf # RF functions. To install this do "conda install -c scikit-rf  scikit-rf" from the command line if you have Anaconda Python installed, otherwise do "pip install scikit-rf"
import math


def test():
	pad_L_s2p_filename = "lpad.s2p"
	pad_2L_s2p_filename = "2lpad.s2p"
	structure_L_s2p_filename = "ltsv.s2p"
	structure_2L_s2p_filename = "2ltsv.s2p"
	
	length_m = 300e-6
	
	(freq, Sri_L, Sri_2L, abcd_L, abcd_2L, Sdb_L, Sdeg_L, Sdb_2L, Sdeg_2L, net_L, net_2L) = l2l_deembed(pad_L_s2p_filename, pad_2L_s2p_filename, structure_L_s2p_filename, structure_2L_s2p_filename)
	(freq, R, L, G, C, gamma, attenuation, losstan, Zc) = distributed_rlgc_from_abcd(length_m, freq, abcd_L)
	write_s_db_deg(Sdb_L, Sdeg_L, freq, "sdb_sdeg.csv")
	write_rlgc(freq, R, L, G, C, "rlgc.csv")
	write_rlgc(freq, gamma.real, gamma.imag, Zc.real, Zc.imag, "gamma.csv")
	
	(freq, R, L, G, C, gamma, attenuation, losstan, Zc) = distributed_rlgc_from_sdb(length_m, freq, Sdb_L, Sdeg_L)
	write_rlgc(freq, R, L, G, C, "rlgc2.csv")
	write_rlgc(freq, gamma.real, gamma.imag, Zc.real, Zc.imag, "gamma2.csv")
	
	(freq, R, L, G, C, gamma, attenuation, losstan, Zc) = distributed_rlgc_from_sri(length_m, freq, Sri_L)
	write_rlgc(freq, R, L, G, C, "rlgc_sri.csv")
	write_rlgc(freq, gamma.real, gamma.imag, Zc.real, Zc.imag, "gamma_sri.csv")
	
	(freq, R, L, G, C, Zdiff, Ycomm, net) = lumped_rlgc_from_Network(net_L)
	write_rlgc(freq, R, L, G, C, "rlgc_lumped.csv")



def distributed_rlgc_from_abcd(length_m, freq, abcd_mat_array, z0_probe=50):
	# length_m:	(m)	Length of structure being measured
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.

	gamma = np.zeros( (len(abcd_mat_array)), dtype=complex)
	Zc = np.zeros( (len(abcd_mat_array)), dtype=complex)
	d_vec = np.zeros( (len(abcd_mat_array)), dtype=complex)
	b_vec = np.zeros( (len(abcd_mat_array)), dtype=complex)
	
	R = np.zeros( (len(abcd_mat_array)) )
	L = np.zeros( (len(abcd_mat_array)) )
	G = np.zeros( (len(abcd_mat_array)) )
	C = np.zeros( (len(abcd_mat_array)) )
	losstan = np.zeros( (len(abcd_mat_array)) )
	attenuation = np.zeros( (len(abcd_mat_array)) )
	
	for idx in range(len(abcd_mat_array)):
		abcd = abcd_mat_array[idx]
		print(freq[idx])
		
		d_vec[idx] = abcd[1][1]
		b_vec[idx] = abcd[1][0] # this is actually the C from abcd
		
		gamma[idx] = 1/length_m*np.arccosh(d_vec[idx])
		Zc[idx] = b_vec[idx]**(-1) * np.sinh(gamma[idx] * length_m)
		
		R[idx] = (gamma[idx] * Zc[idx]).real
		L[idx] = 1/2/math.pi/freq[idx] * ( (gamma[idx]*Zc[idx]).imag)
		G[idx] = (gamma[idx] / Zc[idx]).real
		C[idx] = 1/2/math.pi/freq[idx] * ( (gamma[idx]/Zc[idx]).imag)
		
		losstan[idx] = ( (gamma[idx]/Zc[idx]).real) / ( (gamma[idx]/Zc[idx]).imag)
		attenuation[idx] = 20*np.log10( abs(np.exp( -gamma[idx] * length_m )) )


#	R = ( gamma * Zc).real
#	L = 1/2/math.pi/freq * ((gamma * Zc).imag )
#	G = (gamma/Zc).real
#	C = 1/2/math.pi/freq * (gamma/Zc).imag
#	losstan = (gamma/Zc).real / (gamma/Zc).imag
#
#	attenuation = 20*np.log10( abs(np.exp(-gamma*length_m)) )

	return ( freq, R, L, G, C, gamma, attenuation, losstan, Zc )
	
def distributed_rlgc_from_sdb(length_m, freq, Sdb, Sdeg, z0_probe=50):
	# length_m:	(m)	Length of structure being measured
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.
	
	Sri_mat = np.zeros( (len(Sdb),2,2), dtype=complex)
	
	for idx in range(len(Sdb)):
		Sdb_mat = Sdb[idx]
		Sdeg_mat = Sdeg[idx]
		
		Sri_mat[idx][0][0] = 10**(Sdb_mat[0][0]/20) * ( np.cos(Sdeg_mat[0][0]*np.pi/180) + 1j*np.sin(Sdeg_mat[0][0]*np.pi/180) )
		Sri_mat[idx][0][1] = 10**(Sdb_mat[0][1]/20) * ( np.cos(Sdeg_mat[0][1]*np.pi/180) + 1j*np.sin(Sdeg_mat[0][1]*np.pi/180) )
		Sri_mat[idx][1][0] = 10**(Sdb_mat[1][0]/20) * ( np.cos(Sdeg_mat[1][0]*np.pi/180) + 1j*np.sin(Sdeg_mat[1][0]*np.pi/180) )
		Sri_mat[idx][1][1] = 10**(Sdb_mat[1][1]/20) * ( np.cos(Sdeg_mat[1][1]*np.pi/180) + 1j*np.sin(Sdeg_mat[1][1]*np.pi/180) )
		
	abcd_mat_array = sri2abcd(Sri_mat)
	
	

	gamma = np.zeros( (len(abcd_mat_array)), dtype=complex)
	Zc = np.zeros( (len(abcd_mat_array)), dtype=complex)
	d_vec = np.zeros( (len(abcd_mat_array)), dtype=complex)
	b_vec = np.zeros( (len(abcd_mat_array)), dtype=complex)
	
	#R = np.zeros( (len(abcd_mat_array)) )
	#L = np.zeros( (len(abcd_mat_array)) )
	#G = np.zeros( (len(abcd_mat_array)) )
	#C = np.zeros( (len(abcd_mat_array)) )
	#losstan = np.zeros( (len(abcd_mat_array)) )
	#attenuation = np.zeros( (len(abcd_mat_array)) )
	
	for idx in range(len(abcd_mat_array)):
		abcd = abcd_mat_array[idx]
		
		d_vec[idx] = abcd[1][1]
		b_vec[idx] = abcd[1][0] # this is actually the C from abcd
		gamma[idx] = 1/length_m*np.arccosh(d_vec[idx])
		Zc[idx] = np.sinh(gamma[idx] * length_m)/b_vec[idx]
		#R[idx] = (gamma[idx] * Zc[idx]).real
		#L[idx] = 1/2/math.pi/freq[idx] * ( (gamma[idx] * Zc[idx]).imag)
		#G[idx] = (gamma[idx] / Zc[idx]).real
		#C[idx] = 1/2/math.pi/freq[idx] * ((gamma[idx]/Zc[idx]).imag)
		
		#losstan[idx] = ( (gamma[idx]/Zc[idx]).real) / ( (gamma[idx]/Zc[idx]).imag)
		#attenuation[idx] = 20*np.log10( abs(np.exp( -gamma[idx] * length_m )) )
		
		#print(gamma[idx] * Zc[idx])
		#print(1/2/math.pi/freq[idx] * ( (gamma[idx] * Zc[idx]).imag))


	R = ( gamma * Zc).real
	L = 1/2/math.pi/freq * ((gamma * Zc).imag )
	G = (gamma/Zc).real
	C = 1/2/math.pi/freq * (gamma/Zc).imag
	losstan = (gamma/Zc).real / (gamma/Zc).imag

	attenuation = 20*np.log10( abs(np.exp(-gamma*length_m)) )

	return ( freq, R, L, G, C, gamma, attenuation, losstan, Zc )
	
def distributed_rlgc_from_sri(length_m, freq, Sri, z0_probe=50):
	# length_m:	(m)	Length of structure being measured
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.

	gamma = np.zeros( (len(Sri)), dtype=complex)
	Zc =    np.zeros( (len(Sri)), dtype=complex)
	d_vec = np.zeros( (len(Sri)), dtype=complex)
	b_vec = np.zeros( (len(Sri)), dtype=complex)
	
	
	for idx in range(len(Sri)):
		
		S11 = Sri[idx][0][0]
		S12 = Sri[idx][0][1]
		S21 = Sri[idx][1][0]
		S22 = Sri[idx][1][1]
		
		K = np.sqrt( ((S11**2 - S21**2 + 1)**2 - (2*S11)**2)/(2*S21**2))
		alpha = ( (1 - S11**2 + S21**2)/(2*S21) + K)**(-1)
		gamma[idx] = -1/length_m * np.log(alpha)
		
		Zc[idx] = np.sqrt( z0_probe**2 * ( (1+S11)**2 - S21**2)/( (1-S11)**2 - S21**2) )
		


	R = ( gamma * Zc).real
	L = 1/2/math.pi/freq * ((gamma * Zc).imag )
	G = (gamma/Zc).real
	C = 1/2/math.pi/freq * (gamma/Zc).imag
	losstan = (gamma/Zc).real / (gamma/Zc).imag

	attenuation = 20*np.log10( abs(np.exp(-gamma*length_m)) )

	return ( freq, R, L, G, C, gamma, attenuation, losstan, Zc )


def lumped_rlgc_from_Network(net, z0_probe=50):
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


def abcd2s(abcd_struct, Z01, Z02):
	# convert ABCD matrix to S matrix in real/imag format

	R01 = Z01.real
	R02 = Z02.real
	num_freqs = len(abcd_struct)
	S = np.zeros( (num_freqs, 2, 2), dtype=complex )
	for idx in range(len(abcd_struct)):
		mat = abcd_struct[idx]
		A = mat[0][0]
		B = mat[0][1]
		C = mat[1][0]
		D = mat[1][1]

		denom = (A*Z02 + B + C*Z01*Z02 + D*Z01)

		S11 = ( A*Z02 + B - C*np.conj(Z01)*Z02 - D*np.conj(Z01) ) / denom
		S12 = ( 2*(A*D - B*C)*np.sqrt(R01*R02) ) / denom
		S21 = ( 2*np.sqrt(R01*R02) ) / denom
		S22 = (-A*np.conj(Z02) + B - C*Z01*np.conj(Z02) + D*Z01 ) / denom

		S[idx][0][0] = S11
		S[idx][0][1] = S12
		S[idx][1][0] = S21
		S[idx][1][1] = S22


	return S

def abcd2s_alt(abcd_struct, Z01, Z02):
	# convert ABCD matrix to S matrix in real/imag format

	#R01 = Z01.real
	#R02 = Z02.real
	Z0 = Z01

	num_freqs = len(abcd_struct)
	S = np.zeros( (num_freqs, 2, 2), dtype=complex )
	for idx in range(len(abcd_struct)):
		mat = abcd_struct[idx]
		A = mat[0][0]
		B = mat[0][1]
		C = mat[1][0]
		D = mat[1][1]

		denom = (A + B/Z0 + C*Z0 + D)

		S11 = ( A + B/Z0 - C*Z0 - D ) / denom
		S12 = ( 2*(A*D - B*C) ) / denom
		S21 = ( 2 ) / denom
		S22 = (-A + B/Z0 - C*Z0 + D ) / denom

		S[idx][0][0] = S11
		S[idx][0][1] = S12
		S[idx][1][0] = S21
		S[idx][1][1] = S22


	return S
	

def sri2abcd(s_struct, Z01=50, Z02=50):
	# Convert Sparams in Real/Imag format to ABCD matrix
	R01 = Z01.real
	R02 = Z02.real
	num_freqs = len(s_struct)
	abcd = np.zeros( (num_freqs, 2, 2), dtype=complex )
	for idx in range(len(s_struct)):
		mat = s_struct[idx]
		S11 = mat[0][0]
		S12 = mat[0][1]
		S21 = mat[1][0]
		S22 = mat[1][1]

		denom = 2*S21*np.sqrt(R01*R02)
		
		
		A = ( (np.conj(Z01) + S11*Z01)*(1-S22)+S12*S21*Z01 ) / denom
		B = ( (np.conj(Z01) + S11*Z01)*(np.conj(Z02)+S22*Z02)-S12*S21*Z01*Z02 ) / denom
		C = ( (1-S11)*(1-S22)-S12*S21 ) / denom
		D = ( (1-S11)*(np.conj(Z02)+S22*Z02) + S12*S21*Z02 ) / denom

		abcd[idx][0][0] = A
		abcd[idx][0][1] = B
		abcd[idx][1][0] = C
		abcd[idx][1][1] = D

	return abcd
	
def sri2abcd_alt(s_struct, Z01=50, Z02=50):
	# Convert Sparams in Real/Imag format to ABCD matrix
	R01 = Z01.real
	R02 = Z02.real
	Z0 = Z01
	num_freqs = len(s_struct)
	abcd = np.zeros( (num_freqs, 2, 2), dtype=complex )
	for idx in range(len(s_struct)):
		mat = s_struct[idx]
		S11 = mat[0][0]
		S12 = mat[0][1]
		S21 = mat[1][0]
		S22 = mat[1][1]
		
		
		A = ( (1+S11)*(1-S22) + S12*S21 ) / (2*S21)
		B = Z0*( ((1+S11)*(1+S22) - S12*S21) / (2*S21) )
		C = ( (1-S11)*(1-S22) - S12*S21) / (2*S21*Z0) 
		D = ( (1-S11)*(1+S22) + S12*S21) / (2*S21)

		abcd[idx][0][0] = A
		abcd[idx][0][1] = B
		abcd[idx][1][0] = C
		abcd[idx][1][1] = D

	return abcd


def sdb2sri(sdb_struct, sdeg_struct):
	# convert DB/DEG to real/imag
	num_freqs = len(sdb_struct)
	Sri = np.zeros( (num_freqs, 2, 2), dtype=complex)

	for idx in range(len(sdb_struct)):
		db_mat = sdb_struct[idx]
		S11_db = db_mat[0][0]
		S12_db = db_mat[0][1]
		S21_db = db_mat[1][0]
		S22_db = db_mat[1][1]

		deg_mat = sdeg_struct[idx]
		S11_deg = deg_mat[0][0]
		S12_deg = deg_mat[0][1]
		S21_deg = deg_mat[1][0]
		S22_deg = deg_mat[1][1]

		S11 = 10**(S11_db/20) * np.complex( math.cos(S11_deg*math.pi/180), math.sin(S11_deg*math.pi/180) )
		S12 = 10**(S12_db/20) * np.complex( math.cos(S12_deg*math.pi/180), math.sin(S12_deg*math.pi/180) )
		S21 = 10**(S21_db/20) * np.complex( math.cos(S21_deg*math.pi/180), math.sin(S21_deg*math.pi/180) )
		S22 = 10**(S22_db/20) * np.complex( math.cos(S22_deg*math.pi/180), math.sin(S22_deg*math.pi/180) )

		Sri[idx][0][0] = S11
		Sri[idx][0][1] = S12
		Sri[idx][1][0] = S21
		Sri[idx][1][1] = S22

	return Sri



def sdb2abcd(sdb_struct, sdeg_struct):

	sri = sdb2sri(sdb_struct, sdeg_struct)
	abcd = sri2abcd(sri)

	return abcd
	


def sri2sdb(sri_struct):
	# convert S params from Real/Imag to DB/Deg
	num_freqs = len(sri_struct)
	Sdb = np.zeros( (num_freqs, 2, 2))
	Sdeg = np.zeros( (num_freqs, 2, 2))

	for idx in range(len(sri_struct)):
		ri_mat = sri_struct[idx]
		S11_ri = ri_mat[0][0]
		S12_ri = ri_mat[0][1]
		S21_ri = ri_mat[1][0]
		S22_ri = ri_mat[1][1]


		S11_db = 20*np.log10( np.abs(S11_ri) )
		S12_db = 20*np.log10( np.abs(S12_ri) )
		S21_db = 20*np.log10( np.abs(S21_ri) )
		S22_db = 20*np.log10( np.abs(S22_ri) )

		S11_deg = np.arcsin( S11_ri.imag / np.abs(S11_ri) ) * 180/math.pi
		S12_deg = np.arcsin( S12_ri.imag / np.abs(S12_ri) ) * 180/math.pi
		S21_deg = np.arcsin( S21_ri.imag / np.abs(S21_ri) ) * 180/math.pi
		S22_deg = np.arcsin( S22_ri.imag / np.abs(S22_ri) ) * 180/math.pi
		
#		if ( S11_ri.real < 0 ) and (S11_ri.imag > 0):
#			S11_deg = 180 - S11_deg
#		if ( S12_ri.real < 0 ) and (S12_ri.imag > 0):
#			S12_deg = 180 - S12_deg
#		if ( S21_ri.real < 0 ) and (S21_ri.imag > 0):
#			S21_deg = 180 - S21_deg
#		if ( S22_ri.real < 0 ) and (S22_ri.imag > 0):
#			S22_deg = 180 - S22_deg
		

		Sdb[idx][0][0] = S11_db
		Sdb[idx][0][1] = S12_db
		Sdb[idx][1][0] = S21_db
		Sdb[idx][1][1] = S22_db

		Sdeg[idx][0][0] = S11_deg
		Sdeg[idx][0][1] = S12_deg
		Sdeg[idx][1][0] = S21_deg
		Sdeg[idx][1][1] = S22_deg

	return (Sdb, Sdeg)




def l2l_deembed(pad_L_s2p_filename, pad_2L_s2p_filename, structure_L_s2p_filename, structure_2L_s2p_filename, z0_probe=50):

	net_pad_L = rf.Network(pad_L_s2p_filename, z0=z0_probe) # d
	net_pad_2L = rf.Network(pad_2L_s2p_filename, z0=z0_probe) # c
	net_struct_L = rf.Network(structure_L_s2p_filename, z0=z0_probe) # e
	net_struct_2L = rf.Network(structure_2L_s2p_filename, z0=z0_probe) # f

	# ABCD matrices
	TP_2L = sdb2abcd(net_pad_2L.s_db, net_pad_2L.s_deg)
	TP_L = sdb2abcd(net_pad_L.s_db, net_pad_L.s_deg)
	TS_L = sdb2abcd(net_struct_L.s_db, net_struct_L.s_deg)
	TS_2L = sdb2abcd(net_struct_2L.s_db, net_struct_2L.s_deg)


	TL1 = []
	TL2 = []
	
	for idx, tlp_mat in enumerate(TP_L):
		tlpi_mat = la.inv(tlp_mat)
		tp_2l_mat = TP_2L[idx]
		
		TP_L_inner_pre = np.dot( tlpi_mat, np.dot( tp_2l_mat, tlpi_mat ) ) # TLPI_MAT * TS_2L_MAT * TLPI_MAT matrix multiplication
		TP_L_inner = la.inv(TP_L_inner_pre)
		TP1 = sla.sqrtm(TP_L_inner)
		TL1.append(TP1)

		TP1_inv = la.inv(TP1)
		TS_2L_entry = TS_2L[idx]
		TL2_entry = np.dot( TP1_inv, np.dot( TS_2L_entry, TP1_inv ) )
		TL2.append( TL2_entry )

	abcd_L = TL1
	abcd_2L = TL2
	Sri_L = abcd2s(TL1, z0_probe, z0_probe)
	Sri_2L = abcd2s(TL2, z0_probe, z0_probe)

	(Sdb_L, Sdeg_L) = sri2sdb(Sri_L)
	(Sdb_2L, Sdeg_2L) = sri2sdb(Sri_2L)
	freq = net_pad_L.f
	
	net_L = rf.Network( f=freq*1e-9, s=Sri_L, z0=50)
	net_2L = rf.Network( f=freq*1e-9, s=Sri_2L, z0=50)
	

	return (freq, Sri_L, Sri_2L, abcd_L, abcd_2L, Sdb_L, Sdeg_L, Sdb_2L, Sdeg_2L, net_L, net_2L)


def write_net_db_deg( net, filename):
	outfile = open(filename,'w')

	freq_vec = net.f
	s_db_list = net.s_db
	s_deg_list = net.s_deg

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

		outstr = "{0:.4g},{1:.4g},{2:.4g},{3:.4g},{4:.4g},{5:.4g},{6:.4g},{7:.4g},{8:.4g}\n".format( freq, S11_db, S11_deg, S12_db, S12_deg, S21_db, S21_deg, S22_db, S22_deg)

		outfile.write(outstr)
 


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


def write_rlgc(freq, R, L, G, C, filename):
	outfile = open(filename, 'w')
	
	for idx, f in enumerate(freq):
		r = R[idx]
		l = L[idx]
		g = G[idx]
		c = C[idx]
		
		outstr = "{0:.8g},{1:.8g},{2:.8g},{3:.8g},{4:.8g}\n".format(f, r, l, g, c)
		outfile.write(outstr)
		

