import scipy as sp
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import matplotlib as pl
import skrf as rf # RF functions. To install this do "conda install -c scikit-rf  scikit-rf" from the command line if you have Anaconda Python installed, otherwise do "pip install scikit-rf"
import math


def distributed_rlgc_from_s2p(length_m, s2p_filename, z0_probe=50):
	# length_m:	(m)	Length of structure being measured
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.

	net = rf.Network(s2p_filename, z0=z0_probe)

	freq = net.f
	d_vec = net.a[:,1,1] # the "d" part of the abcd matrix (abcd(2,2))
	b_vec = net.a[:,0,1] # the "b" part of the abcd matrix
	gamma = [1/length_m*math.acosh(dd) for dd in d_vec]	# attenuation vector

	Zc_part = [math.sinh( gg * length_m ) for gg in gamma]
	Zc = 1/b_vec * Zc_part

	R = real( gamma * Zc)
	L = 1/2/math.pi/freq * imag(gamma * Zc)
	G = real(gamma/Zc)
	C = 1/2/pi/freq * imag(gamma/Zc)
	losstan = real(gamma/Zc) / imag(gamma/Zc)

	attenuation = 20*np.log10( abs(np.exp(-gamma*length_m)) )

	return( freq, R, L, G, C, gamma, Zc, net )


def lumped_rlgc_from_s2p(s2p_filename, z0_probe=50):
	# s2p_filename: (str)	s2p filename
	# z0_probe:	(Ohms)	Impedance of network analyzer/probes. 50 Ohm default.

	net = rf.Network(s2p_filename, z0=z0_probe)

	freq = net.f

	Z = s.z
	Y = s.y

	Zdiff = [zz[0,0] - zz[0,1] - zz[1,0] + zz[1,1] for zz in Z]
	Ycomm = [yy[0,0] + yy[0,1] + yy[1,0] + yy[1,1] for yy in Y]

	R = real(Zdiff)
	L = 1/2/pi/freq * imag(Zdiff)
	G = real(Ycomm)
	C = 1/2/pi/freq * imag(Ycomm)

	return (freq, R, L, G, C, Zdiff, Ycomm, net)


def l2l_deembed(pad_L_s2p_filename, pad_2L_s2p_filename, structure_L_s2p_filename, structure_2L_s2p_filename, z0_probe=50):

	net_pad_L = rf.Network(pad_L_s2p_filename, z0=z0_probe) # d
	net_pad_2L = rf.Network(pad_2L_s2p_filename, z0=z0_probe) # c
	net_struct_L = rf.Network(structure_L_s2p_filename, z0=z0_probe) # e
	net_struct_2L = rf.Network(structure_2L_s2p_filename, z0=z0_probe) # f

	# ABCD matrices
	TP_2L = net_pad_2L.a
	TP_L = net_pad_L.a
	TS_L = net_struct_L.a
	TS_2L = net_struct_2L.a

	TP_L_inv = la.inv(TP_L)

	TL1 = []
	TL2 = []

	for idx, tlpi_mat in enumerate(TP_L_inv):
		ts_2l_mat = TS_2L[idx]

		TP_L_inner = np.dot( tlpi_mat, np.dot( ts_2l_mat, tlpi_mat ) ) # TLPI_MAT * TS_2L_MAT * TLPI_MAT matrix multiplication
		TP1 = sla.sqrtm(TP_L_inner)
		TP1_inv = la.inv(TP1)

		TL1.append(TP1)
		TL2_entry = np.dot( TP1_inv, np.dot( TS_2L, TP1_inv ) )
		TL2.append( TL2_entry )

	# This function DOES NOT EXIST! Need to do something about this
	S_final = rf.a2s(TL2, z0_probe)

	# [FIX] Not done yet with the rest of padDeembedding.m -- continue porting next

