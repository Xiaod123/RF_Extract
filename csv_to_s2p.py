import re
import glob

infile_list = glob.glob("*.csv")

begin_data_re = re.compile("^BEGIN CH\d+_DATA")
end_data_re = re.compile("^END")

port_impedance_str = "50"


for infile_name in infile_list:
	infile_name_arr = infile_name.split(".")
	outfile_name = infile_name_arr[0] + ".s2p"
	infile = open(infile_name, 'r')
	outfile = open(outfile_name,'w')


	header_found = False
	header_data_gathered = False
	freq_unit_str = "HZ" # default freq units
	column_index_dict = {}
	for line in infile:

		if not header_found:
			m = begin_data_re.search(line)
			if (m):
				header_found = True
			else:
				if line[0] == "!": # preserve any comments
					outfile.write(line)

		else: # header found
			if not header_data_gathered: # this line will have header data
				# [FIX] For now assume everything has the same format
				#	Freq(Hz),S11(DB),S11(DEG),S12(DB),S12(DEG),S21(DB),S21(DEG),S22(DB),S22(DEG)
				freq_unit_str = "HZ"
				s_param_format = "DB"
				port_impedance_str = "50"
						

				header_str = "# {0:s} S {1:s} R {2:s}\n".format(freq_unit_str, s_param_format, port_impedance_str)
				outfile.write(header_str)
				header_data_gathered = True
			else: # go ahead and write all the data
				# need to remap things a bit, since S2P format requires certain order for S params
				# in S2P, it must be S11, S21, S12, S22
				nline = line.strip()
				line_arr = nline.split(",")
				if (len(line_arr) == 9):
					index_order = [0, 1, 2, 5, 6, 3, 4, 7, 8]
					line_arr_reord = [line_arr[idx] for idx in index_order]

					new_line = " ".join(line_arr_reord)
					outfile.write(new_line + "\n")
				



