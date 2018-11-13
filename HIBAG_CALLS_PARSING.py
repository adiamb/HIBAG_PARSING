"""
scripts to parse HIBAG R readouts and write to a file.
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""
import sys
import argparse
import datetime
import os

parser = argparse.ArgumentParser(description='A filelist of HLA HIBAG calls, one filename per row no quote marks')
parser.add_argument('-F', required=True, help='full path to the filelist.txt')
parser.add_argument('-O', required=True, help='full path to the output file')
args=parser.parse_args()
filein = args.F
fileout = args.O

class hibag_calls(object):
	"""A hibag class that contains the Hibag calls"""
	instance_count = 0

	def __init__(self, name, filein, fileout):
		self.name = name
		self.filein = filein
		self.fileout = fileout
		self.instance_count += 1
		
	def get_attr(self):
		print(' NAME :- {} \n SOURCE FILE :- {} \n INSTANCE COUNT :- {}'.format(self.name, self.filein, self.instance_count))


	def read_filelist(self):
		filelist = []
		with open(self.filein) as f_in:
			for line in f_in:
				parse_line = line.strip()
				filelist.append(parse_line)
		return filelist

	def parse_hla(self):
		filelist = self.read_filelist()
		HLA_CALLS ={}
		n =0
		for file in filelist:
			with open(file) as hla_in:
				for line in hla_in:
					n += 1
					if n > 1:
						parse_line = line.strip().split(',')
						make_key = parse_line[1]
						make_value = parse_line[2]+','+parse_line[3]+','+str("%.3f" % float(parse_line[-2]))
						if make_key in HLA_CALLS:
							get_call = HLA_CALLS.get(make_key)
							HLA_CALLS[make_key] = get_call + ',' +make_value
						else:
							HLA_CALLS[make_key]= make_value
			n =0
		return HLA_CALLS

	#@classmethod
	def process_header(self):
		header = ['HLA_'+i.split('_')[1] for i in self.read_filelist()]
		header_out = ['PLATEID']
		for i in header:
			str_out= i+'_A1'+','+i+'_A2,'+i+'_PROB'
			header_out.append(str_out)
		header_write=','.join(header_out)+'\n'
		return header_write

	def write_out(self):
		hla_dic = self.parse_hla()
		with open(self.fileout, 'w') as hla_out:
			hla_out.write(self.process_header())
			for k, v in hla_dic.items():
				hla_out.write(k+','+v+'\n')

if __name__ == '__main__':
	hibag = hibag_calls(name='HIBAG', filein=filein, fileout=fileout)
	print(hibag.get_attr())
	hibag.write_out()






