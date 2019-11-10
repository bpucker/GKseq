### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python filter_SVIM_VCF.py
					--in <INPUT_VCF>
					--out <OUTPUT_VCF>
					--minsize <MINIMAL_SV_SIZE>[1000]
					"""

import sys, os, glob

# --- end of imports --- #

def check_SVs( vcf, outputfile, size_cutoff ):
	"""! @brief check detected SVs """
	
	with open( outputfile, "w" ) as out:
		with open( vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[6] == "PASS":
						if "SVLEN" in line:
							length = int( parts[7].split('SVLEN=')[1].split(';')[0] )
							if abs( length ) > size_cutoff:
								out.write( line )
						else:
							end = int( parts[7].split('END=')[1].split(';')[0] )
							if abs( int( parts[1] ) - end ) > size_cutoff:
								out.write( line )
				line = f.readline()


def main( arguments ):
	"""! @brief run everything """

	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--minsize' in arguments:
		size_cutoff = int( arguments[ arguments.index('--minsize')+1 ] )
	else:
		size_cutoff=1000
	
	check_SVs( input_file, output_file, size_cutoff )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
