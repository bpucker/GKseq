### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python select_long_reads.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					optional:
					--cut <LENGTH_CUTOFF_IN_BP>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, gzip
import numpy as np

# --- end of imports --- #


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	if '--cut' in arguments:
		cut = int( arguments[ arguments.index( '--cut' )+1 ] )+1	#add one to account for new line character
	else:
		cut = 3001	#add one to account for new line character
	
	# --- open files for reading and writing --- #
	if output_file.split('.')[-1] not in [ "gz", "GZ", "gzip", "GZIP" ]:
		sys.exit( "ERROR: only gziped files are supported!" )
	else:
		if input_file.split('.')[-1] not in [ "gz", "GZ", "gzip", "GZIP" ]:
			sys.exit( "ERROR: only gziped files are supported!" )
		else:
			total_seq_counter = []
			with gzip.open( output_file, "wb" ) as out:
				with gzip.open( input_file, "rb" ) as f:
					line = f.readline()
					while line:
						seq = f.readline()
						x = f.readline()
						qual = f.readline()
						if len( seq ) > cut:
							out.write( line )
							out.write( seq )
							out.write( x )
							out.write( qual )
							total_seq_counter.append( len( seq )-2 )
						line = f.readline()
	
	# --- reports --- #
	print "total length of remaining reads: " + str( sum( total_seq_counter ) )
	print "total number of reads: " + str( len( total_seq_counter ) )
	print "average read length (mean): " + str( np.mean( total_seq_counter ) )
	print "average read length (median): " + str( np.median( total_seq_counter ) )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
