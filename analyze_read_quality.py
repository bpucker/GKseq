### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python analyze_read_quality.py
					--in <FASTQ_INPUT_FILE>
					--info <INFO_FILE>
					--out <OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, gzip, sys
import matplotlib.pyplot as plt
import numpy as np

# --- end of imports --- #


def load_infos( info_file ):
	"""! @brief load infos from given file """
	
	data = []
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			data.append( { 'line': parts[0], 'read': parts[1], 'start': int( parts[2] ), 'end': int( parts[3] ), 'comment': parts[4] } )
			line = f.readline()
	return data


def select_reads( IDs, read_selection_file, input_file ):
	"""! @brief load reads of interest from FASTQ file """
	
	with open( read_selection_file, "w" ) as out:
		with gzip.open( input_file, "r" ) as f:
			line = f.readline()
			while line:
				f.readline()	#seq
				f.readline()	#useless
				qual = f.readline()
				try:
					IDs[ line.split(' ')[0][1:] ]
					out.write( line + qual )
				except KeyError:
					pass
				line = f.readline()


def generate_figure( fig_file, quals, start, end, offset, window_size, step ):
	"""! @brief generate figure per read """
	
	fig, ax = plt.subplots()
	
	x_values = []
	y_values = []
	s = 0
	e = 0 + window_size
	counter = 0
	while s < len( quals ):
		x_values.append( counter )
		y_values.append( ( sum( quals[ s:e ] ) / float( window_size ) ) - offset )
		s += window_size
		e += window_size
		counter += 1
	
	ax.plot( [ start / float( window_size ), end / float( window_size ) ], [ 1, 1 ], color="grey", alpha=0.5 )
	ax.scatter( x_values, y_values, marker=".", s=1, color="black" )
	
	ax.set_ylim( 0, 30 )
	ax.set_xlim( 0, len( x_values ) )
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks( np.arange( 0, len( x_values ), 100) )
	
	final_labels = []
	for value in np.arange( 0, len( x_values ), 100 ):
		final_labels.append( str( ( value*window_size ) / 1000 ) )
	
	ax.set_xticklabels( final_labels )
	ax.set_xlabel( "position in read [kbp]" )
	ax.set_ylabel( "Phred score" )
	
	
	fig.savefig( fig_file, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief run everything """

	input_file = arguments[ arguments.index('--in')+1 ]
	info_file = arguments[ arguments.index('--info')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]

	offset = 33
	window_size = 200
	step = 100

	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	read_selection_file = output_folder + "selection.txt"


	infos = load_infos( info_file )
	IDs = {}
	for each in infos:
		IDs.update( { each['read']: None } )


	if not os.path.isfile( read_selection_file ):
		select_reads( IDs, read_selection_file, input_file )

	# --- load all quality values into dictionary --- #
	read_quals = {}
	with open( read_selection_file, "r" ) as f:
		line = f.readline()
		while line:
			ID = line.split(' ')[0][1:]
			quals = map( ord, f.readline().strip() )
			read_quals.update( { ID: quals } )
			line = f.readline()

	# --- generate one plot per read --- #
	for read in infos:
		fig_file =  output_folder + read['read'] + ".pdf"
		generate_figure( fig_file, read_quals[ read['read'] ], read['start'], read['end'], offset, window_size, step )


if '--in' in sys.argv and '--info' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
