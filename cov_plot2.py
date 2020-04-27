### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###
### reference: https://doi.org/10.3390/genes10090671 ####


__usage__ = """
					python cov_plot2.py
					--in <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--ref <FULL_PATH_TO_REFERENCE_COVERAGE_FILE>
					
					optional:
					--name <NAME>
					"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# --- end of imports --- #


def load_cov( cov_file ):
	"""! @brief load all information from coverage file """
	
	cov = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				cov.update( { header: tmp } )
				header = parts[0]
				tmp = []
			tmp.append( float( parts[-1] ) )
			line = f.readline()
		cov.update( { header: tmp } )
	return cov


def generate_plot( cov, ref_cov, out_file, resolution, saturation, name ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots( figsize=( 10, 7 ) )
	
	ymax = 10
	collected_values = {}
	
	# --- generate list for plotting --- #
	all_data = []
	for idx, key in enumerate( sorted( cov.keys() ) ):
		y = ymax-idx-1
		x = []
		blocks = [ cov[ key ] [ i : i + resolution ] for i in xrange( 0, len( cov[ key ] ), resolution ) ]
		ref_blocks = [ ref_cov[ key ] [ i : i + resolution ] for i in xrange( 0, len( ref_cov[ key ] ), resolution ) ]
		for k, block in enumerate( blocks ):
			b = np.mean( block )
			r = np.mean( ref_blocks[ k ] )
			if r > 0:
				if b > 0:
					if b != r:
						ratio = np.log2( b / r )
					else:
						ratio = 0
				else:
					ratio = -1 * saturation
			else:
				if b > 0:
					ratio = np.log2( b )
				else:
					ratio = 0
			if ratio > 0:
				x.append( min( [ ratio, saturation ] ) )
			else:
				x.append( max( [ ratio, -1 * saturation ] ) )
			all_data.append( { 'chr': key, 'pos': (k+1)*resolution, 'value': ratio } )
		collected_values.update( { key: x } )
	
	# --- plot values --- #
	for idx, key in enumerate( sorted( cov.keys() )[:5] ):
		y = ymax - ( idx*2.3 )
		x = []
		for each in collected_values[ key ]:
			x.append( y + min( [ 1, ( each / saturation ) ] ) )
		
		ax.plot( np.arange( 0, len( x ), 1 ), x, marker="o", markersize=1, linewidth=0, color="lime" )
		
		ax.text( 1, y, key )
		
		ax.plot( [ 0, len( x ) ], [ y, y ], color="black" , linewidth=0.1)
	
		ax.plot( [ 0, len( x ) ], [ y-0.66, y-0.66 ], color="grey" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y-0.33, y-0.33 ], color="grey" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y, y ], color="grey" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+0.33, y+0.33 ], color="grey" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+0.66, y+0.66 ], color="grey" , linewidth=0.1)
		
		ax.plot( [ 0, 0 ], [ y-1, y+1 ], color="black", linewidth=1, markersize=1 )
		ax.text( 0, y+1, "3", ha="right", fontsize=5 )
		ax.text( 0, y+0.66, "2", ha="right", fontsize=5 )
		ax.text( 0, y+0.33, "1", ha="right", fontsize=5 )
		ax.text( 0, y, "0", ha="right", fontsize=5 )
		ax.text( 0, y-0.33, "-1", ha="right", fontsize=5 )
		ax.text( 0, y-0.66, "-2", ha="right", fontsize=5 )
		ax.text( 0, y-1, "-3", ha="right", fontsize=5 )
		
		
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	ax.set_ylabel( "log2( relative coverage )" )
	
	ax.set_xlim( 0, 30500 )
	
	ax.set_title( name )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 10
	
	ax.xaxis.set_ticks( np.arange( 0, 31000, 1000 ) )
	labels = map( str, np.arange( 0, 31, 1 ) )
	ax.set_xticklabels( labels )
	
	plt.subplots_adjust( left=0.03, right=0.999, top=0.95, bottom=0.1 )
	
	fig.savefig( out_file, dpi=300 )
	plt.close( "all" )
	return all_data


def generate_data_output( data_output_file, all_data ):
	"""! @brief generate data outputfile for manual inspection """
	
	all_values = []
	for each in all_data:
		all_values.append( each['value'] )
	
	mean = np.mean( all_values )
	sd = np.std( all_values )
	
	all_data_with_zscore = []
	for each in all_data:
		each.update( { 'z': ( each['value'] - mean ) / sd, 'absz': abs( ( each['value'] - mean ) / sd ) } )
		all_data_with_zscore.append( each )
	
	sorted_data = sorted( all_data_with_zscore, key=itemgetter('absz') )
	with open( data_output_file, "w" ) as out:
		out.write( "Chr\tPos\tLog2(CovRatio)\tZscore\n" )
		for each in sorted_data:
			out.write( "\t".join( map( str, [ each['chr'], each['pos'], each['value'], each['z'] ] ) ) + '\n' )


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file = arguments[ arguments.index( '--in' ) + 1 ]
	output_folder = arguments[ arguments.index( '--out' ) + 1 ]
	ref_cov_file = arguments[ arguments.index( '--ref' ) + 1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	resolution = 1000
	saturation = 3.0
	
	if '--name' in arguments:
		name = arguments[ arguments.index( '--name' ) + 1 ]
	else:
		name = ""
	
	cov = load_cov( cov_file )
	ref_cov = load_cov( ref_cov_file )
	
	# --- generate per chromosome position coveage plot --- #
	out_file = output_folder + name + ".pdf"
	all_data = generate_plot( cov, ref_cov, out_file, resolution, saturation, name )
	
	# --- generate data output file --- #
	data_output_file = output_folder + name + ".txt"
	generate_data_output( data_output_file, all_data )


if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
