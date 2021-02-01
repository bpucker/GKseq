### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python group_to_regions.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					
					optional:
					--win <WINDOW_SIZE>[30000]
					"""

import os, sys
import numpy as np

# --- end of imports --- #

def merge_regions( data, window ):
	"""! @brief merge regions """
	
	updated_data = {}
	for chromosome in data.keys():
		updated_info = []
		info = data[ chromosome ]
		i = 0
		while i < len( info ):
			current_pos = info[ i ]['pos']
			positions = [ info[ i ]['pos'] ]
			counts = [ info[ i ]['count'] ]
			rels = [ info[ i ]['rel'] ]
			genes = [ info[ i ]['gene'] ]
			if i < len( info ) - 1:
				while info[ i+1 ]['pos'] - current_pos < window:
					positions.append( info[ i+1 ]['pos'] )
					counts.append( info[ i+1 ]['count'] )
					rels.append( info[ i+1 ]['rel'] )
					genes.append( info[ i+1 ]['gene'] )
					i += 1
					if i == len( info )-1:
						break
					current_pos = 0 + info[ i ]['pos']
			updated_info.append( { 'pos': int( np.mean( positions ) ), 'count': max( counts ), 'rel': max( rels ), 'genes': ";".join( sorted( list( set( genes ) ) ) ).replace( ";.", "" ).replace( ".;", "" ) } )
			i += 1
		updated_data.update( { chromosome: updated_info } )
	return updated_data


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--win' in arguments:
		window = int( arguments[ arguments.index('--win')+1 ] )
	else:
		window = 30000

	data = {}
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			prev_chr = False
			while line:
				if line[0] == "#":
					out.write( line )
				else:
					parts = line.strip().split('\t')
					try:
						data[ parts[0] ].append( { 'pos': int( parts[1] ), 'count': int( parts[2] ), 'rel': float( parts[3] ), 'gene': parts[4] } )
					except:
						data.update( { parts[0]: [ { 'pos': int( parts[1] ), 'count': int( parts[2] ), 'rel': float( parts[3] ), 'gene': parts[4] } ] } )
				line = f.readline()
		sorted_merged_data = merge_regions( data, window )
		for key in sorted( sorted_merged_data.keys() ):
			content = sorted_merged_data[ key ]
			for each in content:
				out.write( "\t".join( map( str, [ key, each['pos'], each['count'], each['rel'], each['genes'] ] ) ) + "\n" )


if "--in" in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
