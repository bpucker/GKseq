### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v1.2 ###

## based on: Pucker et al., 2019; https://doi.org/10.1371/journal.pone.0216233 ###

__usage__ = """
						python dot_plot_heatmap.py\n
						--in1 <FULL_PATH_TO_FASTA_FILE1>
						--in2 <FULL_PATH_TO_FASTA_FILE2>
						--out <FULL_PATH_TO_OUTPUT_DIRECTORY>[.]\n
						
						
						OPTIONAL:
						--show	dot plot heatmap will be displayed as interactive figure
						--cite	will not run the script, but display the reference to cite for it
						--block <INT, size of blocks represented by one dot [bp]>[1000]
						--namex <STRING, name of first sequence>
						--namey <STRING, name of second sequence>
						--title <STRING, title of plot>[""]
						--name <PREFIX_OF_RESULT_FILES>[""]
						bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
						"""

__reference__ = """Pucker et al., 2019: https://doi.org/10.1371/journal.pone.0216233"""



import sys, os, re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
from operator import itemgetter

# --- end of imports --- #


def chunk_generator( seq, chunk_size):
    """! @brief generates successive n-sized chunks from input sequence """
    
    for i in xrange( 0, len( seq ), chunk_size):
        yield seq[ i:i + chunk_size]


def construct_seq_block_file( input_file, output_file, block_size ):
	"""! @brief construct file with chunked sequence """
	
	counter = 0
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			f.readline() #remove header
			seq = []
			line = f.readline()
			while line:
				if line[0] == ">":
					chunks = chunk_generator( "".join( seq ), block_size)
					for each in chunks:
						out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
						counter += 1
					seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			chunks = chunk_generator( "".join( seq ), block_size)
			for each in chunks:
				out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
				counter += 1


def load_blast_results( blast_result_file ):
	"""! @brief load blast results for plot """
	
	blast_results = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hit = blast_results[ parts[0] ]
				if float( parts[-1] ) > hit['score']:
					del blast_results[ parts[0] ]
					blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			except KeyError:
				blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			line = f.readline()
	return blast_results


def calculate_color( value ):
	"""! @brief calculates color for given value """
		
	r = 255-int( 255*value**4 )
	g = 255-int( 255*value**4 )
	b = 255
		
	color = '#%02x%02x%02x' % (r, g, b)
	
	return color


def construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, namex, namey, title ):
	"""! @brief construct dot plot with score weighted points """
	
	# --- construction of plot --- #
	fig, ax = plt.subplots( figsize=(10,10) )
	
	maxx = 0
	maxy = 0
	for key in self_blast_results.keys():
		try:
			other_point = other_blast_results[ key ]
			self_point = self_blast_results[ key ]
			score = other_point[ 'score' ] / self_point[ 'score' ]
			point_color = calculate_color( score )
			x = ( self_point['pos']  ) / 1000000000.0
			if x > maxx:
				maxx = 0 + x
			y = ( other_point['pos']  ) / 1000000000.0
			if y > maxy:
				maxy = 0 + y
			ax.scatter( x, y, color=point_color,s=1 )
		except KeyError:
			pass
	
	ax.set_xlabel( namex + " [Gbp]" )
	ax.set_ylabel( namey + " [Gbp]" )
	
	ax.set_xlim( 0, maxx )
	ax.set_ylim( 0, maxy )
	
	ax.ticklabel_format( axis='y', style="sci", scilimits=(-2,2) )
	ax.ticklabel_format( axis='x', style="sci", scilimits=(-2,2) )
	
	patches = []
	for i in range( 11 ):
		patches.append( mpatches.Patch(color=calculate_color( i/10.0 ), label=str(i/10.0)  ) )
	ax.legend( handles=patches, loc='upper left' )	#'lower right'
	
	ax.set_title( title )
	
	if show_status:
		plt.show()
	
	fig.savefig( output_figure, dpi=1200 )


def load_all_seqs_from_multiple_fasta_file( filename ):
	"""! @brief load all sequences from multiple fasta file """
	
	data = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				data.update( { header: "".join( seq ) } )
				header = line.strip()[1:].split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		data.update( { header: "".join( seq ) } )
	return data


def generate_concatenated_seq( infile, outfile ):
	"""! @brief generate new sequence file with concatenated sequence """
	
	seq = []
	with open( infile, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '>':
				seq.append( line.strip() )
			line = f.readline()
	
	with open( outfile, "w" ) as out:
		out.write( '>seq\n' + "".join( seq ) + '\n' )


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	in_seq_ref_file1 = arguments[  arguments.index( '--in1' )+1 ]
	in_seq_ref_file2 = arguments[  arguments.index( '--in2' )+1 ]
	
	if '--show' in arguments:
		show_status = True
	else:
		show_status = False
		
	if '--block' in arguments:
		block_size = int( arguments[  arguments.index( '--block' )+1 ] )
	else:
		block_size = 1000
	
	if '--namex' in arguments:
		namex = arguments[  arguments.index( '--namex' )+1 ]
	else:
		namex = "x"
	
	if '--namey' in arguments:
		namey = arguments[  arguments.index( '--namey' )+1 ]
	else:
		namey = "y"
	
	if '--name' in arguments:
		name = arguments[  arguments.index( '--name' )+1 ]
	else:
		name = ""

	if '--title' in arguments:
		title = arguments[  arguments.index( '--title' )+1 ]
	else:
		title = ""
	
	prefix = arguments[  arguments.index( '--out' )+1 ]
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	output_figure = prefix + name +  "_dot_plot_heatmap.pdf"
	seq_ref_file1 = prefix + name + "_seq1.fasta"
	seq_ref_file2 = prefix + name +  "_seq2.fasta"
	
	generate_concatenated_seq( in_seq_ref_file1, seq_ref_file1 )
	generate_concatenated_seq( in_seq_ref_file2, seq_ref_file2 )
	
	
	active = True
		
	# --- set output files --- #
	self_blast_result_file = prefix + "self_blast_result_file.txt"
	other_blast_result_file = prefix + "other_blast_result_file.txt"
	
	if active:
		seq_block_file1 = prefix + "seq_blocks.fasta"
		construct_seq_block_file( seq_ref_file1, seq_block_file1, block_size )
		
		# --- run blast vs. self and vs. other --- #
		self_blast_db = prefix + "self_blast_db"
		other_blast_db = prefix + "other_blast_db"
		if active:
			os.popen( "makeblastdb -in " + seq_ref_file1 + " -out " +  self_blast_db + " -dbtype nucl" )
			os.popen( "makeblastdb -in " + seq_ref_file2 + " -out " +  other_blast_db + " -dbtype nucl" )
		
		if active:
			os.popen( "blastn -query " + seq_block_file1 + " -db " + self_blast_db + " -out " +  self_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
			os.popen( "blastn -query " + seq_block_file1 + " -db " + other_blast_db + " -out " +  other_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
		
	# --- load blast results --- #
	self_blast_results = load_blast_results( self_blast_result_file )
	other_blast_results = load_blast_results( other_blast_result_file )
	
	
	print "Generating dot plot heatmap ..."
	construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, namex, namey, title )


if '--cite' in sys.argv:
	sys.exit( __reference__ )
	
if '--help' in sys.argv or '-h' in sys.argv:
	sys.exit( __usage__ )

elif '--out' in sys.argv and '--in1' in sys.argv and '--in2' in sys.argv:
	main( sys.argv )

else:
	sys.exit( __usage__ )
