### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
__version__ = "v0.14"

__citation__ = "Pucker et al., 2021: Assembly Error Finder (AEF)"

### Assembly Error Finder (AEF) ###

__usage__ = """
					python AEF.py
					--in <INPUT_BED_FILE> (--bed) | --bam <INPUT_BAM_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--fasta <ASSEMLY_FASTA_FILE>
					
					optional:
					--gff <GFF_FILE>
					--res <RESOLUTION>[100]
					--sat <SATURATION>[100]
					--name <NAME>[xxx]
					--relfreq <RELATIVE_ALIGNMENT_END_CUTOFF>[10]
					--factor <IQR_OVER_MEDIAN_CUTOFF>[5]
					--dist <EXCLUDE_DISTANCE_TO_SEQ_ENDS>[1000]
					--tolerance <DIST_TO_GENE>[0]
					"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# --- end of imports --- #


def load_bed( bed_file, resolution ):
	"""! @brief load all information from given BED file """
	
	half = 0.5*resolution
	data = {}
	with open( bed_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			
			# --- include start --- #
			block1 = int( ( int( parts[1] ) + half ) / resolution )
			try:
				data[ parts[0] ][ block1 ]+=1
			except KeyError:
				try:
					data[ parts[0] ].update( { block1: 1 } )
				except KeyError:
					data.update( { parts[0]: { block1: 1 } } )
			
			# --- include end --- #
			block2 = int( ( int( parts[2] ) + half ) / resolution )
			try:
				data[ parts[0] ][ block2 ]+=1
			except KeyError:
				try:
					data[ parts[0] ].update( { block2: 1 } )
				except KeyError:
					data.update( { parts[0]: { block2: 1 } } )
			
			line = f.readline()
	return data


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	assembly_seq_order = []
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					assembly_seq_order.append( header )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
		assembly_seq_order.append( header )
	return sequences, assembly_seq_order


def generate_aln_break_point_plot( aln_end_values, outputfile, saturation, resolution, title, contig_len, factor ):
	"""! @brief generate coverage histogram """
	
	y_values = []
	x_values = []
	for i in range( int( contig_len/ resolution) ):
		x_values.append( i*resolution )
		try:
			y_values.append( min( [ saturation, aln_end_values[ i ] ] ) )
		except KeyError:
			y_values.append( 0 )
	
	
	median = max( [ 1, np.median( y_values ) ] )
	iqr = max( [ 1, stats.iqr( y_values ) ] )
	print "median: " + str( median )
	print "IQR: " + str( iqr )
	
	plot_x_values = []
	plot_y_values = []
	for i, y in enumerate( y_values ):
		if y > median + factor * iqr:
			plot_x_values.append( x_values[ i ] )
			plot_y_values.append( y )
	
	
	# --- construct figure --- #
	fig, ax = plt.subplots()
	
	ax.set_title( title )
	#ax.plot( plot_x_values, plot_y_values, linestyle=":", marker=".", color="black", linewidth=0.01 )
	ax.scatter( plot_x_values, plot_y_values, marker="o", color="black", s=1 )
	
	ax.set_xlabel( "position on reference sequence [bp]" )
	ax.set_ylabel( "number of alignment starts/ends" )
	
	fig.savefig( outputfile, dpi=300 )
	plt.close( "all" )
	
	# --- prepare data for return --- #
	avg_ends = np.mean( y_values )
	final_data = []
	for i, val in enumerate( plot_y_values ):
		final_data.append( { 'pos': plot_x_values[ i ], 'counts': val, 'rel': val / avg_ends } )
	return final_data


def load_gene_pos_per_chr( gff_file, tolerance_space ):
	"""! @brief load gene positions from given GFF file """
	
	genes_per_chr = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[ 2 ] == "gene":
					try:
						ID = parts[-1].split('=')[1]
						if ";" in ID:
							ID = ID.split(';')[0]
					except IndexError:
						ID = parts[-1][:30]
					
					start, end = map( int, parts[3:5] )
					try:
						genes_per_chr[ parts[0] ].append( { 'ID': ID, 'start': start-tolerance_space, 'end': end-tolerance_space } )
					except KeyError:
						genes_per_chr.update( { parts[0]: [ { 'ID': ID, 'start': start-tolerance_space, 'end': end-tolerance_space } ] } )
			line = f.readline()
	return genes_per_chr


def main( arguments ):
	"""! @brief runs everything """
	
	if '--bam' in arguments:
		bam_file = arguments[ arguments.index( '--bam' ) + 1 ]
		bam_state = True
	else:
		bam_state = False
		if "--bed" in arguments:
			bed_file = arguments[ arguments.index( '--bed' ) + 1 ]
		else:
			bed_file = arguments[ arguments.index( '--in' ) + 1 ]
	output_folder = arguments[ arguments.index( '--out' ) + 1 ]
	assembly_file = arguments[ arguments.index( '--fasta' ) + 1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--res' in arguments:
		resolution = int( arguments[ arguments.index( '--res' ) + 1 ] )
	else:
		resolution = 100
	
	if '--sat' in arguments:
		saturation = int( arguments[ arguments.index( '--sat' ) + 1 ] )
	else:
		saturation = 100
	
	if '--name' in arguments:
		name = arguments[ arguments.index( '--name' ) + 1 ]
	else:
		name = "xxx"
	
	if '--factor' in arguments:
		factor = int( arguments[ arguments.index( '--factor' ) + 1 ] )
	else:
		factor=5
	
	if '--dist' in arguments:
		end_dist = int( arguments[ arguments.index( '--dist' ) + 1 ] )
	else:
		end_dist = 1000
	
	if "--tolerance" in arguments:
		tolerance_space = int( arguments[ arguments.index( '--tolerance' ) + 1 ] )
	else:
		tolerance_space = 0	#region around genes to consider when checking for effects
	
	if '--bedtools' in arguments:
		bedtools = arguments[ arguments.index( '--bedtools' ) + 1 ]
	else:
		bedtools = "bedtools"
	
	if "--relfreq" in arguments:
		relfreq_cutoff = int( arguments[ arguments.index( '--relfreq' ) + 1 ] )
	else:
		relfreq_cutoff = 10
	
	if "--gff" in arguments:
		gff_file = arguments[ arguments.index( '--gff' ) + 1 ]
		gene_info = load_gene_pos_per_chr( gff_file, tolerance_space )
	else:
		gene_info = {}
	
	if bam_state:
		bed_file = output_folder + name + ".bed"
		os.popen( bedtools + " bamtobed -i " + bam_file + " > " + bed_file )
	
	data = load_bed( bed_file, resolution )
	assembly, assembly_seq_order = load_sequences( assembly_file )

	# --- generate coverage histograms per chromosome --- #
	potential_error_documentation = output_folder + name +  "_ErrorReport.txt"
	with open( potential_error_documentation, "w" ) as doc:
		doc.write( '## This analysis was performed with Assembly Error Finder (AEF) ' + __version__ + "\n## Please cite: " + __citation__ + "\n" )
		doc.write( "## --res=" + str( resolution ) + "; This value defines the size of blocks to be analyzed.\n" )
		doc.write( "## --sat=" + str( saturation ) + "; This value defines an upper y-axis cutoff in the result visualization.\n" )
		doc.write( "## --name=" + str( name ) + "; This name was attached to all associated files.\n" )
		doc.write( "## --factor=" + str( factor ) + "; This factor defines the minimal difference in the number of read alignment ends in one region compared to the average to flag a potential error.\n" )
		doc.write( "## --dist=" + str( end_dist ) + "; This range [bp] at sequence ends is excluded from the analysis as read alignment ends are enriched in these regions.\n" )
		doc.write( "## --relfreq=" + str( relfreq_cutoff ) + "; Enrichment of alingment ends over average.\n" )
		doc.write( "## input file: " + bed_file + "\n" )
		if "--gff" in arguments:
			doc.write( "## annotation file: " + gff_file + "\n" )
		else:
			doc.write( "## annotation file: n/a\n" )
		doc.write( "#Chr\tPosition\tCountedAlignmentEnds\tEnrichmentOfAlignmentEndsOverAverage\tEffectedGenes\tComments\n" )
		for key in assembly_seq_order:
			outputfile = output_folder + name + "_" + key + ".pdf"
			try:
				content = data[ key ]
				potential_errors = generate_aln_break_point_plot( content, outputfile, saturation, resolution, key, len( assembly[ key ] ), factor )
				for error in potential_errors:
					if  end_dist <= error['pos'] <= len( assembly[ key ] )-end_dist:
						if error['rel'] > relfreq_cutoff:
							new_line = map( str, [ key, error['pos'], error['counts'], error['rel'] ] )
							try:
								relevant_genes = gene_info[ key ]
								hit_genes = []
								for gene in relevant_genes:
									if gene['start'] <= error['pos'] <= gene['end']:
										hit_genes.append( gene['ID'] )
								if len( hit_genes ) > 0:
									new_line.append( ";".join( hit_genes ) )
								else:
									new_line.append( "." )
							except KeyError:
								new_line.append( "." )
							new_line.append( "###" )
							doc.write( "\t".join( new_line ) + "\n" )
			except KeyError:
				pass


if '--bed' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--in' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--bam' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
