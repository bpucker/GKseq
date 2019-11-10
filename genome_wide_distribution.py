### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python genome_wide_distribution.py
					--in <SVIM_VCF_FILE>
					--ref <REFERENCE_SEQUENCE>
					--fig <FIGURE_NAME_WITH_PDF_EXTENSION>
					
					optional:
					--score <INT, minimal SVIM score to consider variant>[10]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, glob, re, sys
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from operator import itemgetter

# --- end of imports --- #

def load_seq_lengths( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	seq_lens = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					seq_lens.update( { header: len( seq ) } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		seq_lens.update( { header: len( seq ) } )
	return seq_lens


def plot_genome_wide_distribution( features_to_plot, chr_lengths, fig_output_file ):
	"""! @brief show genome wide distribution of genes """
	
	# --- construct plot --- #	
	fig, ax = plt.subplots(  )	#figsize=( 5, 10 )
	chr_names = sorted( chr_lengths.keys() )
	y_offset = len( chr_names )
	
	
	# --- adding chromosomes --- #
	for idx, each in enumerate( chr_names ):
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ y_offset-idx, y_offset-idx ] , color="black", linewidth=.5 )
		ax.text( chr_lengths[ each ]/1000000.0 , y_offset-idx, each, fontsize=5 )
	
	
	# --- adding gene positions --- #
	for gene in features_to_plot:
		y = y_offset-chr_names.index( gene['chr'] ) + gene['group']*0.05
		x = gene['pos'] / 1000000.0
		ax.scatter( x, y, s=1, color=gene['color'] )
	
	# --- adding legend --- #
	my_legend = [	mpatches.Patch( color="lime", label="DUP:TANDEM" ),
								mpatches.Patch( color="magenta", label="INV" ),
								mpatches.Patch( color="blue", label="INS:NOVEL" ),
								mpatches.Patch( color="purple", label="DEL" ),
								mpatches.Patch( color="black", label="DUP:INT" )
							]
	ax.legend( handles=my_legend, prop={'size': 5} )
	
	# --- improving overall layout --- #
	ax.set_xlabel( "chromosome position [Mbp]" )
	
	ax.spines["top"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_frame_on(False)
	ax.axes.get_yaxis().set_visible(False)
	
	ax.set_ylim( 0, len( chr_names )+1 )
	
	plt.subplots_adjust( left=0.0, right=0.98, top=1.0, bottom=0.15 )
	
	fig.savefig( fig_output_file, dpi=300 )


def load_features_from_vcf( vcf_file, score_cutoff ):
	"""! @brief load all features from given SVIM VCF file """
	
	features = []
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if int( parts[5] ) >= score_cutoff:
					if parts[ 4 ] == "<DUP:TANDEM>":
						features.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'color': "lime", 'group': 1 } )
					elif parts[ 4 ] == "<INV>":
						features.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'color': "magenta", 'group': 2 } )
					elif parts[ 4 ] == "<INS:NOVEL>":
						features.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'color': "blue", 'group': 3 } )
					elif parts[ 4 ] == "<DEL>":
						features.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'color': "purple", 'group': 4 } )
					elif parts[ 4 ] == "<DUP:INT>":
						features.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'color': "black", 'group': 5 } )
					else:
						print line
			line = f.readline()
	return features


def main( arguments ):
	"""! @brief run everything """
	
	vcf_file = arguments[ arguments.index( '--in' )+1 ]
	fig_output_file = arguments[ arguments.index( '--fig' )+1 ]
	genome_seq_file = arguments[ arguments.index( '--ref' )+1 ]
	
	if '--score' in arguments:
		score_cutoff = int( arguments[ arguments.index( '--score' )+1 ] )
	else:
		score_cutoff = 10
	
	features_to_plot = load_features_from_vcf( vcf_file, score_cutoff )
	
	chr_lengths = load_seq_lengths( genome_seq_file )
	plot_genome_wide_distribution( features_to_plot, chr_lengths, fig_output_file )


if '--in' in sys.argv and '--ref' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
