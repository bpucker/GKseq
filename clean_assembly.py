### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###

__usage__ = """
					python clean_assembly.py
					--out <OUTPUT_FOLDER>
					--assembly <ASSEMBLY_FILE>
					--cov <COVERAGE_FILE>
					--organel <ORGANEL_SEQ_FILE>
					--white <WHITE_LIST_SEQ_FILE>
					--black <BLACK_LIST_SEQ_FILIE>
					
					optional:
					--minlen <INT, MINIMAL_CONTIG_LENGTH>
					--mincov <INT, MINIMAL_READ_MAPPING_COVERAGE_FOR_FILTERING>
					--maxcov <INT, MAXIMAL_READ_MAPPING_COVERAGE_FOR_FILTERING>
					--wordsize <INT, WORD_SIZE_FOR_BLAST>
					--cpus <INT, NUMBER_OF_THREADS_FOR_BLAST_SEARCH>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import os, sys, re

# --- end of imports --- #

def construct_figure( fig_file, cov_values ):
	"""! @brief construct coverage figure per contig """
	
	fig, ax = plt.subplots()
	ax.hist( cov_values, bins=50000 )
	ax.set_xlim( 0, 5000 )
	ax.set_xlabel( "average coverage" )
	ax.set_ylabel( "number of contigs" )
	fig.savefig( fig_file, dpi=300 )
	
	plt.close( "all" )


def coverage_investigation( output_file, sequence_status, cov_file, fig_file ):
	"""! @brief load all coverage values from given coverage file and calculate average """
	
	infos = {}
	with open( output_file, "w" ) as out:
		out.write( "contigID\taverage_coverage\tlength\tBLAST_status\n" )
		average_cov = []
		with open( cov_file, "r" ) as f:
			line = f.readline()
			cov_values = []
			header = line.split('\t')[0]
			while line:
				parts = line.strip().split('\t')
				if parts[0] != header:
					mean = sum( cov_values ) / len( cov_values )
					average_cov.append( mean )
					out.write( header + '\t' + str( mean ) + '\t' + str( len( cov_values ) ) + '\t' + sequence_status[ header ] +  '\n' )
					infos.update( { header: { 'cov': str( mean ), 'len': str( len( cov_values ) ), 'blast': sequence_status[ header ] } } )
					header = parts[0]
					cov_values = []
				cov_values.append( int( float( parts[-1] ) ) )
				line = f.readline()
			mean = sum( cov_values ) / len( cov_values )
			average_cov.append( mean )
			out.write( header + '\t' + str( mean ) + '\t' + str( len( cov_values ) ) + '\t' + sequence_status[ header ] + '\n' )
			infos.update( { header: { 'cov': str( mean ), 'len': str( len( cov_values ) ), 'blast': sequence_status[ header ] } } )
		
		construct_figure( fig_file, average_cov )
	return infos


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'len': int( parts[3] ) } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'len': int( parts[3] ) } } )
			line = f.readline()
	return best_hits


def set_alignment_based_sequence_status( assembly, organel_hits, white_hits, black_hits ):
	"""! @brief define BLAST based status for each sequence in assembly """
	
	seq_status = {}
	for key in assembly.keys():
		try:
			org = organel_hits[ key ]
			if org['len'] > ( len( assembly[ key ] )*0.4 ):
				seq_status.update( { key: "organel" } )
			else:
				try:
					white = white_hits[ key ]
					try:
						black = black_hits[ key ]
						if white['score'] >= black['score']:
							if white['score'] > org['score']:
								seq_status.update( { key: "white" } )
							else:
								seq_status.update( { key: "organel" } )
						else:
							seq_status.update( { key: "black" } )
					except KeyError:
						if white['score'] > org['score']:
							seq_status.update( { key: "white" } )
						else:
							seq_status.update( { key: "organel" } )
				except KeyError:
					try:
						black = black_hits[ key ]
						seq_status.update( { key: "black" } )
					except KeyError:
						seq_status.update( { key: "unknown" } )
		except KeyError:
			try:
				white = white_hits[ key ]
				try:
					black = black_hits[ key ]
					if white['score'] >= black['score']:
						seq_status.update( { key: "white" } )
					else:
						seq_status.update( { key: "black" } )
				except KeyError:
					seq_status.update( { key: "white" } )
			except KeyError:
				try:
					black = black_hits[ key ]
					seq_status.update( { key: "black" } )
				except KeyError:
					seq_status.update( { key: "unknown" } )
	return seq_status


def main( arguments ):
	"""! @brief run everything """
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	assembly_file = arguments[ arguments.index('--assembly')+1 ]
	cov_file = arguments[ arguments.index('--cov')+1 ]
	organel_file = arguments[ arguments.index('--organel')+1 ]
	positive_blast_file = arguments[ arguments.index('--white')+1 ]
	negative_blast_file = arguments[ arguments.index('--black')+1 ]
	
	if '--minlen' in arguments:
		min_len = int( arguments[ arguments.index('--minlen')+1 ] )
	else:
		min_len = 100000
	
	if '--mincov' in arguments:
		min_cov = int( arguments[ arguments.index('--mincov')+1 ] )
	else:
		min_cov = 10
		
	if '--maxcov' in arguments:
		max_cov = int( arguments[ arguments.index('--maxcov')+1 ] )
	else:
		max_cov = 200
	
	if '--wordsize' in arguments:
		word_size = int( arguments[ arguments.index('--wordsize')+1 ] )
	else:
		word_size = 30
	
	if '--cpus' in arguments:
		cpus = int( arguments[ arguments.index('--cpus')+1 ] )
	else:
		cpus = 20
	
	
	blast_status_file = output_folder + "blast_status.txt"
	fig_file = output_folder + "mapping_coverage.pdf"
	cov_info_file = output_folder + "mapping_cov_per_contig.txt"
	clean_assembly_file = output_folder + "clean_assembly_file.fasta"
	final_doc_file = output_folder + "FINAL_DOCS.txt"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- run sequence comparison via BLAST --- #
	organel_db = output_folder + "organel_db"
	organel_result_file = output_folder + "organel_result_file.txt"
	if not os.path.isfile( organel_result_file ):
		os.popen( "makeblastdb -in " + organel_file + " -out " + organel_db + " -dbtype nucl" )
		os.popen( "blastn -query " + assembly_file + " -db " + organel_db + " -out " + organel_result_file + " -outfmt 6 -evalue 0.0000000001 -num_threads "+str( cpus )+" -word_size "+str( word_size ) )
	organel_hits = load_best_blast_hit( organel_result_file )
	
	white_db = output_folder + "white_db"
	white_result_file = output_folder + "white_result_file.txt"
	if not os.path.isfile( white_result_file ):
		os.popen( "makeblastdb -in " + positive_blast_file + " -out " + white_db + " -dbtype nucl" )
		os.popen( "blastn -query " + assembly_file + " -db " + white_db + " -out " + white_result_file + " -outfmt 6 -evalue 0.0000000001 -num_threads "+str( cpus )+" -word_size "+str( word_size ) )
	white_hits = load_best_blast_hit( white_result_file )
	
	black_db = output_folder + "black_db"
	black_result_file = output_folder + "black_result_file.txt"
	if not os.path.isfile( black_result_file ):
		os.popen( "makeblastdb -in " + negative_blast_file + " -out " + black_db + " -dbtype nucl" )
		os.popen( "blastn -query " + assembly_file + " -db " + black_db + " -out " + black_result_file + " -outfmt 6 -evalue 0.0000000001 -num_threads "+str( cpus )+" -word_size "+str( word_size ) )
	black_hits = load_best_blast_hit( black_result_file )
	
	assembly = load_sequences( assembly_file )
	
	sequence_status = set_alignment_based_sequence_status( assembly, organel_hits, white_hits, black_hits )
	with open( blast_status_file, "w" ) as out:
		for key in sequence_status.keys():
			out.write( key + '\t' + sequence_status[ key ] + '\n' )
	
	# --- run coverage analysis --- #
	infos = coverage_investigation( cov_info_file, sequence_status, cov_file, fig_file )
	
	with open( final_doc_file, "w" ) as doc:
		with open( clean_assembly_file, "w" ) as out:
			with open( cov_info_file, "r" ) as f:
				f.readline()	#remove header
				line = f.readline()
				while line:
					parts = line.strip().split('\t')
					if min_cov <= int( parts[1] ) <= max_cov:
						if int( parts[2] ) > min_len:
							if sequence_status[ parts[0] ] == "white" or sequence_status[ parts[0] ] == "unknown":
								out.write( '>' + re.findall( "tig\d+", parts[0] )[0] + "\n" + assembly[ parts[0] ] + '\n' )
								doc.write( "\t".join( map( str, [ re.findall( "tig\d+", parts[0] )[0], infos[ parts[0] ]['cov'], infos[ parts[0] ]['len'], infos[ parts[0] ]['blast'] ] ) ) + '\n' )
							if sequence_status[ parts[0] ] == "unknown":
								print "WARNING: unknown sequence incluced - " + parts[0]
					line = f.readline()


if '--out' in sys.argv and '--assembly' in sys.argv and '--cov' in sys.argv and '--organel' in sys.argv and '--white' in sys.argv and '--black' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
