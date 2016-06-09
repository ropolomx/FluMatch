#!/usr/bin/env python3

###################################################################
##### FluMatch: a script to automate annotation and search of #####
############## the closest strain match of contigs ################
###################################################################

from Bio import SeqIO
from Bio.Blast import NCBIXML
from multiprocessing import cpu_count
from io import StringIO
import subprocess
import re
import os
import argparse
import csv

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--blast-db', required=True,
                        help='Path to local database to BLAST your contigs against')

    parser.add_argument('-t', '--top-hits', type=int, default=10,
                        help='Number of top BLAST hits to report for each contig [10]')

    parser.add_argument('-p', '--prokka-dir', default='./', help = 'Directory for prokka files')
    
    parser.add_argument('-r', '--report-out', default='./TopBLASTHits.txt', help='Output filename')
    
    parser.add_argument('--cores', type=int, default=cpu_count(),
                        help='Number of CPU cores to use [all]')

    parser.add_argument('contigs', help = 'FASTA formatted file containing contig sequences')
 
    return parser.parse_args()

def strain_name(fasta):
    return os.path.splitext(os.path.basename(fasta))[0]

def prokka_annotate(prokka_dir, fasta, cores):

    strain = strain_name(fasta)

    prokka = ('prokka',
              '--kingdom', 'Viruses',
              '--outdir', os.path.join(prokka_dir, strain),
              '--prefix', strain,
              '--locustag', strain,
              '--cpus', str(cores),
              fasta)

    subprocess.call(prokka)
    
    return os.path.join(prokka_dir, strain, strain + '.ffn')

def blast(query, blastdb, cores):

    blastn = ('blastn',
              '-query', query, 
              '-db', blastdb,
              '-outfmt', '5',
              '-num_threads', str(cores))

    blastn_out = subprocess.check_output(blastn)

    return NCBIXML.parse(StringIO(blastn_out))

def tabulate_hsp_xml(result, num_top_results):

    for res in result:
        for aln in res.alignments[0:num_top_results]:
            for hsp in aln.hsps:

                query = res.query

                strain = aln.hit_def

                identities = hsp.identities

                aln_length = hsp.align_length

                query_length = res.query_length

                sbjct_length = len(hsp.sbjct)

                identity_perc = 100. * identities / aln_length

                query_cov =  abs(100. * aln_length / query_length)

                query_start = hsp.query_start

                query_end = hsp.query_end

                sbjct_start = hsp.sbjct_start

                sbjct_end = hsp.sbjct_end

                e_value = hsp.expect

                yield (query, strain, query_cov, identity_perc, identities, aln_length, query_start, query_end, sbjct_start, sbjct_end, query_length, sbjct_length, e_value)

def blast_report(reportpath, result, num_top_results):

    headers = ['Query','Matching Strain','Query Coverage', 'Percent ID', 'Identities', 'Alignment Length', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'Query Length', 'Subject Length', 'e-value']
    
    with open(reportpath,'w') as f:
        
        out = csv.writer(f, delimiter='\t')
        
        out.writerow(headers)

        for hsp in tabulate_hsp_xml(result, num_top_results):
            out.writerow(hsp)

def main():

    args = arguments()
    
    annotations = prokka_annotate(args.prokka_dir, args.contigs, args.cores)
    
    blast_result = blast(annotations, args.blast_db, args.cores)

    blast_report(args.report_out, blast_result, args.top_hits)

if __name__ == '__main__':
    main()
