#!/usr/bin/env python

###################################################################
##### Diagflu: a script to automate annotation and search of #####
############## the closest strain match of AIV contigs ############
###################################################################

from Bio import SeqIO
from Bio.Blast import NCBIXML
# from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing import cpu_count
from StringIO import StringIO
import subprocess
import re
import os
import argparse
import glob
import csv

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--blast-database', required=True,
                        help='Path to database to BLAST your contigs against')

    parser.add_argument('--top-hits', type=int, default=10,
                        help='Number of top BLAST hits to report for each contig [10]')

    parser.add_argument('--prokka-dir', default='./', help = 'Directory for prokka files')
    
    parser.add_argument('--report-out', default='./TopBLASTHits.txt', help='Output filename')
    
    parser.add_argument('--cores', type=int, default=cpu_count(),
                        help='Number of CPU cores to use [all]')

    parser.add_argument('contigs', help = 'FASTA formatted file with contig sequences')
 
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

def blast(query, blastdb, cores):

    blastn = ('blastn',
              '-query', query,
              '-db', blastdb,
              '-outfmt', '5',
              '-num_threads', str(cores))

    blastn_out = subprocess.check_output(blastn)

    return NCBIXML.parse(StringIO(blastn_out))


def tablulate_hsp_xml(result, num_top_results):

    for res in result:
        for aln in res.alignments[:num_top_results]:
            for hsp in aln.hsps:

                query = hsp.query

                strain = re.sub(r'\s+','_', re.sub(r'gnl\|(\w*|\W*)\|\d*\s','',
                                aln.title.strip()))

                sbjct_start = hsp.sbjct_start

                sbjct_end = hsp.sbjct_end

                identities = hsp.identities

                length = hsp.align_length

                slen = hsp.sbjct_end - (hsp.sbjct_start - 1)
                qlen = hsp.query_end - (hsp.query_start - 1)


                query_cov =  abs(100. * slen / qlen)

                identity_perc = 100. * identities / length

                e_value = hsp.expect

                yield (query, strain, sbjct_start, sbjct_end,
                       identities, query_cov, identity_perc, length, e_value)

def blast_report(reportpath, result, num_top_results):

    headers = ['Query','Strain','SubjStart','SubjEnd', 'Identities',
               'QueryCoverage','PercentID','Length','e-value']
    
    with open(reportpath,'w') as f:
        
        out = csv.writer(f, delimiter='\t')
        
        out.writerow(headers)

        for hsp in tablulate_hsp_xml(result, num_top_results):
            out.writerow(hsp)

def main():

    args = arguments()
    
    annotations = prokka_annotate(args.prokka_dir, args.contigs, args.cores)
    
    blast_result = blast(args.contigs, args.blast_database, args.cores)

    blast_report(args.report_out, blast_result, args.top_hits)

if __name__ == '__main__':
    main()
