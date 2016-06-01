#!/usr/bin/env python

###################################################################
##### Diagflu: a script to automate annotation and search of #####  
############## the closest strain match of AIV contigs ############
###################################################################

from Bio import SeqIO
from Bio.Blast import NCBIXML
# from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import re
import os
import argparse
import glob

def arguments():

    parser = argparse.ArgumentParser() 

    parser.add_argument('--blast-database', required = True, default = None, help = 'Path to database to BLAST your contigs against')

    parser.add_argument('contigs', help = 'FASTA formatted file with contig sequences')

    parser.add_argument('--top-hits', required = True, default = None, help = 'Number of top BLAST hits to report for each contig')

    parser.add_argument('--prokka-dir', help = 'Directory for prokka files')

    return parser.parse_args()

def prokka_annotate(prokka_dir, fasta):

    prokka = subprocess.call(['prokka', '--kingdom', 'Viruses',\
            '--outdir', prokka_dir, fasta]) # outdir should be strain name or sequence run ID; also set prefix and locustag to strain name

    return prokka 


def get_top_blastn_hit(num_top_results, blastdb):

    # prokka_ffn = glob.iglob('./prokka_annotation/*.ffn')
    # os.walk
    blastn = ('blastn',
              '-query', prokka_dir,
              '-db', blastdb,
              '-outfmt', '5')

    blastn_out = subprocess.check_output(blastn)
    
    results = NCBIXML.read(StringIO(blastn_out))

    with open('TopBLASTHits.txt','w') as f:

        f.write('\t'.join(['Query','Strain','SubjStart','SubjEnd',\
                'Identities','QueryCoverage','PercentID','Length',\
                 'e-value']))
        f.write('\n')
        for res in results:
            for aln in res.alignments[0:int(num_top_results)]:
                for hsp in aln.hsps:
                    res_query = res.query
                    strain=aln.title
                    strain=re.sub(r'gnl\|(\w*|\W*)\|\d*\s','',strain)
                    sbjct_start=hsp.sbjct_start
                    sbjct_end=hsp.sbjct_end
                    identities=hsp.identities
                    query_cover=(float(hsp.sbjct_end)-\
                            float(hsp.sbjct_start-1))/\
                            (float(hsp.query_end)-float(hsp.query_start\
                            -1))*100
                    identity_perc=float(hsp.identities)/\
                            float(hsp.align_length)*100
                    length=hsp.align_length
                    e_value=hsp.expect
        f.write('\t'.join([str(res_query),str(strain),str(sbjct_start),\
                str(sbjct_end),str(identities),str(query_cover),\
                str(identity_perc),str(e_value)]))
        f.write('\n')

    report = f

    return report

def main():
 
    args = arguments()
    annotations = prokka_annotate(args.prokka_dir, args.contigs)
    top_hits = get_top_blastn_hit(args.top_hits, args.blast_database)

if __name__ == '__main__':
    main()
