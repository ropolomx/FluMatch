#!/usr/bin/env python

###################################################################
##### Diagflu: a program to automate annotation and search of #####  
############## the closest strain match of AIV CDS ################
###################################################################

from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import os
import argparse
import sys

def arguments():

    parser = argparse.ArgumentParser() 
    #parser.add_argument('--blast-database', default = None, help = 'Path to database to BLAST your contigs against')

    parser.add_argument('contigs', help = 'FASTA formatted file with contig sequences')

#    parser.add_argument('--top-results', default = None, help = 'Number of top BLAST hits to report for each contig')

    return parser.parse_args()

def handle_input(filepath):
    
#    if sample_name:
#        name = sample_name 
#    elif filepath != '-':
    name = os.path.splitext(os.path.basename(filepath))[0]
#    else:
#        name = 'stdin'

    f = open (filepath)

    with f as o:
        return o.read(), name

def prokka_annotate(handle):

    prokka = ('prokka', '--kingdom', 'Viruses', '--outdir',
    'prokka_annotation_script')
    return subprocess.Popen(prokka, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            shell=True)

#def blastn_top(contig):

#def results():
#
# open("TopHitsFromBLASTXML.txt", "w") as f:
#
#    f.write("\t".join(["Query", "Strain", "SubjStart", "SubjEnd", "Identities", "QueryCoverage", "PercentID", "e-value"]))
#    f.write("\n")
#    for rec in NCBIXML.parse(open("flublast.xml")):
#        for aln in rec.alignments[0:2]:
#            for hsp in aln.hsps:
#                rec_query = rec.query
#                strain = aln.title
#		strain = re.sub(r'gnl\|(\w*|\W*)\|\d*\s','',strain)
#		sbjct_start = hsp.sbjct_start
#		sbjct_end = hsp.sbjct_end
#		identities = hsp.identities
#		query_cover = (float(hsp.sbjct_end) - float(hsp.sbjct_start - 1))/(float(hsp.query_end)-float(hsp.query_start - 1))*100
#		identity_perc = float(hsp.identities) / float(hsp.align_length) * 100
#		e_value = hsp.expect
#                f.write("\t".join([str(rec_query), str(strain), str(sbjct_start), str(sbjct_end), str(identities), str(query_cover), str(identity_perc), str(e_value)]))
#                f.write("\n")

def main():
 
    args = arguments()
    handle = handle_input(args.contigs)
    annotations = prokka_annotate(handle)

if __name__ == '__main__':
    main()
