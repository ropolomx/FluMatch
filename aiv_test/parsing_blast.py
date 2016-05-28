# coding: utf-8
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
        for h in qresult.hits[0:2]:
                experiment={{"query":h.query_id}, {"strain":h.id}, {"desc":h.description}}
        
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
import csv
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:3]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:1]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:1]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[1:2]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.map
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.hsps
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.seq_len
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.filter, h.fragments
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.filter
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.fragments
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description_all
        
get_ipython().magic(u'pinfo SearchIO.Hit.mro')
get_ipython().magic(u'pinfo SearchIO.Hit.description_all')
get_ipython().magic(u'pinfo SearchIO.Hit.filter')
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.filter
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.description
        
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    print qresult.hsps
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    print qresult.id
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    print qresult.param_evalue_threshold
    
result_handle = open("flublast.xml")
blast_records = NCBIXML.parse(result_handle)
for aln in blast_records.alignments:
    print(alignment.title)
    
for rec in blast_records:
    for r in rec.alignments:
        print r.title
        
for rec in blast_records:
    for r in rec.alignments:
        print r.hit_id
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.hit_id
        
for rec in blast_records:
    for r in rec.alignments:
        print r.hit_def
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        r.accession
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        r.title
        
for rec in blast_records:
    for aln in aln.alignments[0:2]:
        print r.accession
        
for rec in blast_records:
    for aln in aln.alignments[0:2]:
        print(r.accession)
        
for rec in blast_records:
    for aln in aln.alignments[0:2]:
        print(r.title)
        
for rec in blast_records:
    for aln in aln.alignments:
        print(r.title)
        
for rec in blast_records:
    for r in rec.alignments:
        print r.title
        
for rec in blast_records:
    for r in rec.alignments:
        print r.title
        
for rec in blast_records:
    for r in rec.alignments:
        print r.length
        
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments:
        print r.length
        
result_handle = open("flublast.xml")
for rec in blast_records:
    for r in rec.alignments:
        print r.length
        
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments:
        print r.length
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.length
        
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.length
        
result_handle = open("flublast.xml")
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.length
        
result_handle = open("flublast.xml")
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.title, r.accession, r.hit_id, r.hit_def, r.hsps
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.title, r.accession, r.hit_id, r.hit_def
        
result_handle = open("flublast.xml")
blast_records = NCBIXML.parse(result_handle)
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.title, r.accession, r.hit_id, r.hit_def
        
for rec in blast_records:
    for r in rec.alignments[0:2]:
        print r.accession, r.hit_id, r.hit_def
        
result_handle = open("flublast.xml")
for rec in NCBIXML.parse(result_handle):
    for r in rec.alignments[0:2]:
        print r.accession, r.hit_id, r.hit_def
        
for rec in NCBIXML.parse(result_handle):
    for r in rec.alignments[0:2]:
        print r.hit_id, r.hit_def
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        print r.hit_id, r.hit_def
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        print r.hit_def
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        print r.hsps
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        print r.accession
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.expect
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            sub_start = hsp.sbjct_start
   sub_end = hsp.sbjct_end
   
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            sub_start = hsp.sbjct_start
            sub_end = hsp.sbjct_end
            
sub_start
print sub_start
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            sub_start[sub_end] = {hsp.subjct_start: hsp.subjct_end}
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data = hsp.sbjct_start, hsp.sbjct_end
            
subject_data
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data = {hsp.sbjct_start, hsp.sbjct_end}
            
subject_data
subject_data = {}
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data.append = {hsp.sbjct_start, hsp.sbjct_end}
            
subject_data['subject'] = {}
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data['subject'].append = {hsp.sbjct_start, hsp.sbjct_end}
            
subject_data = []
subject_data = {}
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data[hsp].append(hsp.sbjct_start, hsp.sbjct_end)
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data[subject]  = hsp.sbjct_start, hsp.sbjct_end
            
subject_data = {}
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data[subject]  = hsp.sbjct_start, hsp.sbjct_end
            
for record in NCBIXML.parse(open("flublast.xml")):
    for rec in record.alignments[0:2]:
        for hsp in rec.hsps:
            subject_start  = hsp.sbjct_start
            subject_end = hsp.sbjct_end
            subject_data[subject_end] = 1
            subject_data[subject_start] = 1
            
subject_data
for record in NCBIXML.parse(open("flublast.xml")):
    for rec in record.alignments[0:2]:
        for hsp in rec.hsps:
            subject_start  = hsp.sbjct_start
            subject_end = hsp.sbjct_end
            subject_data[subject_start] = 1
            subject_data[subject_end] = 1
subject_data
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            subject_data[subject]  = hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end
            
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:
        print h.query_id,h.id, h.description
        
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, hsp.sbjct
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, rec.query_id
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, rec.title
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, r.title
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, r.title, r.length
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.sbjct_start, hsp.sbjct_end, r.title, r.hsps
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.match[0:10]
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print r.accession, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query_id, r.accession, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print r.query_id, r.accession, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for r in rec.alignments[0:2]:
        for hsp in r.hsps:
            print r.title, hsp.sbjct_start, hsp.sbjct_end
            
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits
    print h.description, h.seq_len
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits:
    print h.description, h.seq_len
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits:
    print h.description, h.seq_len
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits:print h.description, h.seq_len
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult.hits[0:2]:print h.description, h.seq_len
    
qresult
key_func = lambda qresult: qresult.id.split('_')[1]
key_func
print key_func
search_dict = SearchIO.to_dict(qresult, key_func)
key_func = qresult.id
search_dict = SearchIO.to_dict(qresult, key_func)
key_func = lambda qresult: qresult.id[1]
search_dict = SearchIO.to_dict(qresult, key_func)
key_func = lambda qresult: qresult.id[0]
search_dict = SearchIO.to_dict(qresult, key_func)
key_func = lambda qresult: qresult.id
search_dict = SearchIO.to_dict(qresult, key_func)
search_dict['PROKKA_0009']
qresult
search_dict['PROKKA_00009']
search_dict
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult:print h.id
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult:print h.id[0:2], h.description
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult:print h.id[0:2], h.description[0:2]
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.id, h.description
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.query_id, h.id, h.description
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:for q in qresult.hits:print h.query_id, h.id, h.description
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:
        for q in qresult.hits:
            print h.query_id
            
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:
        for q in qresult.hits[0:2]:
            print h.query_id
            
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.id_all
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.id
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.fragments
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.index
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.map
    
for qresult in SearchIO.parse('flublast.xml', 'blast-xml'):
    for h in qresult[0:2]:print h.hsps.extend
    
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.id
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title, hsp.query.title
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.query.title
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.identities
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.align_length
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.positives
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title, hsp.sbjct_start, hsp.sbjct.end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print hsp.num_alignments, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title, hsp.query, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print aln.title, rec.query, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, hsp.sbjct_start, hsp.sbjct_end, hsp.strand.index
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query_id, rec.hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query_id, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, rec.reference, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end,
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hct_end, hsp.align_length, hsp.identities, rec.
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hct_end, hsp.align_length, hsp.identities, hsp.align_length
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            print rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.align_length
            
save_file = open("Pythonic_BLAST_hits.txt", "w")
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write('>%s\n' % (rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.id))entities, hsp.align_length
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write('>%s\n' % (rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.align_length))
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write('>%s\t' % (rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.align_length))
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write(rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.align_length))
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write(rec.query, aln.title, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.align_length)
            
for rec in NCBIXML.parse(open("flublast.xml")):
    for aln in rec.alignments[0:2]:
        for hsp in r.hsps:
            save_file.write(aln.title)
            
save_file.close()
