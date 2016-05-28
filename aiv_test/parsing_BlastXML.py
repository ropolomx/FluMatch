with open("TopHitsFromBLASTXML.txt", "w") as f:
    f.write("\t".join(["Query", "Strain", "SubjStart", "SubjEnd", "Identities", "QueryCoverage", "PercentID", "e-value"]))
    f.write("\n")
    for rec in NCBIXML.parse(open("flublast.xml")):
        for aln in rec.alignments[0:2]:
            for hsp in aln.hsps:
                rec_query = rec.query
                strain = aln.title
		strain = re.sub(r'gnl\|(\w*|\W*)\|\d*\s','',strain)
		sbjct_start = hsp.sbjct_start
		sbjct_end = hsp.sbjct_end
		identities = hsp.identities
		query_cover = (float(hsp.sbjct_end) - float(hsp.sbjct_start - 1))/(float(hsp.query_end)-float(hsp.query_start - 1))*100
		identity_perc = float(hsp.identities) / float(hsp.align_length) * 100
		e_value = hsp.expect
                f.write("\t".join([str(rec_query), str(strain), str(sbjct_start), str(sbjct_end), str(identities), str(query_cover), str(identity_perc), str(e_value)]))
                f.write("\n")
