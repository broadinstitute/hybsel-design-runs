In the the viral accession list (updated in commit 1699608), there are
very few Influenza genomes (only RefSeq genomes). So only a few
Influenza genomes (namely, three genomes - 8 segments each - one for
A/B/C) are downloaded from that list. (In the previous version of the
accession list, there similarly were only RefSeq genomes; however, it
had included other Influenza genomes in the output because I had
merged the downloaded sequences from the accession list with sequences
that I had previously downloaded, and I had previously downloaded a
selection of Influenza sequences.) NCBI has confirmed that it does not
include Influenza sequences as part of Genome Neighbors; thus, they
would be left out of the viral accession list. See
<https://support.ncbi.nlm.nih.gov/ics/support/KBAnswer.asp?questionID=1958&folderID=0&subscribe=1>,
which says:
    Exception to note: Genome neighbors are not calculated for the
    influenza sequences and you will not be able to obtain the sequences
    for this taxon as described above. The GenBank influenza sequences are
    distributed separately through the INFLUENZA FTP site maintained
    through the Influenza Virus Resource page. 

To download Influenza sequences, I used the Influenza Virus Resource
at <https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database>.
For each of Type A, B, and C, I added a query that had host 'Human',
left all the other selections at 'any', and checked the box next to
'Full-length only' sequences. I selected 'Nucleotide' sequence type
at the top. This left the following number of sequences:
    Type A: 194737
    Type B: 57256
    Type C: 743
     TOTAL: 252736
Here's a permanent link for the search:
<https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?cmd=&go=database&accfile=&acc=&sequence=N&search=&searchin=strain&type=c&host=Human&country=any&segment=any&subtype_h=any&subtype_n=any&flen=&tlen=&sonly=on&fyear=&fmonth=&fday=&tyear=&tmonth=&tday=&fpyear=&fpmonth=&fpday=&tpyear=&tpmonth=&tpday=&showfilters=&swine=include&niaid=include&lab=exclude&vac_strain=include&lineage=include&subtype_mix=include&query_1_line=on&query_1_line_num=1&query_1_country=any&query_1_host=Human&query_1_lineage=include&query_1_searchin=strain&query_1_segment=any&query_1_sequence=N&query_1_subtype_h=any&query_1_subtype_mix=include&query_1_subtype_n=any&query_1_type=a&query_1_vac_strain=include&query_1_count=194737&query_1_count_genome_sets=0&query_1_query_key=1&query_1_sonly=on&query_1_swine=include&query_1_niaid=include&query_1_lab=exclude&query_1_vac_strain=include&query_1_lineage=include&query_1_subtype_mix=include&query_2_line=on&query_2_line_num=2&query_2_country=any&query_2_host=Human&query_2_lineage=include&query_2_searchin=strain&query_2_segment=any&query_2_sequence=N&query_2_subtype_h=any&query_2_subtype_mix=include&query_2_subtype_n=any&query_2_type=b&query_2_vac_strain=include&query_2_count=57256&query_2_count_genome_sets=0&query_2_query_key=1&query_2_sonly=on&query_2_swine=include&query_2_niaid=include&query_2_lab=exclude&query_2_vac_strain=include&query_2_lineage=include&query_2_subtype_mix=include&query_3_line=on&query_3_line_num=3&query_3_country=any&query_3_host=Human&query_3_lineage=include&query_3_searchin=strain&query_3_segment=any&query_3_sequence=N&query_3_subtype_h=any&query_3_subtype_mix=include&query_3_subtype_n=any&query_3_type=c&query_3_vac_strain=include&query_3_count=743&query_3_count_genome_sets=0&query_3_query_key=1&query_3_sonly=on&query_3_swine=include&query_3_niaid=include&query_3_lab=exclude&query_3_vac_strain=include&query_3_lineage=include&query_3_subtype_mix=include&download-select=tab&defline_saved=%3E%7Baccession%7D%20%7Bstrain%7D%20%7Byear%7D/%7Bmonth%7D/%7Bday%7D%20%7Bsegname%7D&defline=%3E%7Baccession%7D%20%7Bstrain%7D%20%7Byear%7D/%7Bmonth%7D/%7Bday%7D%20%7Bsegname%7D&>
Note that many of these sequences are identical to each other (I did not
click 'Collapse identical sequences'); there are 155839 unique sequences
among the 252736.
I downloaded the results in table format (influenza.virus-resource-download.txt).

I converted this to a format like the viral accession list using
the script convert_to_acc_list.sh - this output influenza.acc-list.txt.

Then, using the Influenza accession list created above, I downloaded
fastas like I do as usual with the viral accession list. In particular,
I ran download_influenza_seqs.sh to call the download script. This
places the downloaded fastas in ../all/.
