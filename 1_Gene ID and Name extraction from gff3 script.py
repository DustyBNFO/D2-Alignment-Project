# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 13:11:44 2023

@author: zelif
"""

import re
from datetime import datetime


today = datetime.today().date()
today_underscores = today.strftime('%Y_%m_%d')


gene_id_pattern = r"ID=gene:([^;]+)"
xref_pattern = r"xref=([^;]+)"
parent_transcript_pattern = r"parent_transcript=([^;]+)"
gene_line_no_gene_name = False
#This is to move to the next line when a gene name isn't found. This is because the parent
#transcript ID isn't found in the same line as the gene information, its found in the next line
#this is consistent throughout the entire gff3 file so this method works.
gene_names_file = "gene_names_from_annotation_",today_underscores,".txt"

with open("Annotation_File.gff3", "r") as file, open(gene_names_file, "w") as output_file:
    gene_count = 0
    no_gene_name_count = 0

    for line in file:
        if gene_line_no_gene_name:
            parent_transcript_match = re.search(parent_transcript_pattern, line)
            if parent_transcript_match:
                parent_transcript = parent_transcript_match.group(1)
                output_file.write(f"Transcript ID: {parent_transcript}\n")
                gene_line_no_gene_name = False
                no_gene_name_count += 1

        columns = line.strip().split("\t")

        if len(columns) >= 3:

            if "gene" in columns[2]:

                gene_id_match = re.search(gene_id_pattern, line)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    xref_match = re.search(xref_pattern, line)
                    if xref_match:
                        xref = xref_match.group(1)
                    else:
                        xref = "Gene Name not found"
                        gene_line_no_gene_name = True
                    output_file.write(f"{gene_id}\t{xref}\n")
                    gene_count += 1


    output_file.write(f"\nSummary:\n")
    output_file.write(f"Total genes: {gene_count}\n")
    output_file.write(f"Genes with names not found: {no_gene_name_count}\n")
    output_file.write(f"Percentage of genes with names not found: {no_gene_name_count / gene_count * 100:.2f}%\n")
        
