from Bio import SeqIO

# Replace with your actual GenBank file path
genbank_file = "CP046355.gb"
output_file = "CP046355_gene_sizes_and_cds_products_with_notes.txt"

# Parse the file and write output to a file 
with open(genbank_file, "r") as input_handle, open(output_file, "w") as output_handle:
    output_handle.write("Gene/LocusTag\tSize\tCDS Product\tNote\n")
    for record in SeqIO.parse(input_handle, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                # Get the gene name or locus tag
                gene_name = feature.qualifiers.get("gene", ["[unknown]"])[0]
                if gene_name == "[unknown]":
                    locus_tag = feature.qualifiers.get("locus_tag", ["[no locus tag]"])[0]
                    gene_name = locus_tag if locus_tag != "[no locus tag]" else "[unknown]"

                start = feature.location.start.position
                end = feature.location.end.position
                gene_length = end - start

                # Default CDS product and note
                cds_product = "[no CDS product]"
                note = "[no note]"

                # Check for corresponding CDS feature for product extraction
                for sub_feature in record.features:
                    if sub_feature.type == "CDS":
                        # Match based on gene name or locus tag
                        sub_gene_name = sub_feature.qualifiers.get("gene", ["[unknown]"])[0]
                        sub_locus_tag = sub_feature.qualifiers.get("locus_tag", ["[no locus tag]"])[0]

                        if sub_gene_name == gene_name or sub_locus_tag == gene_name:
                            cds_product = sub_feature.qualifiers.get("product", ["[no CDS product]"])[0]
                            note = sub_feature.qualifiers.get("note", ["[no note]"])[0]
                            break

                output_handle.write(f"{gene_name}\t{gene_length}\t{cds_product}\t{note}\n")
