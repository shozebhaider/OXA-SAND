from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import pandas as pd
import os
import shutil


def upload_and_save_query():
    """Handles the upload of a query FASTA file and copies it as query.fasta."""
    file_name = input("Please enter the path to your query FASTA file: ")

    if not os.path.exists(file_name):
        raise FileNotFoundError(f"File not found: {file_name}")

    import shutil
    shutil.copy(file_name, "query.fasta")
    print("File uploaded and saved as query.fasta")

    return "query.fasta"

def align_sequences():
    """Concatenates reference and query FASTA files, aligns them to the HMM profile, and converts to Clustal format."""
    # Concatenate reference and query FASTA files
    with open("concatenated_output.fasta", "w") as output_handle:
        with open("./files/OXA-48.fasta", "r") as file1:
            for record in SeqIO.parse(file1, "fasta"):
                SeqIO.write(record, output_handle, "fasta")
        with open("query.fasta", "r") as file2:
            for record in SeqIO.parse(file2, "fasta"):
                SeqIO.write(record, output_handle, "fasta")

    # Align concatenated FASTA file to HMM profile
    os.system("hmmalign ./files/DBL_profile.hmm concatenated_output.fasta > query_ref_aligned.sto")

    # Convert alignment from Stockholm to Clustal format
    with open("query_ref_aligned.sto", "r") as input_handle, open("query_ref_aligned.aln", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "stockholm")
        SeqIO.write(sequences, output_handle, "clustal")

    print("Query sequence aligned to the HMM profile.")

def parse_alignment(aligned_file):
    """Parses the alignment file and identifies the reference and query sequences."""
    alignments = list(SeqIO.parse(aligned_file, "clustal"))

    reference = None
    query = None

    for alignment in alignments:
        if "OXA-48" in alignment.id:
            reference = alignment
        else:
            query = alignment

    if reference is None or query is None:
        raise ValueError("Could not identify reference (OXA-48) or query in the alignment file.")

    return reference, query

def load_secondary_structure_annotations(file_path):
    """Loads the secondary structure annotations from a CSV file."""
    ss_df = pd.read_csv(file_path)

    # Drop rows with missing 'Reference residue number' or 'Secondary Structure Annotation'
    ss_df = ss_df.dropna(subset=['Reference residue number', 'Secondary Structure Annotation'])

    # Ensure 'Reference residue number' is numeric and clean up invalid rows
    ss_df['Reference residue number'] = pd.to_numeric(ss_df['Reference residue number'], errors='coerce')
    ss_df = ss_df.dropna(subset=['Reference residue number'])  # Remove rows with NaN after conversion

    # Convert to dictionary, ensuring no NaN keys or values
    ss_dict = {
        int(row['Reference residue number']): row['Secondary Structure Annotation']
        for _, row in ss_df.iterrows()
        if pd.notna(row['Secondary Structure Annotation'])
    }
    return ss_dict

def map_and_save_csv(reference, query, output_file_col, ss_dict):
    """Maps residues between reference and query sequences and saves the result to a CSV file."""
    ref_index = 0
    query_index = 0
    suffix = ''

    data = [
        ['Reference residue number'],
        ['Reference Secondary Structure Annotation'],
        ['Reference residue name'],
        ['Query residue name'],
        ['Query original numbering'],
        ['Query Standard numbering'],
        ['Comments']
    ]

    for ref_base, query_base in zip(reference.seq, query.seq):
        if ref_base == "-" and query_base == "-":
            continue

        comment = ""
        new_number = ""

        if ref_base != "-":
            ref_index += 1
            suffix = ''

        if query_base != "-":
            query_index += 1
            if ref_base != "-":
                new_number = f"{ref_index}"
            else:
                if not suffix:
                    suffix = 'a'
                else:
                    suffix = chr(ord(suffix) + 1)
                new_number = f"{ref_index}{suffix}"
                comment = "Insertion"
        else:
            comment = "Gap" if ref_base == "-" else "Deletion"

        ss_annotation = ss_dict.get(ref_index, "") if ref_base != "-" else ""

        data[0].append(ref_index if ref_base != "-" else "-")
        data[1].append(ss_annotation)
        data[2].append(ref_base if ref_base != "-" else "-")
        data[3].append(query_base if query_base != "-" else "-")
        data[4].append(query_index if query_base != "-" else "-")
        data[5].append(new_number if query_base != "-" else "-")
        data[6].append(comment)

    transposed_data = zip(*data)
    cleaned_data = [[cell if cell != "NaN" else "" for cell in row] for row in transposed_data]
    with open(output_file_col, 'w', newline='') as csvfile_col:
        writer = csv.writer(csvfile_col)
        writer.writerows(cleaned_data)

    print(f"Mapped CSV file saved as {output_file_col} (column format).")

def print_mapping_table(reference, query, ss_dict):
    """Displays a mapping table between the reference and query sequences."""
    ref_index = 0
    query_index = 0
    suffix = ''

    data = {
        'Reference residue number': [],
        'Reference Secondary Structure Annotation': [],
        'Reference residue name': [],
        'Query residue name': [],
        'Query original numbering': [],
        'Query Standard numbering SAND': [],
        'Comments': []
    }

    for ref_base, query_base in zip(reference.seq, query.seq):
        if ref_base == "-" and query_base == "-":
            continue

        comment = ""
        new_number = ""

        if ref_base != "-":
            ref_index += 1
            suffix = ''

        if query_base != "-":
            query_index += 1
            if ref_base != "-":
                new_number = f"{ref_index}"
            else:
                if not suffix:
                    suffix = 'a'
                else:
                    suffix = chr(ord(suffix) + 1)
                new_number = f"{ref_index}{suffix}"
                comment = "Insertion"
        else:
            comment = "Gap" if ref_base == "-" else "Deletion"

        ss_annotation = ss_dict.get(ref_index, "") if ref_base != "-" else ""

        data['Reference residue number'].append(ref_index if ref_base != "-" else "-")
        data['Reference Secondary Structure Annotation'].append(ss_annotation)
        data['Reference residue name'].append(ref_base if ref_base != "-" else "-")
        data['Query residue name'].append(query_base if query_base != "-" else "-")
        data['Query original numbering'].append(query_index if query_base != "-" else "-")
        data['Query Standard numbering SAND'].append(new_number if query_base != "-" else "-")
        data['Comments'].append(comment)

    df = pd.DataFrame(data)
    df = df.fillna("")  # Replace NaN with empty strings
    display(df)

def cleanup_files():
    """Deletes temporary files created during the process."""
    temp_files = ["concatenated_output.fasta", "query.fasta", "query_ref_aligned.sto", "query_ref_aligned.aln"]
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)