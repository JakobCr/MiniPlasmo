# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:35:11 2024

@author: Jakob Cronshagen
"""
import json
import os
import sys
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
import re
from colorama import Fore, Style, init
global gene_data  # Ensure gene_data is global
global unspliced_file
genome_file = None
strain = None

# Initialize colorama
init(autoreset=True)

# Get Python script path
def get_resource_path(relative_path):
    # Get the absolute path to the resource, works for both bundled and non-bundled execution.
    try:
        base_path = sys._MEIPASS  # PyInstaller's temp folder for bundled files
    except AttributeError:
        base_path = os.path.abspath(".")  # Script's directory

    return os.path.join(base_path, relative_path)


# file containing literature informations about all genes of the 3D7 strain
lit_file_3D7 = get_resource_path(os.path.join('ref_files', '3D7_Literature'))

# Define a mapping of user choices to file details; strain ID in case needed in further development
# Path to the reference files file paths:
    # Files contianing information about every gene in the genome of the corresponding strain downloaded from plasmoDB and ToxoDB
file_mapping = {
    '1':  {'name': 'Plasmodium falciparum 3D7', 'input_file': get_resource_path(os.path.join('ref_files', '3D7.json')), 'strain': 1, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_Pfalciparum3D7_Genome.fasta'))},
    '2':  {'name': 'Plasmodium falciparum IT4', 'input_file': get_resource_path(os.path.join('ref_files', 'IT4.json')), 'strain': 2, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumIT_Genome.fasta'))},
    '3':  {'name': 'Plasmodium falciparum 7G8', 'input_file': get_resource_path(os.path.join('ref_files', '7G8.json')), 'strain': 4, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_Pfalciparum7G8_Genome.fasta'))},
    '4':  {'name': 'Plasmodium falciparum CD01', 'input_file': get_resource_path(os.path.join('ref_files', 'CD01.json')), 'strain': 5, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumCD01_Genome.fasta'))},
    '5':  {'name': 'Plasmodium falciparum DD2', 'input_file': get_resource_path(os.path.join('ref_files', 'DD2.json')), 'strain': 6, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumDd2_Genome.fasta'))},
    '6':  {'name': 'Plasmodium falciparum GA01', 'input_file': get_resource_path(os.path.join('ref_files', 'GA01.json')), 'strain': 7, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumGA01_Genome.fasta'))},
    '7':  {'name': 'Plasmodium falciparum GB4', 'input_file': get_resource_path(os.path.join('ref_files', 'GB4.json')), 'strain': 8, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumGB4_Genome.fasta'))},
    '8':  {'name': 'Plasmodium falciparum GN01', 'input_file': get_resource_path(os.path.join('ref_files', 'GN01.json')), 'strain': 9, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumGN01_Genome.fasta'))},
    '9':  {'name': 'Plasmodium falciparum HB3', 'input_file': get_resource_path(os.path.join('ref_files', 'HB3.json')), 'strain': 10, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumHB3_Genome.fasta'))},
    '10': {'name': 'Plasmodium falciparum KE01', 'input_file': get_resource_path(os.path.join('ref_files', 'KE01.json')), 'strain': 11, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumKE01_Genome.fasta'))},
    '11': {'name': 'Plasmodium falciparum KH01', 'input_file': get_resource_path(os.path.join('ref_files', 'KH01.json')), 'strain': 12, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumKH01_Genome.fasta'))},
    '12': {'name': 'Plasmodium falciparum KH02', 'input_file': get_resource_path(os.path.join('ref_files', 'KH02.json')), 'strain': 13, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumKH02_Genome.fasta'))},
    '13': {'name': 'Plasmodium falciparum ML01', 'input_file': get_resource_path(os.path.join('ref_files', 'ML01.json')), 'strain': 14, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumML01_Genome.fasta'))},
    '14': {'name': 'Plasmodium falciparum NFS4', 'input_file': get_resource_path(os.path.join('ref_files', 'NFS4.json')), 'strain': 15, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumNFS4_Genome.fasta'))},
    '15': {'name': 'Plasmodium falciparum NF135', 'input_file': get_resource_path(os.path.join('ref_files', 'NF135.json')), 'strain': 16, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumNF135C10_Genome.fasta'))},
    '16': {'name': 'Plasmodium falciparum NF166', 'input_file': get_resource_path(os.path.join('ref_files', 'NF166.json')), 'strain': 17, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumNF166_Genome.fasta'))},
    '17': {'name': 'Plasmodium falciparum SD01', 'input_file': get_resource_path(os.path.join('ref_files', 'SD01.json')), 'strain': 22, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumSD01_Genome.fasta'))},
    '18': {'name': 'Plasmodium falciparum SN01', 'input_file': get_resource_path(os.path.join('ref_files', 'SN01.json')), 'strain': 23, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumSN01_Genome.fasta'))},
    '19': {'name': 'Plasmodium falciparum TG01', 'input_file': get_resource_path(os.path.join('ref_files', 'TG01.json')), 'strain': 24, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PfalciparumTG01_Genome.fasta'))},
    '20': {'name': 'Plasmodium berghei ANKA', 'input_file': get_resource_path(os.path.join('ref_files', 'PbANKA.json')), 'strain': 3, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PbergheiANKA_Genome.fasta'))},
    '21': {'name': 'Plasmodium vivax P01', 'input_file': get_resource_path(os.path.join('ref_files', 'PvP01.json')), 'strain': 18, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PvivaxP01_Genome.fasta'))},
    '22': {'name': 'Plasmodium vivax PAM', 'input_file': get_resource_path(os.path.join('ref_files', 'PvPAM.json')), 'strain': 19, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PvivaxPAM_Genome.fasta'))},
    '23': {'name': 'Plasmodium vivax W1', 'input_file': get_resource_path(os.path.join('ref_files', 'PvW1.json')), 'strain': 20, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PvivaxPvW1_Genome.fasta'))},
    '24': {'name': 'Plasmodium vivax X', 'input_file': get_resource_path(os.path.join('ref_files', 'PvX.json')), 'strain': 21, 'genome_file': get_resource_path(os.path.join('ref_files', 'PlasmoDB-68_PvivaxSal1_Genome.fasta'))},
    '25': {'name': 'Toxoplasma gondii ME49', 'input_file': get_resource_path(os.path.join('ref_files', 'TgME49.json')), 'strain': 25, 'genome_file': get_resource_path(os.path.join('ref_files', 'ToxoDB-68_TgondiiME49_Genome.fasta'))},
}

# Plasmodium falciparum codon usage table (frequency per thousand codons)
pf_codon_usage = {
    'F': [('TTT', 83.7), ('TTC', 16.3)],
    'L': [('TTA', 47.4), ('TTG', 10.5), ('CTT', 8.7), ('CTC', 1.8), ('CTA', 6.0), ('CTG', 1.5)],
    'I': [('ATT', 36.0), ('ATC', 6.3), ('ATA', 50.2)],
    'M': [('ATG', 100.0)],
    'V': [('GTT', 15.1), ('GTC', 2.4), ('GTA', 15.6), ('GTG', 4.8)],
    'S': [('TCT', 14.6), ('TCC', 5.1), ('TCA', 16.6), ('TCG', 3.0), ('AGT', 20.4), ('AGC', 3.9)],
    'P': [('CCT', 7.8), ('CCC', 2.1), ('CCA', 9.0), ('CCG', 0.9)],
    'T': [('ACT', 10.5), ('ACC', 4.8), ('ACA', 21.6), ('ACG', 3.8)],
    'A': [('GCT', 8.1), ('GCC', 2.1), ('GCA', 8.3), ('GCG', 1.1)],
    'Y': [('TAT', 50.8), ('TAC', 6.2)],
    'H': [('CAT', 20.6), ('CAC', 3.5)],
    'Q': [('CAA', 23.9), ('CAG', 3.7)],
    'N': [('AAT', 124.1), ('AAC', 20.1)],
    'K': [('AAA', 95.7), ('AAG', 21.6)],
    'D': [('GAT', 55.9), ('GAC', 8.7)],
    'E': [('GAA', 60.9), ('GAG', 10.3)],
    'C': [('TGT', 15.3), ('TGC', 2.3)],
    'W': [('TGG', 100.0)],
    'R': [('CGT', 3.0), ('CGC', 0.4), ('CGA', 2.4), ('CGG', 0.3), ('AGA', 16.0), ('AGG', 4.4)],
    'G': [('GGT', 11.7), ('GGC', 1.3), ('GGA', 12.4), ('GGG', 2.8)],
    '*': [('TAA', 91.0), ('TGA', 27.0), ('TAG', 13.0)]
}

# List of most important categories
important_categories = [
    'gene_source_id', 'gene_name', 'gene_product', 'gene_previous_ids', 'gene_orthomcl_name',
    'gene_type', 'molecular_weight', 'isoelectric_point', 'tm_count', 'signalp_peptide', 'exon_count',
    'genomic_sequence_length', 'cds_length', 'transcript_length', 'protein_length', 'location_text',
    'interpro_description', 'gene_paralog_number', 'gene_ortholog_number', 'annotated_go_function', 'annotated_go_process',
    'cds', 'protein_sequence'
]

greeting_message = r"""
        _       _   ___ _                           
  /\/\ (_)_ __ (_) / _ \ | __ _ ___ _ __ ___   ___  
 /    \| | '_ \| |/ /_)/ |/ _` / __| '_ ` _ \ / _ \ 
/ /\/\ \ | | | | / ___/| | (_| \__ \ | | | | | (_) |
\/    \/_|_| |_|_\/    |_|\__,_|___/_| |_| |_|\___/            
..............::..............
.......:=%%%%%%%%%%%%#:.......
.....+%%%%%#%#######%#%%#=.... 
...=%%%%%%%###%#########%%%...  
.:#%%%%%%%%##%%**++++*#%%%%%=.  
.+%%%%####%%+*---::::=++**%%%:  
.%%%%#%%%%*=+##*-::*%*+=--#%%=  
:%%%#####++=+%%%#:-%%*=-::#%%+  
.#%%%####%**=-=-=-::##***#%%%-  
.-%%%######%*=--%%+::-=#%%#%*.
..:%%###########%%%%######%*..  Data provided by VEuPathDB (PlasmoDB.com and ToxoDB.com)
....+%%#################%#:...  for non-commercial use only
......:#%############%%*:.....  version 0.31
..........:+*#%%%#*=..........  created by Jakob Cronshagen       

Welcome to the genome analysis tool for gathering information about malaria genes
    """

# Function to greet the user with an ASCII art
def greet_user():
    print(greeting_message)
    input("Press Enter to continue...")

# Function to load JSON data from the selected file
def load_json_data(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Load the orthology data into a dictionary
def load_orthology_data(ortho_file):
    ortho_data = {}
    with open(ortho_file, 'r') as file:
        for line in file:
            if line.strip():
                ortho_name, genes = line.split(':')
                ortho_data[ortho_name.strip()] = genes.strip().split()
    return ortho_data

# File containing all orthogroups
ortho_file = get_resource_path(os.path.join('ref_files', 'groups_OrthoMCL-CURRENT.txt'))
orthology_data = load_orthology_data(ortho_file)

# Function to get information by Gene ID and category; "Parser"
def get_gene_info(gene_id, category):
    attributes = gene_data.get(gene_id)
    if attributes:
        return attributes.get(category, None)
    return None

# Function to display important categories for a given Gene ID
def show_important_categories(gene_id):
    print(f"\nImportant categories for Gene ID {gene_id}:")
    result = {}
    for category in important_categories:
        info = get_gene_info(gene_id, category)
        if info is not None:
            print(f"{category}: {info}")
            result[category] = info
        else:
            result[category] = 'N/A'
    return result


# Function to search by gene_name, previous_id, or part of gene_product
def search_gene(query):
    results = []
    query_lower = query.lower()
    for gene_id, attributes in gene_data.items():
        if (query_lower in (attributes.get('gene_name', '') or '').lower() or
            query_lower in (attributes.get('gene_previous_ids', '') or '').lower() or
            query_lower in (attributes.get('gene_product', '') or '').lower() or
            query_lower == gene_id.lower()):
            
            # Include all attributes of the gene in the result
            gene_info = {'gene_source_id': gene_id}
            gene_info.update(attributes)  # Add all attributes to the result
            results.append(gene_info)
    return results


# Recodonize based on Plasmodium falciparum codon usage (or standard codon table)
def recodonize(dna_seq):
    # Use the codon table to optimize for Plasmodium falciparum usage
    recodonized_seq = ''

    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            # Use the codon to find the corresponding amino acid
            amino_acid = Seq(codon).translate(table="Standard", to_stop=False)

            if amino_acid != '*':  # Ignore stop codons for recodonization
                # Get all codons for this amino acid based on Plasmodium falciparum usage
                if amino_acid and amino_acid in pf_codon_usage:
                    # Sort codons by decreasing frequency to prioritize the most common ones
                    sorted_codons = sorted(pf_codon_usage[amino_acid], key=lambda x: x[1], reverse=True)
                    # Select the most frequent codon that differs from the original
                    for candidate_codon, _ in sorted_codons:
                        if candidate_codon != codon:
                            recodonized_seq += candidate_codon
                            break
                    else:
                        # If no alternative codon is found, use the most frequent one
                        recodonized_seq += sorted_codons[0][0]  # The codon with the highest frequency
            else:
                # Retain stop codons
                recodonized_seq += codon
        else:
            recodonized_seq += codon  # Retain partial codons if any

    return recodonized_seq


# Function to handle recodonization and display relevant information
def recodonize_menu():
    while True:
        print("\nEnter 'exit' to go back to the main menu.")
        
        # Ask the user if they want to input a Gene ID or a DNA sequence
        input_type = input("Do you want to enter a Gene ID or a DNA sequence? (id/seq): ").strip().lower()
        
        if input_type == 'exit':
            return True
        
        if input_type == 'id':
            choice = input("Enter a Gene ID to recodonize: ").strip()
            
            if choice.lower() == 'exit':
                return True
            
            if choice in gene_data:
                # Get the original CDS sequence
                cds_seq = get_gene_info(choice, 'cds')
                if cds_seq:
                    print(f"\nOriginal CDS sequence for Gene ID {choice}: {cds_seq}")
                    
                    # Recodonize the CDS sequence
                    recodonized_cds_seq = recodonize(cds_seq)
                    
                    # Translate the recodonized CDS sequence back to amino acids
                    recodonized_amino_acid_seq = Seq(recodonized_cds_seq).translate()
                    
                    # Display the original CDS, recodonized CDS, and recodonized amino acid sequence
                    print(f"\nRecodonized CDS sequence (Plasmodium falciparum optimized): {recodonized_cds_seq}")
                    print(f"\nRecodonized CDS translated to Amino Acids: {recodonized_amino_acid_seq}\n")
                    
                else:
                    print(f"No CDS sequence found for Gene ID {choice}.")
            else:
                print(f"Gene ID {choice} not found.")
        
        elif input_type == 'seq':
            dna_seq = input("Enter a DNA sequence to recodonize: ").strip().upper()
            
            if dna_seq.lower() == 'exit':
                break
            
            # Check if the DNA sequence length is divisible by 3
            if len(dna_seq) % 3 != 0:
                print("The DNA sequence length is not divisible by 3. Please enter a valid sequence.")
            else:
                # Display the original DNA sequence
                print(f"\nOriginal DNA sequence: {dna_seq}")
                
                # Recodonize the DNA sequence
                recodonized_dna_seq = recodonize(dna_seq)
                
                # Translate the recodonized DNA sequence to amino acids
                recodonized_amino_acid_seq = Seq(recodonized_dna_seq).translate()
                
                # Display the original, recodonized, and translated sequences
                print(f"\nRecodonized DNA sequence (Plasmodium falciparum optimized): {recodonized_dna_seq}")
                print(f"\nRecodonized DNA sequence translated to Amino Acids: {recodonized_amino_acid_seq}\n")
        
        else:
            print("Invalid input. Please enter 'id' or 'sequence'.")


# Function to extract a specific region from a genomic record
def extract_genomic_region(chromosome_id, start, end, strand='+'):
    # Parse the genome FASTA file
    genome_records = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    # Find the chromosome/contig record
    if chromosome_id in genome_records:
        record = genome_records[chromosome_id]
        # Extract the region -500 and +500 of transcripted region
        region_seq = record.seq[start-50:end+500]
        
        # Check strand orientation
        if strand == '-':
            # If the strand is negative, reverse complement the sequence
            region_seq = region_seq.reverse_complement()
        
        # Return the extracted sequence
        return region_seq
    else:
        raise ValueError(f"Chromosome {chromosome_id} not found in genome FASTA.")

# Function to parse the location text
def parse_location(location_text):
    # Use a regular expression to extract chromosome, start, end, and strand
    match = re.match(r'(?P<chromosome>[\w\d_]+):(?P<start>\d+)\.\.(?P<end>\d+)\((?P<strand>[+-])\)', location_text)
    if match:
        chromosome = match.group('chromosome')
        start = int(match.group('start'))
        end = int(match.group('end'))
        strand = match.group('strand')
        return chromosome, start, end, strand
    else:
        raise ValueError(f"Invalid location format: {location_text}")

# Shows the CDS in the context of the genome sequence to show intron sequence, upstream und downstream sequence; Helpfull for primer design
def align_cds_to_genome(gene_id, cds_sequence, location_text):
    chromosome, start_position, end_position, strand_orientation = parse_location(location_text)
    extracted_sequence = extract_genomic_region(chromosome, start_position, end_position, strand_orientation)
    
    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2   # Match score
    aligner.mismatch_score = -1000  # Mismatch penalty, super high as there shouldn´t be any mismatches
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    # Forces global aligner in semi global alignment
    aligner.left_gap_score = 0    # No penalty for starting gaps
    aligner.left_extend_gap_score = 0
    aligner.right_gap_score = 0
    aligner.right_extend_gap_score = 0
    alignments = aligner.align(cds_sequence, extracted_sequence)
    
    # Get the best alignment
    best_alignment = alignments[0]
    aligned_cds = best_alignment[0]  # Access the first sequence
    aligned_genomic = best_alignment[1]  # Access the second sequence
    # Highlight matched bases
    highlighted_genomic = []

    for i in range(len(aligned_genomic)):
        if aligned_genomic[i] == aligned_cds[i] and aligned_genomic[i] != '-':
            # Use green for matched bases
            highlighted_genomic.append(Fore.GREEN + aligned_genomic[i].upper() + Style.RESET_ALL)
        else:
            # Use red for unmatched bases
            highlighted_genomic.append(Fore.RESET + aligned_genomic[i].lower() + Style.RESET_ALL)
    print(f"\nShowing exons (capital letters) for Gene ID {gene_id}:")
    return print(''.join(highlighted_genomic))



# Function to display all categories for a given Gene ID
def show_all_categories(gene_id):
    print(f"\nAll categories for Gene ID {gene_id}:")
    attributes = gene_data.get(gene_id)
    if attributes:
        for category, info in attributes.items():
            print(f"{category}: {info}")
    else:
        print("No data found.")

# Function to display available categories
def show_categories():
    # Get the first gene entry in gene_data to access its attributes
    first_gene_id = next(iter(gene_data))
    
    # Access the attributes of the first gene to retrieve available categories
    categories = sorted(gene_data[first_gene_id].keys())
    
    print("\nAvailable categories (sorted alphabetically):")
    for category in categories:
        print(f"{category}", category)
        

# Function to show orthologues with the help of the orthogroup file
def show_orthologues(gene_id):
    # Name of the orthogroup is available in the gene infos
    orthomcl_name = get_gene_info(gene_id, 'gene_orthomcl_name')
    if orthomcl_name:
        # Gets all orthologues from for the orthogroup, if not available returns an empty list
        orthologues = orthology_data.get(orthomcl_name, [])
        if orthologues:
            print(f"\nOrthologues for {orthomcl_name}:")
            for orthologue in orthologues:
                print(orthologue)
        else:
            print(f"No orthologues found for {orthomcl_name}.")
    else:
        print(f"No orthologue data found for gene ID {gene_id}.")




def show_literature(gene_id):
    
    pass




def search_hit_options(gene_id):
    while True:
        command = input(
    "\nEnter:\n"
    "category name to explore more details\n"
    "'cat' to list category names\n"
    "'all' to show all categories\n"
    "'ortho' to show orthologues\n"
    "'exon' to align and highlight exons in the genome sequence\n"
    "'new' to search again\n"
    "'exit' to quit: "
).strip().lower()
        
        if command == 'all':
            show_all_categories(gene_id)
        elif command == 'ortho':
            show_orthologues(gene_id)
        elif command == 'cat':
            show_categories()
        elif command == 'exit':
            return False
        elif command == 'new':
            return True
        elif command == 'lit':
            pass
        elif command == 'exon':
            # Print the CDS sequence first
            cds_sequence = gene_data[gene_id].get('cds')
            location_text = gene_data[gene_id].get('location_text')
            if cds_sequence:
                print(f"CDS Sequence for {gene_id}: {cds_sequence}")
            else:
                print(f"No CDS sequence found for Gene ID {gene_id}.")
                continue
            # Align the CDS to the matched genomic sequence
            align_cds_to_genome(gene_id, cds_sequence, location_text)
        elif command in important_categories:
            info = get_gene_info(gene_id, command)
            if info:
                print(f"{command}: {info}")
            else:
                print(f"No data found for category '{command}'.")
        else:
            print(f"Invalid input '{command}'. Type 'help' to see available categories.")

# Function to handle picking a number from multiple matches
def pick_gene_from_multiple_matches(search_results):
    while True:
        print("\nMultiple matches found:")
        for idx, result in enumerate(search_results, start=1):
            print(f"{idx}. Gene ID: {result['gene_source_id']}")
            print(f"   Gene Name: {result['gene_name']}")
            print(f"   Gene Product: {result['gene_product']}")
            print(f"   {result['gene_previous_ids']}\n")
        # Get user input
        choice = input("Enter the listed number of the gene you want to explore, or 'exit' to go back: ").strip().lower()

        if choice == 'exit':  # Exit condition
            return None

        # Try to convert the input to an integer
        try:
            choice = int(choice)  # Convert the input to an integer
            if 1 <= choice <= len(search_results):  # Validate if the choice is within the range
                return search_results[choice - 1]['gene_source_id']  # Return the selected gene ID
            else:
                print("Invalid choice. Please enter a number from the list.")
        except ValueError:
            print("Invalid input. Please enter a number or type 'exit' to go back.")

def search_menu():
    while True:
        query = input("Enter 'exit' to return\nEnter Gene ID(s) (separated by commas), Gene Name, Previous ID, or part of Gene Product: ").strip()
        if query.lower() == 'exit':
            break
        if query.lower() == '' or query == ',':
            continue
        search_querys = [q.strip() for q in query.split(',')]
        results = {}
        
        # Handling multiple Gene IDs
        if len(search_querys) > 1:
            for search_query in search_querys:
                search_results = search_gene(search_query)
                
                if not search_results:
                    print(f"No matches found for '{search_query}'.")
                    continue
                
                if len(search_results) == 1:
                    gene_id = search_results[0]['gene_source_id']
                    results[gene_id] = show_important_categories(gene_id)
                else:
                    print(f"\nMultiple matches found for '{search_query}':")
                    chosen_gene_id = pick_gene_from_multiple_matches(search_results)
                    if chosen_gene_id:
                        results[chosen_gene_id] = show_important_categories(chosen_gene_id)
        # Handling a single Gene ID or search term
        elif len(search_querys) == 1:
            search_query = search_querys[0]
            search_results = search_gene(search_query)
            
            if not search_results:
                print(f"No matches found for '{search_query}'.")
                continue
            
            if len(search_results) == 1:
                gene_id = search_results[0]['gene_source_id']
                results[gene_id] = show_important_categories(gene_id)
                if not search_hit_options(gene_id):  # If search_hit_options returns False, exit the loop
                    break
            else:
                print(f"\nMultiple matches found for '{search_query}':")
                chosen_gene_id = pick_gene_from_multiple_matches(search_results)
                if chosen_gene_id:
                    results[chosen_gene_id] = show_important_categories(chosen_gene_id)
                    search_hit_options(chosen_gene_id)

# Second layer menu
def main_menu_2():
    global genome_file
    while True:
        print("\nMain Menu:")
        print("1. Search Genes")
        print("2. Recodonize DNA Sequences")
        print("3. Change reference genome")
        choice = input("Please choose an option (1-3): ").strip()
        
        if choice == '1':
            search_menu()
        elif choice == '2':
            recodonize_menu()
        elif choice == '3' or choice.lower() == 'exit':
            break
        else:
            print("Invalid option. Please select 1-3.")

# Main menu function
def main_menu():
    greet_user()
    global gene_data  # Ensure gene_data is global
    global strain
    global genome_file  # Ensure genome_sequences is global
    
    
    while True:
        # Display the menu options
        print("\nChoose reference genome:")
        for key, value in file_mapping.items():
            print(f"{key}. {value['name']}")
    
        file_choice = input("Please choose a reference file: ").strip()
    
            # Handle the user selection
        if file_choice in file_mapping:
            selected_file = file_mapping[file_choice]
            strain = selected_file['strain']
            gene_data = {record['attributes']['gene_source_id']: record['attributes'] for record in load_json_data(selected_file['input_file'])['records']}
            genome_file = selected_file['genome_file']
            main_menu_2()
        else:
            print("Invalid choice. Please select a valid option.")


if __name__ == "__main__":
    main_menu()

