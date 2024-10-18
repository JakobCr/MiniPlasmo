# MiniPlasmo
A Python-based tool for exploring and analyzing Plasmodium falciparum gene data retrieved from PlasmoDB and ToxoDB, including recodonization, genome alignment, and orthologue retrieval.

**Author**: Jakob Cronshagen

**Date Created**: August 14, 2024

## Overview

MiniPlasmo is a Python-based application designed to facilitate the exploration and analysis of genetic data related to various species and strains of *Plasmodium* (the causative agent of malaria). This tool enables users to search for gene-related information, perform codon recoding based on species-specific codon usage, and analyze exonic and genomic sequences for primer design and other research applications.

The tool was originally developed during a period when PlasmoDB was temporarily offline due to funding issues. Since PlasmoDB's return, the tool has continued to evolve, incorporating additional features and improvements.

The tool utilizes external resources like PlasmoDB and ToxoDB and integrates sequence alignment and recodonization capabilities, providing users with comprehensive information for studying the genomic data of these organisms.

## Key Features

1. **Search and Explore Gene Data:**
   - Search for gene information using *gene ID*, *gene name*, or *gene product*.
   - Access important gene-related categories, including protein sequences, orthologs, gene length, and more.

2. **Codon Recodonization:**
   - Recodonize DNA sequences based on *Plasmodium falciparum* codon usage to optimize for genetic engineering and other purposes.
   - Handle both manual DNA sequence input and gene ID-based recodonization using species-specific codon frequency tables.

3. **Gene-to-Genome Alignment:**
   - Align coding DNA sequences (CDS) to the genomic context, highlighting exons and identifying intron-exon boundaries.
   - Useful for primer design and analysis of gene structure.

4. **Orthologue Retrieval:**
   - Find and display orthologues from a pre-downloaded orthogroup file for comparative genomic studies.

5. **Multiple Strain Support:**
   - Supports multiple strains of *Plasmodium falciparum*, *Plasmodium vivax*, and *Toxoplasma gondii*.
   - Easily switch between strains using a predefined file mapping system.

6. **Interactive Menu Interface:**
   - The tool uses a simple command-line interface with easy-to-navigate menus, including searching genes, recodonizing sequences, and switching between different genomic reference files.

## Dependencies

The tool uses several external libraries:

- **Biopython**: For sequence handling, translation, and alignment (`Bio.Seq`, `Bio.Align`, `Bio.SeqIO`).
- **Colorama**: For colored terminal output, enhancing the user experience with visual cues for matched and unmatched bases in alignments.
- **JSON**: To handle genomic and gene data stored in JSON files.
- **Regular Expressions (re)**: For parsing genomic coordinates and location text.

You can install dependencies using:
```bash
pip install biopython colorama
