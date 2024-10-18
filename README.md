
# MiniPlasmo

A Python-based tool for exploring and analyzing Plasmodium falciparum gene data retrieved from PlasmoDB and ToxoDB, including recodonization, genome alignment, and orthologue retrieval.

**Author**: Jakob Cronshagen  
**Date Created**: August 14, 2024

## Overview

MiniPlasmo is a Python-based application designed to facilitate the exploration and analysis of genetic data related to various species and strains of *Plasmodium* (the causative agent of malaria). This tool enables users to search for gene-related information, perform codon recoding based on species-specific codon usage, and analyze exonic and genomic sequences for primer design and other research applications.

The tool was originally developed during a period when PlasmoDB was expected to go offline due to funding issues. Since PlasmoDB's return, the tool has continued to evolve, incorporating additional features and improvements.

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
```

## Project Structure

The repository is structured as follows:

```
├── MiniPlasmo.py        # Main script file for running the tool
├── ref_files                        # Directory containing reference files
│   ├── 3D7.json                     # Example gene data JSON file for P. falciparum 3D7 strain
│   ├── PlasmoDB-68_*.fasta          # Genomic FASTA files for various strains
│   ├── groups_OrthoMCL-CURRENT.txt  # Orthogroup file for ortholog data
└── README.md                        # This file, providing an overview and instructions
```

## How to Use

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yourusername/miniplasmo.git
   ```

2. **Navigate to the Project Directory:**
   ```bash
   cd miniplasmo
   ```

3. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
4. **Download the reference files and copy them into the ref_files folder**
   See INSTRUCTIONS_ref_files_Download.txt in the ref_file folder
   
5. **Run the Tool:**
   ```bash
   python malaria_gene_info_tool.py
   ```

6. **Follow the Interactive Command-Line Instructions.**

## License

This project is licensed under the [MIT License](LICENSE).

## References

Amos B, Aurrecoechea C, Barba M, Barreto A, Basenko EY, Bażant W, Belnap R, Blevins AS, Böhme U, Brestelli J, Brunk BP, Caddick M, Callan D, Campbell L, Christensen MB, Christophides GK, Crouch K, Davis K, DeBarry J, Doherty R, Duan Y, Dunn M, Falke D, Fisher S, Flicek P, Fox B, Gajria B, Giraldo-Calderón GI, Harb OS, Harper E, Hertz-Fowler C, Hickman MJ, Howington C, Hu S, Humphrey J, Iodice J, Jones A, Judkins J, Kelly SA, Kissinger JC, Kwon DK, Lamoureux K, Lawson D, Li W, Lies K, Lodha D, Long J, MacCallum RM, Maslen G, McDowell MA, Nabrzyski J, Roos DS, Rund SSC, Schulman SW, Shanmugasundram A, Sitnik V, Spruill D, Starns D, Stoeckert CJ, Tomko SS, Wang H, Warrenfeltz S, Wieck R, Wilkinson PA, Xu L, Zheng J. (2022). VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics resource center. Nucleic acids research, 50(D1), D898–D911. https://doi.org/10.1093/nar/gkab929

Fischer, S., Brunk, B. P., Chen, F., Gao, X., Harb, O. S., Iodice, J. B., Shanmugam, D., Roos, D. S., & Stoeckert, C. J., Jr (2011). Using OrthoMCL to assign proteins to OrthoMCL-DB groups or to cluster proteomes into new ortholog groups. Current protocols in bioinformatics, Chapter 6, Unit–6.12.19. https://doi.org/10.1002/0471250953.bi0612s35

