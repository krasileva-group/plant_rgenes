plant_rgenes
============

Set of scripts to annotate Pfam domains and extract NLR plant immune receptors and their architectures

Our basic pipeline

0) Obtain protein sequences of species of interest and organise them into a directory. 

We follow the Phytozome organisation of `master_dir/species/annotation/species_version_proteins.fa` where each species is denoted by the first letter of the genus name and all letters in the species names, for example `Athaliana`

1) Pfam-based annotation of domains

Usage: bash run_pfam_scan.sh dir

Dependencies:
1. HMMER software (http://hmmer.janelia.org/) including pfam_scan.pl (part of HMMER) Move in same directory as this script or set path at command string.
2. Pfam database (http://pfam.xfam.org/)
3. File names should be consistent with Phytozome and include Species_*_protein.fa

2) Parsing the pfamscan output with K-parse_Pfam_domains_v3.1.pl 

The script parses the output of pfam_scan.pl

The script extracts all Domains for each proteins and removes redundant nested hits with larger e-values. Domains are printed out in the order of apprearance in the query. By default, Pfam_B domains are skipped.

usage: perl K-parse_Pfam_domains_v3.1.pl  <options>
-p|--pfam <pfamscan.out>
-e|--evalue <evalue cutoff>
-o|--output
-v|--verbose <T/F> default F. Display more information about each domain (start, stop, evalue)

We usually parse all pfam outputs of interest in parallel using `xargs` 

3) Identification of non-canonical NLR-ID domain combinations with K-parse_Pfam_domains_NLR-fusions-v2.2.pl

This script is configured to find any parsed pfam files in specified directory or its sub-directories

usage: perl K-parse_Pfam_domains_NLR-fusions-v2.2.pl <options>

-i|--indir directory for batch retrieval of input *pfamscan*.parsed.verbose files
-e|--evalue evalue cutoff for determining domain fusions [default 1e-3]
-o|--output output directory
-d|--db_description description of datasets used in the analyses [Organism Species_ID NCBI_taxon_ID Family Database Date_aquired Restrictions Version Common_Name Source Reference] for example of this dataset see Additional file 1 in Sarris et al BMC Biology 2016

The script will parse the output of K-parse_Pfam_domains-v3.1.pl

NLR proteins are identified based on the presence of NB-ARC domain. Fusions are identified based on the presence of non-NBS non-LRR domains with specified evalue cutoff (default 1e-3). It will generate several outputs:

Summary of the number of NLRs and NLR-IDs identified in each species (such as Additional file 2 in Sarris et al BMC Biology 2016)
Summary of integrated domains with species list for each domain (such as Additional file 3 in Sarris et al BMC Biology 2016)
Abundance list  of integrated domains (counted once for each family) that can be used to generate a Wordcloud (such as Figure 2 in Sarris et al BMC Biology 2016)
Contingency tables (per ID domain) for each species as well as for all species and Fisher's Exact left test 
 

 