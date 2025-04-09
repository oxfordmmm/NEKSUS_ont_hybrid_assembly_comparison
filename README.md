# Nanopore long-read only genome assembly is complete and accurate for surveillance of Enterobacterales bloodstream infections in England

## Summary 

This repo contains the extended Methods (below) and Nextflow pipeline (assembly_comparison_pipeline_v2.tar.gz) used to compare Nanopore long-read only with hybrid bacterial genome assembly for Enterobacterales. This analysis was a pilot study for the NEKSUS study (National *E. coli* and *Klebsiella* Bloodstream Infection and CPE Surveillance Study), and was presented at the 2025 Congress of the European Society of  Clinical Microbiology and Infectious Diseases as an ePoster. 


## Contents
- Extended methods
- Acknowledgements
- Citation
- References
 
## Extended Methods

### Isolate collection

9 English hospitals from 8 English regions were recruited through a partnership between the University of Oxford and UK Health Secuirty Agency (UKHSA), and consecutive bloodstream infection (BSI) and rectal screening isolates were collected between October 2023 and March 2024. Isolate collection was part of routine clinical practice under existing UKHSA permissions for extended surveillance. Isolates were prepared by participating sites and shipped to Oxford, catalogued, and stored in brain-heart infusion (BHI) broth with 10% Glycerol at -70◦C unitl preparation for DNA extraction and sequencing by an external provider.

### DNA extraction & sequencing

DNA extraction, library preparations, and long- and short-read sequencing were conducted at GENEWIZ Germany GmbH (Leipzig, Germany). Briefly, DNA extracts were both long-read (PromethION using R10.4.1/kit 14 flowcells/chemistry [Oxford Nanopore Technologies]) and short-read sequenced (Illumina NovaSeq X Plus, 2x150 paired-end reads). 

### Basecalling 

Long-read sequences were basecalled with Dorado (v5.0.0 super high accuracy (sup) simplex DNA models). The genome assembly/annotation pipeline (assembly_comparison_pipeline_v2.tar.gz) was written in Nextflow (v24.04.3) based on an erlier Nanopore assembly comparison [1]. Raw-read and assembly quality was assessed with SeqKit (v2.9.0). Long-reads were subsampled to 40x and short-reads to 100x with Rasusa (v2.1.0), with 5.6Mb estimated genome size.  I used four long-read-only, and 4 hybrid methods to generate assemblies. The long-read-only methods I evaluated were Flye (v2.9.5) without and with one round of Medaka (v2.0.1) polishing with the subsampled or full-set of reads, or Hybracter (‘long’ mode; v0.8.0). I generated hybrid assemblies with Unicycler (v0.5.0) in normal and bold mode for short-read-first assemblies, and Hybracter (hybrid mode), or Flye-Polypolish (v0.6.0) for long-read-first hybrid assemblies. Dnadiff from MUMmer4 (v4.0.0) was used to align all assemblies against the Flye-Polypolish hybrid assembly (chosen as the “gold-standard”) for each isolate. mlst (v2.23.0) was used to identify multi locus sequence type (MLST). Key resistance, stress and virulence genes were identified using AMRFinder Plus (v4.0.3) with the species flag inferred from Kraken2 (v2.1.3). Pilot analyses were carried out on the Oxford Biomedical Research Computer Cluster. Statistical analysis and visualisation were done in R (v4.4.1; R Core Team 2024).

### Assembly

### Annotation 

### Statistical analysis and visualisation

## Acknowledgements

This work was funded by the UK Health Security Agency (UKHSA) and supported by the National Institute for Health Research (NIHR) Health Protection Research Unit in Healthcare Associated Infections and Antimicrobial Resistance (NIHR200915), a partnership between the UKHSA and the University of Oxford, the NIHR Oxford Biomedical Research Centre (BRC) and the UKHSA [UKHSA PhD Funding Competition]. The cloud infrastructure used in this study was donated by Oracle Corporation.


## Citation

Nagy, D., Pennetta, V., Rodger, G., Hopkins, K., Jones, C., The NEKSUS Study Group, Hopkins, S., Crook, D., Walker, A.S., Robotham, J., Hopkins, K.L., Ledda, A., Williams, D., Hope, R., Brown, C.S., Stoesser, N., Lipworth, S. (2025). *Nanopore long-read only genome assembly is complete and accurate for surveillance of Enterobacterales bloodstream infections in England.* [Conference ePoster]. ESCMID Global 2025, Vienna, Austria. https://registration.escmid.org//AbstractList.aspx?e=30&header=0&preview=1&aig=-1&ai=33860 

## References
