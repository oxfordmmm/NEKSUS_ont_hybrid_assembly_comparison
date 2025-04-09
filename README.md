# Nanopore long-read only genome assembly is complete and accurate for surveillance of Enterobacterales bloodstream infections in England

## Summary 

This repo contains the extended Methods (below) and Nextflow pipeline (nanopore_assembly_comparison_pipeline_v2.tar.gz) used to compare Nanopore long-read only with hybrid bacterial genome assembly for Enterobacterales. This analysis was a pilot study for the NEKSUS study (National *E. coli* and *Klebsiella* Bloodstream Infection (BSI) and Carbapenemase-Producing Enterobacterales (CPE) UK Surveillance Study), and was presented at the 2025 Congress of the European Society of  Clinical Microbiology and Infectious Diseases as an ePoster. 


## Contents
- Extended methods
- Acknowledgements
- Citation
- References
 
## Extended Methods
### Isolate collection
9 English hospitals from 8 English regions were recruited through a partnership between the University of Oxford and UK Health Secuirty Agency (UKHSA). Consecutive bloodstream infection (BSI) and positive rectal Carbapenemase-Producing Enterobacterales (CPE) screening isolates were collected between October 2023 and March 2024. The first 92 Enterobacterales isolates collected, from 3 English regions, were included in the pilot analysis, predominantly *Escherichia coli* and *Klebsiella pneumoniae*. Isolate collection was part of routine clinical practice under existing UKHSA permissions for extended surveillance. Isolates were prepared by participating sites and shipped to Oxford, catalogued, and stored in brain-heart infusion (BHI) broth with 10% Glycerol at -70◦C unitl preparation for DNA extraction and sequencing by an external provider. 

### DNA extraction & sequencing
DNA extraction, library preparations, and long- and short-read sequencing were conducted at GENEWIZ Germany GmbH (Leipzig, Germany). Briefly, DNA extracts were both long-read (PromethION using R10.4.1/kit 14 flowcells/chemistry [Oxford Nanopore Technologies]) and short-read sequenced (Illumina NovaSeq X Plus, 2x150 paired-end reads). 

### Bioinformatics 
#### Basecalling
Long-read sequences were basecalled with Dorado (v5.0.0 super high accuracy simplex DNA models). The genome assembly/annotation pipeline (nanopore_assembly_comparison_pipeline_v2.tar.gz) was written in Nextflow (v24.04.3) based on an earlier Nanopore assembly comparison<sup>1</sup>. 
#### Quality control
Raw-read and assembly quality was assessed with SeqKit<sup>2</sup> v2.9.0.
#### Subsampling
Long-reads were subsampled to 60x using the built-in subsampling and genome size estimation scripts from Autocycler<sup>3</sup> v0.2.1, and short-reads to 100x with Rasusa<sup>4</sup> v2.1.0. 
#### Assembly
Genomes were assembled using three long-read only assemblers (Flye<sup>5</sup> v2.9.5, Hybracter<sup>6</sup> long v0.11.2, and the consensus assembler Autocycler<sup>3</sup> v0.2.1), and three hybrid assemblers (Hybracter<sup>6</sup> hybrid, Unicycler<sup>7</sup> v0.5.1 normal and bold modes). The input long-read assemblies used for Autocycler were 4 assemblies each of Canu v2.2, Flye, Raven v1.8.3, Miniasm v0.3, and Hybracter long, where each of the 4 assemblies was derived from an independently subsampled set of reads. The Flye an Hybracter long assemblies from the first subsampled reads set was used in downstream analyses.
#### Polishing
Three polishing modalities were investigated:
- long-read polishing with one round of Medaka<sup>8</sup> v2.0.1 with subsampled reads
- long-read polishing with one round of Medaka<sup>8</sup> with a full-set of reads
- short-read polishing with Polypolish<sup>9</sup> v0.6.0 and Pypolca<sup>10</sup> v0.3.1 
#### Assembly evaluation
Assembly completeness was assessed by the proportion of chromosomes and plasmids reconstructed. Chromosomes were defined as 'fully reconstructed' where a contig was >4Mb and fully circularised. Plasmids were either fully reconstructed, or partly reconstructed/misassembled. Fully reconstructed plasmids were defined as contigs that aligned to eachother (dnadiff from MUMmer4<sup>11</sup> v4.0.0), had the same replicon/Inc type assignment (PlasmidFinder<sup>12</sup> v2.1.6), had the same circulairity (either all circular, or all linear) and were within 200bp in length of eachother, from at least 2 assemblers. Misassembled plasmids were defined as fully reconstructed, but where one of the Inc typ, circularity, or length criteria were not met. 

Accuracy was assessed by substitutions and indels corrected by realigning Illumina short-reads to assemblies using Pypolca<sup>10</sup>, MLST assignment (mlst<sup>13</sup> v2.23.0), and the recovery of key resistance and virulence genes (AMRFinderPlus<sup>14</sup> v4.0.3 with the species flag inferred from Kraken2<sup>15</sup> v2.1.3).   
  
### Statistical analysis and visualisation
Statistical analysis and visualisation were done in R<sup>16</sup> v4.4.1. Proportions were compared with Fisher’s exact test, and counts with Wilcoxon signed-rank tests.


## Acknowledgements

This work was funded by the UK Health Security Agency (UKHSA) and supported by the National Institute for Health Research (NIHR) Health Protection Research Unit in Healthcare Associated Infections and Antimicrobial Resistance (NIHR200915), a partnership between the UKHSA and the University of Oxford, the NIHR Oxford Biomedical Research Centre (BRC) and the UKHSA [UKHSA PhD Funding Competition]. The cloud infrastructure used in this study was donated by Oracle Corporation.


## Citation

Nagy, D., Pennetta, V., Rodger, G., Hopkins, K., Jones, C., The NEKSUS Study Group, Hopkins, S., Crook, D., Walker, A.S., Robotham, J., Hopkins, K.L., Ledda, A., Williams, D., Hope, R., Brown, C.S., Stoesser, N., Lipworth, S. (2025). *Nanopore long-read only genome assembly is complete and accurate for surveillance of Enterobacterales bloodstream infections in England.* [Conference ePoster]. ESCMID Global 2025, Vienna, Austria. https://registration.escmid.org//AbstractList.aspx?e=30&header=0&preview=1&aig=-1&ai=33860 

## References
<sub>
1. Sanderson, N. D., Hopkins, K. M. V., Colpus, M., Parker, M., Lipworth, S., Crook, D., & Stoesser, N. (2024). Evaluation of the accuracy of bacterial genome reconstruction with Oxford Nanopore R10.4.1 long-read-only sequencing. Microb Genom, 10(5). https://doi.org/10.1099/mgen.0.001246
<br>2. Shen, W., Le, S., Li, Y., Hu, F. (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
<br>3. Wick, R.R. Autocycler. https://github.com/rrwick/Autocycler. 2025. doi:10.5281/zenodo.14642607 
<br>4. Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, https://doi.org/10.21105/joss.03941
<br>5. Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. (2019). Assembly of Long Error-Prone Reads Using Repeat Graphs. Nature Biotechnology. https://doi.org/doi:10.1038/s41587-019-0072-8 
<br>6. Bouras, G., Houtak, G., Wick, R. R., Mallawaarachchi, V., Roach, M. J., Papudeshi, B., Judd, L. M., Sheppard, A. E., Edwards, R. A., & Vreugde, S. (2024). Hybracter: Enabling Scalable, Automated, Complete and Accurate Bacterial Genome Assemblies. bioRxiv. https://doi.org/10.1101/2023.12.12.571215 
<br>7. Wick, R.R., Judd, L.M., Gorrie, C.L., Holt, K.E. (2017). Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLOS Computational Biology 13(6): e1005595. https://doi.org/10.1371/journal.pcbi.1005595 
<br>8. Oxford Nanopore Technologies (ONT). Medaka, Github. https://github.com/nanoporetech/medaka. [Accessed 20/11/2024]
<br>9. Wick, R.R., Holt, K.E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLOS Computational Biology 18(1): e1009802. https://doi.org/10.1371/journal.pcbi.1009802
<br>10. Bouras, G., Judd, L.M., Edwards, R.A., Vreugde, S., Stinear, T.P., Wick, R.R. (2024). How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics. doi: https://doi.org/10.1099/mgen.0.001254.
<br>11. Marçais, G, Delcher, AL, Phillippy, AM, Coston, R, Salzberg, SL, et al. (2018) MUMmer4: A fast and versatile genome alignment system. PLOS Computational Biology 14(1): e1005944. https://doi.org/10.1371/journal.pcbi.1005944
<br>12. Carattoli, A., Zankari, E., Garcia-Fernandez, A., Voldby Larsen, M., Lund, O., Villa, L., Moller Aarestrup, F., & Hasman, H. (2014). In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. 
<br>13. Seemann, T., mlst Github https://github.com/tseemann/mlst [Accessed 20/11/2024]
<br>14. Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J.G., Haendiges, J., Haft, D.H., Hoffmann, M., Pettengill, J.B., Prasad, A.B., Tillman, G.E., Tyson, G.H., Klimke, W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984.
<br>15. Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
<br>16. R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
</sub>
