#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// include modules
include {RAWQC_ONT} from './modules/assemble.nf'
include {RAWQC_ILLUMINA} from './modules/assemble.nf'
include {RAWQC_ONT_SUBSAMPLED} from './modules/assemble.nf'
include {AUTOCYCLER_SUBSAMPLE} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_CANU} from './modules/assemble.nf' 
include {AUTOCYCLER_ASSEMBLE_FLYE} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_MINIASM} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_NECAT} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_RAVEN} from './modules/assemble.nf'
include {AUTOCYCLER_COMPRESS} from './modules/assemble.nf'
include {SUBSAMPLE_ONT} from './modules/assemble.nf'
include {SUBSAMPLE_ILLUMINA} from './modules/assemble.nf'
include {FLYE as FLYE} from './modules/assemble.nf'
include {MEDAKA} from './modules/assemble.nf'
include {MEDAKA_FULL} from './modules/assemble.nf'
include {POLYPOLISH} from './modules/assemble.nf'
include {PYPOLCA} from './modules/assemble.nf'
include {UNICYCLER} from './modules/assemble.nf'
include {UNICYCLER_BOLD} from './modules/assemble.nf'
include {HYBRACTER_HYBRID_INDIVIDUAL} from './modules/assemble.nf'
include {HYBRACTER_LONG_INDIVIDUAL} from './modules/assemble.nf'
include {PYPOLCA_ANALYSE} from './modules/assemble.nf'
include {SEQKIT} from './modules/assemble.nf'
include {BAKTA_DOWNLOAD} from './modules/assemble.nf'
include {BAKTA} from './modules/assemble.nf'
include {MLST} from './modules/assemble.nf'
include {DNADIFF} from './modules/assemble.nf'
include {DNADIFF_REV} from './modules/assemble.nf'
include {KRAKEN2} from './modules/assemble.nf'
include {SPECIES_FROM_KRAKEN2} from './modules/assemble.nf'
include {AMRFINDERPLUS} from './modules/assemble.nf'
include {PLASMIDFINDER} from './modules/assemble.nf'


workflow {
    
    // INPUT CHANNELS: 
    Channel.fromPath( "${params.inputFastq_ONT}/*.fastq" )
           .map{ file -> tuple(file.simpleName, 'raw_ont', file) }
           .set{ ont_labeled_ch }
    //ont_labeled_ch.view()
     
         
    Channel.fromFilePairs("${params.inputFastq_illumina}/*{_R1,_R2}_001.fastq.gz")
           .map{ row -> tuple(row[0], 'raw_illumina', row[1][0], row[1][1])}
           .set{illumina_labeled_ch}
    //illumina_labeled_ch.view()    

    main:
    
    AUTOCYCLER_SUBSAMPLE(ont_labeled_ch)
    //AUTOCYCLER_SUBSAMPLE.out.fq.view() 
    //AUTOCYCLER_SUBSAMPLE.out.fol.view()
    SUBSAMPLE_ILLUMINA(illumina_labeled_ch)
   
    RAWQC_ONT_SUBSAMPLED(AUTOCYCLER_SUBSAMPLE.out.fqs)
    RAWQC_ONT(ont_labeled_ch)
    RAWQC_ILLUMINA(illumina_labeled_ch.mix(SUBSAMPLE_ILLUMINA.out.fq))
    
    //AUTOCYCLER_ASSEMBLE.out.fol.view()
    

    subsampled_fqs_mixed = AUTOCYCLER_SUBSAMPLE.out.fq01.mix(AUTOCYCLER_SUBSAMPLE.out.fq02, AUTOCYCLER_SUBSAMPLE.out.fq03, AUTOCYCLER_SUBSAMPLE.out.fq04)

    AUTOCYCLER_ASSEMBLE_CANU(subsampled_fqs_mixed)
    //AUTOCYCLER_ASSEMBLE_CANU.out.fa.view()
    AUTOCYCLER_ASSEMBLE_FLYE(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_MINIASM(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_RAVEN(subsampled_fqs_mixed)
    HYBRACTER_LONG_INDIVIDUAL(subsampled_fqs_mixed)
    
    hybracter_long_01_ch = HYBRACTER_LONG_INDIVIDUAL.out.fasta
          .filter { sample, assembler, number, fasta_path -> number == '01'}
          .map { sample, assembler, number, fasta_path -> tuple(sample, assembler, fasta_path) }
      

    
    HYBRACTER_HYBRID_INDIVIDUAL(AUTOCYCLER_SUBSAMPLE.out.fq01.combine(SUBSAMPLE_ILLUMINA.out.fq, by:0))
    
    UNICYCLER(AUTOCYCLER_SUBSAMPLE.out.fq01.combine(SUBSAMPLE_ILLUMINA.out.fq, by:0))
    UNICYCLER_BOLD(AUTOCYCLER_SUBSAMPLE.out.fq01.combine(SUBSAMPLE_ILLUMINA.out.fq, by:0))
    

    assemblies_for_autocycler =  AUTOCYCLER_ASSEMBLE_CANU.out.mix(AUTOCYCLER_ASSEMBLE_FLYE.out.fa, AUTOCYCLER_ASSEMBLE_MINIASM.out.fa, AUTOCYCLER_ASSEMBLE_RAVEN.out.fa, HYBRACTER_LONG_INDIVIDUAL.out.fasta)
    assemblies_for_autocycler_grouped = assemblies_for_autocycler.groupTuple(sort: true)
    assemblies_for_autocycler_grouped = assemblies_for_autocycler_grouped.map {sample, assemblers, numbers, paths ->
             return tuple(sample, paths)
             }


    //assemblies_for_autocycler_grouped.view()
    flattened_assemblies = assemblies_for_autocycler_grouped.map {sample, paths -> tuple(sample, *paths) }
    //flattened_assemblies.view()

   AUTOCYCLER_COMPRESS(flattened_assemblies)
   //AUTOCYCLER_COMPRESS.out.fa.view()
   
    flye_01_ch = AUTOCYCLER_ASSEMBLE_FLYE.out.fa
          .filter { sample, assembler, number, fasta_path -> number == '01'}
          .map { sample, assembler, number, fasta_path -> tuple(sample, assembler, fasta_path) }
    //flye_01_ch.view()
    
    long_read_assemblies_for_polishing_ch = flye_01_ch.mix(AUTOCYCLER_COMPRESS.out.fa)

    polishing_input_long_read_ch =  long_read_assemblies_for_polishing_ch.combine(AUTOCYCLER_SUBSAMPLE.out.fq01, by: 0)
    MEDAKA(polishing_input_long_read_ch)
    flye_medaka_ch = MEDAKA.out.fasta
           .filter {sample, source, assembly_path -> source == 'flye_medaka'}
    //flye_medaka_ch.view()
    autocycler_medaka_ch = MEDAKA.out.fasta
           .filter {sample, source, assembly_path -> source == 'autocycler_consensus_medaka'}
    //autocycler_medaka_ch.view()
    
    medaka_full_polishing_input_long_read_ch =  long_read_assemblies_for_polishing_ch.combine(ont_labeled_ch, by: 0)
    MEDAKA_FULL(medaka_full_polishing_input_long_read_ch)
   //MEDAKA_FULL.out.view()
    flye_medaka_full_ch = MEDAKA_FULL.out.fasta
           .filter {sample, source, assembly_path -> source == 'flye_medaka_full'}
    //flye_medaka_full_ch.view()
    autocycler_medaka_full_ch = MEDAKA_FULL.out.fasta
           .filter {sample, source, assembly_path -> source == 'autocycler_consensus_medaka_full'}
    //autocycler_medaka_full_ch.view()

    polishing_input_short_read_ch =  long_read_assemblies_for_polishing_ch.combine(illumina_labeled_ch, by: 0)
    POLYPOLISH(polishing_input_short_read_ch)
    //POLYPOLISH.out.view()
    PYPOLCA(POLYPOLISH.out)
    //PYPOLCA.out.fasta.view()
     flye_pypolca_ch = PYPOLCA.out.fasta
           .filter {sample, source, assembly_path -> source == 'flye_polypolish_pypolca'}
    //flye_pypolca_ch.view()
    autocycler_pypolca_ch = PYPOLCA.out.fasta
           .filter {sample, source, assembly_path -> source == 'autocycler_consensus_polypolish_pypolca'}
    //autocycler_pypolca_ch.view()  
    
    assemblies_to_compare = flye_01_ch.mix(flye_medaka_ch, flye_medaka_full_ch, flye_pypolca_ch, AUTOCYCLER_COMPRESS.out.fa, autocycler_medaka_ch, autocycler_medaka_full_ch, UNICYCLER.out.fasta, UNICYCLER_BOLD.out.fasta, hybracter_long_01_ch, HYBRACTER_HYBRID_INDIVIDUAL.out.fasta) 
    //assemblies_to_compare.view()

    all_assemblies = assemblies_to_compare.mix(autocycler_pypolca_ch)
    //all_assemblies.view()
    SEQKIT(all_assemblies)
    //SEQKIT.out.view()
 
    MLST(all_assemblies)
    //BAKTA_DOWNLOAD()

    
    //all_assemblies.combine(BAKTA_DOWNLOAD.out).view()
    BAKTA(all_assemblies)
            
    DNADIFF(autocycler_pypolca_ch.combine(assemblies_to_compare, by: 0))
    //DNADIFF_REV(autocycler_pypolca_ch.combine(assemblies_to_compare, by: 0))
    


    pypolca_analyse_input_ch = all_assemblies.combine(illumina_labeled_ch, by: 0)
    PYPOLCA_ANALYSE(pypolca_analyse_input_ch)


    KRAKEN2(all_assemblies)
    SPECIES_FROM_KRAKEN2(KRAKEN2.out.fasta)
    AMRFINDERPLUS(SPECIES_FROM_KRAKEN2.out)

    PLASMIDFINDER(all_assemblies)
}

