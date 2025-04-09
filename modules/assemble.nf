process AUTOCYCLER_SUBSAMPLE {
    tag {sample}
    cpus 32
    publishDir "subsamples/autocycler_subsample/${sample}", pattern: "subsampled_reads", mode: "copy"
    publishDir "genome_sizes", pattern: "genome_size.txt", mode: 'copy', saveAs: { filename -> "${sample}_genome_size.txt"}

    input:
    tuple val(sample), val(source), path('reads.fastq.gz')

    output:
    tuple val(sample), val('subsampled_ont'), path("subsampled_reads"), emit: fol
    tuple val(sample), val('subsampled_ont'), path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), path("subsampled_reads/sample_02.fastq"), path("subsampled_reads/sample_03.fastq"), path("subsampled_reads/sample_04.fastq"), emit: fqs
    tuple val(sample), val('subsampled_ont'), val("01"), path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), emit: fq01
    tuple val(sample), val('subsampled_ont'), val("02"), path("genome_size.txt"), path("subsampled_reads/sample_02.fastq"), emit: fq02
    tuple val(sample), val('subsampled_ont'), val("03"), path("genome_size.txt"), path("subsampled_reads/sample_03.fastq"), emit: fq03
    tuple val(sample), val('subsampled_ont'), val("04"), path("genome_size.txt"), path("subsampled_reads/sample_04.fastq"), emit: fq04
    tuple val(sample), val('subsampled_ont'), path("genome_size.txt"), emit: txt

    script:
    
    """
    genome_size=\$(genome_size_raven.sh "reads.fastq.gz" $task.cpus)  # can set this manually if you know the value

    echo "\$genome_size" > genome_size.txt
    
    autocycler subsample --reads reads.fastq.gz --out_dir subsampled_reads --genome_size "\$genome_size"
    """


}

process SUBSAMPLE_ILLUMINA {

    tag {sample}

    label 'short'

    publishDir "subsamples/${task.process.replaceAll(":","_")}", mode: 'copy' //, saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    stageInMode 'copy'

    output:
    tuple val(sample), val('subsampled_illumina'), path("${sample}_r1.fastq.gz"), path("${sample}_r2.fastq.gz"), emit: fq

    script:
    """
    rasusa reads reads1.fastq.gz reads2.fastq.gz --output ${sample}_r1.fastq.gz ${sample}_r2.fastq.gz --coverage 100 --genome-size 5.6m --output-type g
    """
    stub:
    """
    touch ${sample}.fastq.gz
    """
}


process RAWQC_ONT{
    
    tag {sample}
    publishDir "raw_QC_CSVs/${source}", mode: 'copy'

    input:
    tuple val(sample), val(source), path('reads.fastq.gz')

    stageInMode 'copy'

    output:
    tuple val(sample), val(source), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all reads.fastq.gz > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
    }

process RAWQC_ILLUMINA {
    tag {sample}
    publishDir "raw_QC_CSVs/${source}", mode: 'copy'

    input:
    tuple val(sample), val(source), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    stageInMode 'copy'

    output:
    tuple val(sample), val(source), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all reads1.fastq.gz reads2.fastq.gz > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
    }


process RAWQC_ONT_SUBSAMPLED {

    tag {sample}
    publishDir "raw_QC_CSVs/${source}", mode: 'copy'

    input:
    tuple val(sample), val(source),  path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), path("subsampled_reads/sample_02.fastq"), path("subsampled_reads/sample_03.fastq"), path("subsampled_reads/sample_04.fastq")


    output:
    tuple val(sample), val(source), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all subsampled_reads/sample_01.fastq subsampled_reads/sample_02.fastq subsampled_reads/sample_03.fastq subsampled_reads/sample_04.fastq > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
    }




process AUTOCYCLER_ASSEMBLE {
    tag {sample}
    cpus 16

    publishDir "autocycler_assemblies/${sample}", pattern: "autocycler_assemblies", mode: "copy"
    publishDir "autocycler_assemblies/canu/${sample}", pattern: "autocycler_assemblies/canu_*.fasta", mode: "copy"
    publishDir "autocycler_assemblies/flye/${sample}", pattern: "autocycler_assemblies/flye_*.fasta", mode: "copy"
    publishDir "autocycler_assemblies/miniasm/${sample}", pattern: "autocycler_assemblies/miniasm_*.fasta", mode: "copy"
    publishDir "autocycler_assemblies/necat/${sample}", pattern: "autocycler_assemblies/necat_*.fasta", mode: "copy"
    publishDir "autocycler_assemblies/raven/${sample}", pattern: "autocycler_assemblies/raven_*.fasta", mode: "copy"
    
    input:
    tuple val(sample), val(source),  path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), path("subsampled_reads/sample_02.fastq"), path("subsampled_reads/sample_03.fastq"), path("subsampled_reads/sample_04.fastq")

    output:
    tuple val(sample), val('autocycler_assemblies'), path("genome_size.txt"), path("autocycler_assemblies"), emit: fol
    tuple val(sample), val('canu'), path("autocycler_assemblies/canu_01.fasta"), emit: canu
    tuple val(sample), val('flye'), path("autocycler_assemblies/flye_01.fasta"), emit: flye
    tuple val(sample), val('miniasm'), path("autocycler_assemblies/miniasm_01.fasta"), emit: miniasm
    tuple val(sample), val('necat'), path("autocycler_assemblies/necat_01.fasta"), emit: necat
    tuple val(sample), val('raven'), path("autocycler_assemblies/raven_01.fasta"), emit: raven

    script:
    """
    genome_size=\$(<genome_size.txt)

    mkdir -p autocycler_assemblies
    for assembler in canu flye miniasm necat raven; do
        for i in 01 02; do
            "\$assembler".sh subsampled_reads/sample_"\$i".fastq autocycler_assemblies/"\$assembler"_"\$i" "${task.cpus}" "\$genome_size"
        done
    done
    """
}

process AUTOCYCLER_ASSEMBLE_CANU {
    tag {sample + ' ' + number}
    cpus 8

    publishDir "autocycler_assemblies/canu/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/canu_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}





    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('canu'), val(number), path("canu_${number}.fasta"), emit: fa

    script:
    """
    genome_size=\$(<genome_size.txt)

    canu.sh subsampled_reads.fastq canu_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_FLYE {
    tag {sample + ' ' + number}
    cpus 8

    publishDir "autocycler_assemblies/flye/${sample}", mode: "copy"
    publishDir "assemblies/flye_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('flye'), val(number), path("flye_${number}.fasta"), emit: fa
    tuple val(sample), val('flye'), val(number), path("flye_${number}.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)

    flye.sh subsampled_reads.fastq flye_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_MINIASM {
    tag {sample + ' ' + number}
    cpus 8

    publishDir "autocycler_assemblies/miniasm/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/miniasm_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.fasta"), emit: fa
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.gfa"), optional: true, emit: gfa
    script:
    """
    genome_size=\$(<genome_size.txt)

    miniasm.sh subsampled_reads.fastq miniasm_${number} "${task.cpus}" "\$genome_size"
    """
}
process AUTOCYCLER_ASSEMBLE_NECAT {
    tag {sample + ' ' + number}
    cpus 8

    publishDir "autocycler_assemblies/necat/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/necat_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/necat_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/necat_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/necat_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('necat'), val(number), path("necat_${number}.fasta"), emit: fa

  script:
    """
    genome_size=\$(<genome_size.txt)

    necat.sh subsampled_reads.fastq necat_${number} "${task.cpus}" "\$genome_size"
    """
}
process AUTOCYCLER_ASSEMBLE_RAVEN {
    tag {sample + ' ' + number}
    cpus 8

    publishDir "autocycler_assemblies/raven/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/raven_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('raven'), val(number), path("raven_${number}.fasta"), emit: fa
    tuple val(sample), val('raven'), val(number), path("raven_${number}.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)

    raven.sh subsampled_reads.fastq raven_${number} "${task.cpus}" "\$genome_size"
    """
}



process HYBRACTER_LONG_INDIVIDUAL {
    scratch true
    errorStrategy 'retry'
    cpus 1
    tag {sample + ' ' + number}

    publishDir "assemblies/hybracter_long/${sample}/${number}", mode: 'copy'
 
    input:
    tuple val(sample), val(source), val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_${number}_final.fasta"), emit: fasta
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out"), emit: fol
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_${number}_plasmid.fasta"), optional: true, emit: plasmids_candidate
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/incomplete/${sample}_${number}_final.fasta"), optional: true, emit: plasmids_fallback
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_${number}_chromosome.fasta"), optional: true, emit: chromosomes
 
    script:
    """
    hybracter long-single -l subsampled_reads.fastq -s ${sample}_${number} -c 4000000 -o hybracter_long_out -t ${task.cpus} --databases "${projectDir}/bin/hybracter_databases" #--medakaModel ${params.inputmodel}
    """
}


process AUTOCYCLER_COMPRESS { 
    tag {sample} 
    cpus 4 
    publishDir "assemblies/autocycler_consensus/${sample}", mode: "copy" 
    publishDir "assemblies/autocycler_consensus", mode: "copy", pattern: "autocycler_out/consensus_assembly.fasta", saveAs: { filename -> "${sample}.fasta" } 
    publishDir "assemblies/autocycler_consensus", mode: "copy", pattern: "autocycler_out/consensus_assembly.gfa", saveAs: { filename -> "${sample}.gfa" }
  
  
  
    input:
    tuple val(sample), path('assembly_01.fasta'), path('assembly_02.fasta'), path('assembly_03.fasta'), path('assembly_04.fasta'), path('assembly_05.fasta'), path('assembly_06.fasta'), path('assembly_07.fasta'), path('assembly_08.fasta'), path('assembly_09.fasta'), path('assembly_10.fasta'), path('assembly_11.fasta'), path('assembly_12.fasta'), path('assembly_13.fasta'), path('assembly_14.fasta'), path('assembly_15.fasta'), path('assembly_16.fasta'), path('assembly_17.fasta'), path('assembly_18.fasta'), path('assembly_19.fasta'), path('assembly_20.fasta')

    output:
    tuple val(sample), val("autocycler_consensus"), path("autocycler_out/consensus_assembly.fasta"), emit: fa
    tuple val(sample), val("autocycler_consensus"), path("autocycler_out/consensus_assembly.gfa"), emit: gfa
    tuple val(sample), val("autocycler_consensus"), path("autocycler_out"), emit: fol

    script:
    """
 
    mkdir -p "autocycler_assemblies"
    
    for assembly_file in assembly_01.fasta assembly_02.fasta assembly_03.fasta assembly_04.fasta assembly_05.fasta assembly_06.fasta assembly_07.fasta assembly_08.fasta assembly_09.fasta assembly_10.fasta assembly_11.fasta assembly_12.fasta assembly_13.fasta assembly_14.fasta assembly_15.fasta assembly_16.fasta assembly_17.fasta assembly_18.fasta assembly_19.fasta assembly_20.fasta; do
       echo "Processing \${assembly_file}"
       cp "\${assembly_file}" "autocycler_assemblies/\${assembly_file}"
    done

    
    # Step 3: compress the input assemblies into a unitig graph
    autocycler compress -i autocycler_assemblies -a autocycler_out --max_contigs 80

    # Step 4: cluster the input contigs into putative genomic sequences
    autocycler cluster -a autocycler_out --max_contigs 80 --min_assemblies 5

    # Steps 5 and 6: trim and resolve each QC-pass cluster
    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c"
        autocycler resolve -c "\$c"
    done

    # Step 7: combine resolved clusters into a final assembly
    autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa


    """
}


process SUBSAMPLE_ONT {
    
    tag {sample}

    label 'short'

    publishDir "subsamples/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fastq.gz"}

    input:
    tuple val(sample), val(source), path('reads.fastq.gz')

    stageInMode 'copy'    

    output:
    tuple val(sample), val('ont_subsampled'), path("${sample}.fastq.gz"), emit: fq

    script:
    """
    rasusa reads reads.fastq.gz --coverage 40 --genome-size 5.6m | gzip > ${sample}.fastq.gz
    ##rasusa reads reads.fastq.gz --output ${sample}.fastq.gz --coverage 40 --genome-size 5.6m --output-type g 
    """ 
    stub:
    """
    touch ${sample}.fastq.gz
    """
}


    
process FLYE {
    
    tag {sample}
    cpus 4

    label 'short'

    publishDir "assemblies/${task.process.replaceAll(":","_")}", pattern: "assembly/assembly.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/${task.process.replaceAll(":","_")}", pattern: "assembly/assembly_graph.gfa", mode: 'copy', saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/${task.process.replaceAll(":","_")}", pattern: "assembly/assembly_info.txt", mode: 'copy', saveAs: { filename -> "${sample}.txt"}

    input:
    
    tuple val(sample), val(source),  file('reads.fastq.gz')

    output:
    tuple val(sample), val('flye'), path('assembly/assembly.fasta'), emit: fasta
    tuple val(sample), val('flye'), path('assembly/assembly_graph.gfa'), emit: gfa
    tuple val(sample), val('flye'), path('assembly/assembly_info.txt'), emit: txt


    script:
    """
    flye -o assembly --meta --threads ${task.cpus} --nano-hq reads.fastq.gz
    """
    stub:
    """
    mkdir assembly
    touch assembly/assembly.fasta
    """
}

process MEDAKA {
    
    label 'medaka'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "assemblies/${source_contigs}_medaka", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source_contigs), path("contigs.fasta"), val(source_reads), val("01"), path("genome_size.txt"), path("reads.fastq")
    
    output:
    tuple val(sample), val("${source_contigs}_medaka"), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus -i reads.fastq -d contigs.fasta -o output -t ${task.cpus}  --bacteria -m r1041_e82_400bps_bacterial_methylation
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}

process MEDAKA_FULL {

    label 'medaka'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "assemblies/${source_contigs}_medaka_full", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source_contigs), path('contigs.fasta'), val(source_reads), path('reads.fastq.gz')
  
    output:
    tuple val(sample), val("${source_contigs}_medaka_full"), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus -i reads.fastq.gz -d contigs.fasta -o output -t ${task.cpus}  --bacteria -m r1041_e82_400bps_bacterial_methylation
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}

process POLYPOLISH {
    label 'polypolish'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "assemblies/${source_contigs}_polypolish", pattern: "polished.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source_contigs), path('contigs.fasta'), val(source_illumina), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    output:
    tuple val(sample), val(source_contigs), val('polypolish'), path('polished.fasta'), path('reads1.fastq.gz'), path('reads2.fastq.gz'), emit: fasta

    script:
    """
    bwa index contigs.fasta
    bwa mem -t ${task.cpus} -a contigs.fasta reads1.fastq.gz > alignments_1.sam
    bwa mem -t ${task.cpus} -a contigs.fasta reads2.fastq.gz > alignments_2.sam
    polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish polish contigs.fasta filtered_1.sam filtered_2.sam > polished.fasta
    rm *.amb *.ann *.bwt *.pac *.sa *.sam
    """
    stub:
    """
    touch polished.fasta
    """
}


process PYPOLCA {
    label 'pypolca'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "assemblies/${source_contigs}_polypolish_pypolca", pattern: "pypolca_out/pypolca_corrected.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 
    publishDir "assemblies/${source_contigs}_polypolish_pypolca/${sample}", mode: 'copy' 

    input: 
    tuple val(sample), val(source_contigs), val('polypolish'), path('polypolish_polished.fasta'), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    output:
    tuple val(sample), val("${source_contigs}_polypolish_pypolca"), path('pypolca_out/pypolca_corrected.fasta'), emit: fasta
    tuple val(sample), val("${source_contigs}_polypolish_pypolca"), path('pypolca_out'), emit: fol

    script:
    """
    pypolca run -a polypolish_polished.fasta -1 reads1.fastq.gz -2 reads2.fastq.gz -t ${task.cpus} -o pypolca_out --careful
   
    """
}



process UNICYCLER {
    label 'unicycler'
    tag {sample}
    cpus 4 
    
    publishDir "assemblies/unicycler", pattern: "output/assembly.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/unicycler", pattern: "output/assembly.gfa", mode: 'copy', saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(subsampled_ont), val(number_01), path("genome_size.txt"), path('reads.fastq.gz'), val(subsampled_illumina), path('reads1.fastq.gz'), path('reads2.fastq.gz')


    output:
    tuple val(sample), val('unicycler'), path('output/assembly.fasta'), emit: fasta
    tuple val(sample), val('unicycler'), path('output/assembly.gfa'), emit: gfa

    script:
    """
    unicycler -1 reads1.fastq.gz -2 reads2.fastq.gz -l reads.fastq.gz -o output -t ${task.cpus} --mode normal
    """
    stub:
    """
    mkdir output
    touch output/assembly.fasta
    touch output/assembly.gfa
    """

}

process UNICYCLER_BOLD {
    label 'unicycler'
    tag {sample}
    cpus 4
    label 'short'

    publishDir "assemblies/unicycler_bold", pattern: "output/assembly.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/unicycler_bold", pattern: "output/assembly.gfa", mode: 'copy', saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(subsampled_ont), val(number_01), path("genome_size.txt"), path('reads.fastq.gz'), val(subsampled_illumina), path('reads1.fastq.gz'), path('reads2.fastq.gz')


    output:
    tuple val(sample), val('unicycler_bold'), path('output/assembly.fasta'), emit: fasta
    tuple val(sample), val('unicycler_bold'), path('output/assembly.gfa'), emit: gfa

    script:
    """
    unicycler -1 reads1.fastq.gz -2 reads2.fastq.gz -l reads.fastq.gz -o output -t ${task.cpus} --mode bold
    """
    stub:
    """
    mkdir output
    touch output/assembly.fasta
    """
}



process HYBRACTER_HYBRID_INDIVIDUAL{
    tag {sample}
    scratch true
    label 'online'
    cpus 4

    publishDir "assemblies/hybracter_hybrid/${sample}", mode: 'copy'
 
    input:
    tuple val(sample), val(subsampled_ont), val(number_01), path("genome_size.txt"), path('reads.fastq.gz'), val(subsampled_illumina), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    output:
    tuple val(sample), val('hybracter_hybrid'), path("hybracter_hybrid_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_final.fasta"), emit: fasta 
    tuple val(sample), val('hybracter_hybrid'), path("hybracter_hybrid_out"), emit: fol
 
    script:
    """
    hybracter hybrid-single -l reads.fastq.gz -1 reads1.fastq.gz -2 reads2.fastq.gz -s ${sample} --auto -o hybracter_hybrid_out -t ${task.cpus}  --databases "${projectDir}/bin/hybracter_databases" #--medakaModel ${params.inputmodel} 
    """
}


process PYPOLCA_ANALYSE {
    label 'pypolca'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "pypolca/${source_contigs}", pattern: "pypolca_out/pypolca_corrected.fasta", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}
    publishDir "pypolca/${source_contigs}", pattern: "pypolca_out/pypolca.report", mode: 'copy', saveAs: { filename -> "${sample}.report"}
    publishDir "pypolca/${source_contigs}", pattern: "pypolca_out/pypolca.vcf", mode: 'copy', saveAs: { filename -> "${sample}.vcf"}
    publishDir "pypolca/${source_contigs}/${sample}", mode: 'copy'

    input:
    tuple val(sample), val(source_contigs), path('polished_assembly.fasta'), val(source_illumina), path('reads1.fastq.gz'), path('reads2.fastq.gz')

    output:
    tuple val(sample), val("${source_contigs}_pypolca"), path('pypolca_out/pypolca_corrected.fasta'), emit: fasta
    tuple val(sample), val("${source_contigs}_pypolca"), path('pypolca_out/pypolca.report'), emit: report
    tuple val(sample), val("${source_contigs}_pypolca"), path('pypolca_out/pypolca.vcf'), emit: vcf
    tuple val(sample), val("${source_contigs}_pypolca"), path('pypolca_out'), emit: fol

    script:
    """
    pypolca run -a polished_assembly.fasta -1 reads1.fastq.gz -2 reads2.fastq.gz -t ${task.cpus} -o pypolca_out --careful

    """
}


process BAKTA_DOWNLOAD {
    label 'online'

    output:
    path('db')

    script:
    """
    bakta_db download --output db --type full
    """

}

process BAKTA {
    tag {sample + ' ' + assembler}
    errorStrategy 'ignore'
    
    publishDir "bakta/${assembler}", mode: 'copy'
    publishDir "bakta_gff/${assembler}", pattern: "*.gff3", mode: 'copy', saveAs: { filename -> "${sample}.gff3"}
    
    cpus 16

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')
    
        
    output:
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.gff3"), emit: gff3
    tuple val(sample), val(assembler), path("${sample}_bakta"), emit: fol


    script:
    """
    bakta --db ${projectDir}/bin/db/db -t ${task.cpus} --prefix ${sample} -o ${sample}_bakta/ assembly.fasta
    """
    stub:
    """
    mkdir ${sample}_bakta
    touch ${sample}_bakta/${sample}.gbk 
    """
}

process SEQKIT {
    tag {sample + ' ' + assembler}

    publishDir "assembly_QC_CSVs/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all assembly.fasta > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
}

process MLST {
    tag {sample + ' ' + assembler}
    cpus 4

    label 'short'

    publishDir "MLST/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path("assembly.fasta")
  
    output:
    path("${sample}.tsv"), optional: true

    script:
    """
    mlst assembly.fasta > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
}

process DNADIFF {
    tag {sample + ' ' + assembler}

    publishDir "dnadiff/${assembler}", mode: 'copy'
    publishDir "snps/${assembler}", pattern: '*.snps', mode: 'copy', saveAs: { filename -> "${sample}.snps"}

    input:
    tuple val(sample), val(ref_assembler),  path('ref.fa'), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}.report"), path("${sample}.delta"), path("${sample}.1coords"), emit: all
    tuple val(sample), val(assembler), path("${sample}.1coords"), emit: coords
    tuple val(sample), val(assembler), path("${sample}.snps"), emit: snps


    script:
    """
    dnadiff ref.fa assembly.fasta -p ${sample}
    """
    stub:
    """
    touch ${sample}.report ${sample}.delta ${sample}.1coords
    """
}

process DNADIFF_REV {
    tag {sample + ' ' + assembler}

    publishDir "dnadiff_rev/${assembler}", mode: 'copy'
    publishDir "snps_rev/${assembler}", pattern: '*.snps', mode: 'copy', saveAs: { filename -> "${sample}.snps"}

    input:
    tuple val(sample), val(ref_assembler),  path('ref.fa'), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}.report"), path("${sample}.delta"), path("${sample}.1coords"), emit: all
    tuple val(sample), val(assembler), path("${sample}.1coords"), emit: coords
    tuple val(sample), val(assembler), path("${sample}.snps"), emit: snps


    script:
    """
    dnadiff  assembly.fasta ref.fa -p ${sample}
    """
    stub:
    """
    touch ${sample}.report ${sample}.delta ${sample}.1coords
    """
}


process KRAKEN2 {
    scratch true
    tag {sample + ' ' + assembler}
    cpus 32

    publishDir "kraken2/${assembler}", pattern: "*_k2_report.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_report.tsv"}
    publishDir "kraken2/${assembler}", pattern: "*_k2_output.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_output.tsv"}

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}_k2_output.tsv"),  path('species.txt'), emit: output
    tuple val(sample), val(assembler), path("${sample}_k2_report.tsv"),  path('species.txt'), emit: report
    tuple val(sample), val(assembler), path('assembly.fasta'), path("${sample}_k2_output.tsv"), path("${sample}_k2_report.tsv"), emit: fasta

    script:
    """
    kraken2 --db /mnt/nfs-gram-negative-study/assembly/flye-medaka-pipeline_v2/kraken2_db assembly.fasta  --use-names --report ${sample}_k2_report.tsv --output ${sample}_k2_output.tsv  --threads ${task.cpus}

    species=\$(awk '\$4 == "S" {print substr(\$0, index(\$0,\$6)) ; exit}' ${sample}_k2_report.tsv)

    if [[ -z "\$species" ]]; then
        echo "Error: No matching row with 4th column 'S' found in output.tsv" >&2
        exit 1
    fi

    echo \$species > species.txt
    """

    stub:
    """
    touch ${sample}_k2_report.tsv
    touch ${sample}_k2_output.tsv
    """
}


process SPECIES_FROM_KRAKEN2 {
    tag {sample + ' ' + assembler}
    publishDir "kraken2/${assembler}", pattern: "*.txt", mode: 'copy', saveAs: { filename -> "${sample}_chromosome_species.txt"}
    

    input:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('output.tsv'), path('report.tsv')

    output:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('main_chromosome_species.txt')

    script:
    """
    awk -F'\t' '\$4 > 1000000' output.tsv | sort -t\$'\t' -k4,4nr > filtered.tsv

    if [ -s filtered.tsv ]; then
         awk -F'\t' '{print \$3}' filtered.tsv | awk '{print \$1, \$2}' | paste -sd, - > main_chromosome_species.txt

    else
        awk '\$4 == "S"' report.tsv > filtered.tsv
        max_value=\$(awk '{print \$1}' filtered.tsv | sort -nr | head -n1)
        awk -v max="\$max_value" '\$1 == max {print substr(\$0, index(\$0,\$6))}' filtered.tsv | paste -sd, - > main_chromosome_species.txt
    fi
    """

}

process AMRFINDERPLUS {
    errorStrategy 'ignore'    
    tag {sample + ' ' + assembler}
    cpus 32

    label 'short'

    publishDir "AMRFinderPlus/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('species.txt')

    output:
    tuple val(sample), val(assembler), path("${sample}.tsv"), path('species.txt')

    script:
    """
    species=\$(<species.txt)

    if [[ \$species == Escherichia* ]]; then
        amrfinder -n assembly.fasta -O Escherichia -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1

    elif [[ \$species == "Klebsiella pneumoniae" ]]; then
        amrfinder -n assembly.fasta -O Klebsiella_pneumoniae -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Klebsiella oxytoca" ]]; then
        amrfinder -n assembly.fasta -O Klebsiella_oxytoca -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Enterobacter cloacae" ]]; then
        amrfinder -n assembly.fasta -O Enterobacter_cloacae -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Enterobacter asburiae" ]]; then
        amrfinder -n assembly.fasta -O Enterobacter_absuriae -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Citrobacter freundii" ]]; then
        amrfinder -n assembly.fasta -O Citrobacter_freundii -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Serratia marcescens" ]]; then
        amrfinder -n assembly.fasta -O Serratia_marcescens -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    else
        amrfinder -n assembly.fasta -o ${sample}.tsv --plus -d ${projectDir}/bin/amrfinder_db/2024-12-18.1
    fi
    """
}



process PLASMIDFINDER {
    tag {sample + ' ' + assembler}
    cpus 4

    label 'short'

    publishDir "PlasmidFinder/${assembler}", mode: 'copy', pattern: "*.txt"
    publishDir "PlasmidFinder/${assembler}", mode: 'copy', pattern: "plasmidfinder_out"

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')
    
    output:
    tuple val(sample), val(assembler), path("plasmidfinder_out"), emit: fol
    tuple val(sample), val(assembler), path("${sample}_plasmidfinder.txt"), emit: txt

    script:
    """
    mkdir -p plasmidfinder_out 
    plasmidfinder.py -i assembly.fasta -o plasmidfinder_out -p ${projectDir}/bin/plasmidfinder_db/database > ${sample}_plasmidfinder.txt 
    """
}