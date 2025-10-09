/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main' 
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main' 
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { UNZIPPER               } from '../modules/local/unzipper/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { CUSTOM_GETCHROMSIZES   } from '../modules/nf-core/custom/getchromsizes/main'
include { SALMON_INDEX           } from '../modules/nf-core/salmon/index/main' 
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'  
include { CUSTOM_TX2GENE         } from '../modules/nf-core/custom/tx2gene/main' 
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseqy_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow RNASEQY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    //ch_samplesheet.dump()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // MODULE: TrimGalore (trim reads first)
    //
    TRIMGALORE_OUT = TRIMGALORE(ch_samplesheet)
    // add TrimGalore reports to MultiQC and versions
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE_OUT.zip.map { it[1] })
    ch_versions = ch_versions.mix(TRIMGALORE_OUT.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        TRIMGALORE_OUT.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnaseqy_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    

   //
   // MODULE: STAR Genome Index Generation
    //
    ch_reference_fasta = Channel.of([ [ id: 'genome' ], file(params.fasta_genome) ])
    ch_annotation_gtf  = Channel.of([ [ id: 'genome' ], file(params.gtf) ])

    STAR_GENOMEGENERATE_OUT = STAR_GENOMEGENERATE(
        ch_reference_fasta,
        ch_annotation_gtf
    )

    ch_star_index = STAR_GENOMEGENERATE_OUT.index
    

    //
    // MODULE: STAR Alignment
    //
    UNZIPPER_OUT = UNZIPPER(TRIMGALORE_OUT.reads)


    STAR_ALIGN_OUT = STAR_ALIGN(
        UNZIPPER_OUT.reads,
        ch_star_index.collect(),
        ch_annotation_gtf.collect(),
        false,
        [],
        []
    )

    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_OUT.log_final.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_OUT.log_out.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_OUT.log_progress.map { it[1] })
    ch_versions = ch_versions.mix(STAR_ALIGN_OUT.versions.first())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    ch_bam = STAR_ALIGN_OUT.bam

    CUSTOM_GETCHROMSIZES_OUT = CUSTOM_GETCHROMSIZES (
        ch_reference_fasta
    )

    SAMTOOLS_SORT_OUT = SAMTOOLS_SORT(
        ch_bam,
        ch_reference_fasta.collect(),
        []
    )

    
    PICARD_MARKDUPLICATES_OUT = PICARD_MARKDUPLICATES(
        SAMTOOLS_SORT_OUT.bam,
        ch_reference_fasta.collect(),
        CUSTOM_GETCHROMSIZES_OUT.fai.collect()
    )

    //
    // MODULE: Salmon Index
    //
    SALMON_INDEX_OUT = SALMON_INDEX(
        Channel.of(file(params.fasta_genome)),
        Channel.of(file(params.fasta_transcripts))
    )

    //
    // MODULE: Salmon Quant
    //
    SALMON_QUANT_OUT = SALMON_QUANT(
        TRIMGALORE_OUT.reads,
        SALMON_INDEX_OUT.index,         // <- output from SALMON_INDEX
        Channel.of(file(params.gtf)),
        Channel.of(file(params.fasta_transcripts)),
        false,
        []
    )

    /*CUSTOM_TX2GENE_OUT = CUSTOM_TX2GENE(
        Channel.of(file(params.gtf)),
        [],
        [],
        [],
        []
    )*/



    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
