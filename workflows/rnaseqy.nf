/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseqy_pipeline'

// Included for our pipeline
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main' 
include { UNZIPPER as UNZIP_FASTQ} from '../modules/local/unzipper/main'
include { UNZIPPER as UNZIP_GTF  } from '../modules/local/unzipper/main'
include { UNZIPPER as UNZIP_GFF  } from '../modules/local/unzipper/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main' 
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { CUSTOM_GETCHROMSIZES   } from '../modules/nf-core/custom/getchromsizes/main'  
include { SUBREAD_FEATURECOUNTS  } from '../modules/nf-core/subread/featurecounts/main' 
include { STRINGTIE_STRINGTIE    } from '../modules/nf-core/stringtie/stringtie/main'    
include { STRINGTIE_MERGE        } from '../modules/nf-core/stringtie/merge/main'  
include { MERGY                  } from '../modules/local/mergy/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow RNASEQY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_reference_fasta = Channel.of([ [ id: 'genome' ], file(params.fasta) ])
    ch_reference_gtf  = Channel.of([ [ id: 'genome' ], file(params.gtf) ])
    ch_reference_gff = Channel.of([ [ id: 'genome' ], file(params.gff) ])

    //
    // MODULE: Unzipper (selfwritten to unzip our files)
    //
    UNZIP_GTF(ch_reference_gtf)
    UNZIP_GFF(ch_reference_gff)
    // add Unzip to versions
    ch_versions = ch_versions.mix(UNZIP_GTF.out.versions).mix(UNZIP_GFF.out.versions)

    //
    // MODULE: TrimGalore (trim reads first)
    //
    TRIMGALORE(
        ch_samplesheet
    )
    // add TrimGalore reports to MultiQC and versions
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.map { it[1] })
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        TRIMGALORE.out.reads
    )
    // add FASTQC report to MultiQC and versions
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: STAR Genome Index Generation
    //
    STAR_GENOMEGENERATE(
        ch_reference_fasta,
        UNZIP_GTF.out.unzipped
    )

    ch_star_index = STAR_GENOMEGENERATE.out.index
    // add star_genomegenerate report to versions
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    
    //
    // MODULE: STAR Alignment
    //
    // for this we need to unzip the trimmed reads
    UNZIP_FASTQ(
        TRIMGALORE.out.reads
    )
    // add Unzip to versions
    ch_versions = ch_versions.mix(UNZIP_FASTQ.out.versions)

    STAR_ALIGN(
        UNZIP_FASTQ.out.unzipped,
        ch_star_index.collect(),
        UNZIP_GTF.out.unzipped.collect(),
        false,
        [],
        []
    )
    ch_bam = STAR_ALIGN.out.bam
    // add star_align reports to MultiQC and versions
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_out.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_progress.map { it[1] })
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // MODULE: Picard MarkDuplicates
    //
    // for this module we need two helper tools
    CUSTOM_GETCHROMSIZES (
        ch_reference_fasta
    )
    // add custom_getchromosomesize to versions
    ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    SAMTOOLS_SORT(
        ch_bam,
        ch_reference_fasta.collect(),
        []
    )
    // add samtools_sort to versions
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    
    PICARD_MARKDUPLICATES(
        SAMTOOLS_SORT.out.bam,
        ch_reference_fasta.collect(),
        CUSTOM_GETCHROMSIZES.out.fai.collect()
    )
    // add picard_markduplikates to versions
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    //
    // Start feature count
    //
    ch_gtf_path = UNZIP_GTF.out.unzipped.map { meta, gtf -> gtf }
    ch_gff_path = UNZIP_GFF.out.unzipped.map { meta, gff -> gff }
    ch_for_featurecounts = PICARD_MARKDUPLICATES.out.bam.combine(ch_gtf_path)
        .map { meta, bam, gtf ->
            tuple(meta, bam, gtf)}

    SUBREAD_FEATURECOUNTS(
        ch_for_featurecounts
    )
    // add subread_featurecounts to versions
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)


    STRINGTIE_STRINGTIE(
        PICARD_MARKDUPLICATES.out.bam,
        ch_gff_path.collect()
    )
    ch_transcript_gtfs = STRINGTIE_STRINGTIE.out.transcript_gtf.map { meta, gtf -> gtf }
    // add stringtie to versions
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

    STRINGTIE_MERGE(
        ch_transcript_gtfs.collect(), 
        ch_gff_path.collect()
    )
    // add stringtie_merge to versions
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)

    //
    // MODULE: Merge single count files and compute PCA of samples
    //
    // therfore two channels are necessary to run this module. 
    ch_stringtie = Channel.fromPath("${params.outdir}/stringtie")
    ch_outdir = Channel.fromPath("${params.outdir}")

    MERGY (
        ch_stringtie,
        ch_outdir,
        STRINGTIE_MERGE.out.gtf // this is a space holder to make sure that MERGY runs after all samples went through stringtie (discussed this with Mark as a workaround)
    )
    // add mergy to versions
    ch_versions = ch_versions.mix(MERGY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MERGY.out.table)
    ch_multiqc_files = ch_multiqc_files.mix(MERGY.out.plot)

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

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    // add MultiQC to versions
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    //
    // Save MultiQC to software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnaseqy_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    
    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
