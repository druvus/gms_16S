#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed locally
//

include { MERGE_BARCODES              } from '../modules/local/merge_barcodes/main'
include { MERGE_BARCODES_SAMPLESHEET  } from '../modules/local/merge_barcodes_samplesheet/main'
include { GENERATE_INPUT              } from '../modules/local/generate_input/main'
include { EMU_ABUNDANCE               } from '../modules/local/emu/abundance/main'

//
// MODULE: Installed directly from nf-core/modules
//

include { NANOPLOT as NANOPLOT1       } from '../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT2       } from '../modules/nf-core/nanoplot/main'
include { PORECHOP_ABI                } from '../modules/nf-core/porechop/abi/main'
include { FILTLONG                    } from '../modules/nf-core/filtlong/main'
include { KRONA_KTIMPORTTAXONOMY      } from '../modules/nf-core/krona/ktimporttaxonomy/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'

include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow GMSEMU {
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Handle input based on parameters
    if (params.merge_fastq_pass && !params.barcodes_samplesheet) {
        MERGE_BARCODES(params.merge_fastq_pass)
        GENERATE_INPUT(MERGE_BARCODES.out.fastq_dir_merged)
        ch_input = GENERATE_INPUT.out.sample_sheet_merged
    } else if (params.merge_fastq_pass && params.barcodes_samplesheet) {
        MERGE_BARCODES_SAMPLESHEET(params.barcodes_samplesheet, params.merge_fastq_pass)
        GENERATE_INPUT(MERGE_BARCODES_SAMPLESHEET.out.fastq_dir_merged)
        ch_input = GENERATE_INPUT.out.sample_sheet_merged
    } else {
        ch_input = file(params.input)
    }

    // Input check
    INPUT_CHECK(ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    // QC
    NANOPLOT1(INPUT_CHECK.out.reads)
    ch_versions = ch_versions.mix(NANOPLOT1.out.versions.first())

    FASTQC(INPUT_CHECK.out.reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())



    // Processing
    if (params.adapter_trimming && !params.quality_filtering) {
        PORECHOP_ABI(INPUT_CHECK.out.reads)
        ch_processed_reads = PORECHOP_ABI.out.reads.map { meta, reads -> [ meta + [single_end: true], reads ] }
        ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_ABI.out.log)
    } else if (!params.adapter_trimming && params.quality_filtering) {
        ch_processed_reads = FILTLONG(INPUT_CHECK.out.reads.map { meta, reads -> [meta, [], reads] }).reads
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FILTLONG.out.log)
    } else if (!params.adapter_trimming && !params.quality_filtering) {
        ch_processed_reads = INPUT_CHECK.out.reads
    } else {
        PORECHOP_ABI(INPUT_CHECK.out.reads)
        ch_clipped_reads = PORECHOP_ABI.out.reads.map { meta, reads -> [ meta + [single_end: true], reads ] }
        ch_processed_reads = FILTLONG(ch_clipped_reads.map { meta, reads -> [ meta, [], reads ] }).reads
        ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_ABI.out.log)
        ch_multiqc_files = ch_multiqc_files.mix(FILTLONG.out.log)
    }



    NANOPLOT2(ch_processed_reads)

    // Main analysis
    EMU_ABUNDANCE(ch_processed_reads)
    ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions.first())

    if (params.run_krona) {
        KRONA_KTIMPORTTAXONOMY(EMU_ABUNDANCE.out.report, file(params.krona_taxonomy_tab, checkIfExists: true))
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
    }

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))


/*
    //
    // MODULE: MultiQC Preproccessed
    //
    workflow_summary    = WorkflowGmsemu.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowGmsemu.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // testing other tools
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT1.out.txt.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT2.out.txt.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

*/

/*
    // Prepare inputs for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

    // Generate workflow summary for MultiQC
    def workflow_summary = paramsSummaryMultiqc(workflow, params)
    Channel.value(workflow_summary)
        .collectFile(name: 'workflow_summary_mqc.yaml')
        .set { ch_workflow_summary_yaml }

    // Collect all MultiQC inputs
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary_yaml)
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect())
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT1.out.txt.collect())

    // Run MultiQC
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo
    )
    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions

*/
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
