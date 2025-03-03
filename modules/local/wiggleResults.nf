#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process wiggleResults {
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

    publishDir params.outDir, mode: 'copy', pattern: "*.bw"

    input:
    tuple val(meta), path(smoothed)


    output:
    tuple val(meta), path("*.bw")

    script:
    """
    cat $smoothed |  sort -k1,1 -k2,2n >smoothed.bed
    bedGraphToBigWig smoothed.bed ${params.seqSizeFile} ${meta.id}.bw
    """

}
