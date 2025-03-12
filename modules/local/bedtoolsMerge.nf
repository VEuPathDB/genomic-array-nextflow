#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process bedtoolsMerge {
    container "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"

    input:
    tuple val(meta), path(bed)


    output:
    tuple val(meta), path("merged.bed")

    script:
    """
    cat $bed |  sort -k1,1 -k2,2n -k3,3n >sorted.bed

    bedtools merge -i sorted.bed -c 4 -d 0 -o max > merged.bed
    """


}
