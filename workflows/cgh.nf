#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { rawToGenomeCoordinates } from '../modules/local/rawToGenomeCoordinates.nf'
include { wiggleResults } from '../modules/local/wiggleResults.nf'
include { bedtoolsMerge } from '../modules/local/bedtoolsMerge.nf'

workflow CGH {
    take:
    samples


    main:
    rawToGenomeCoordinates(samples)

    bedtoolsMerge(rawToGenomeCoordinates.out)
    wiggleResults(bedtoolsMerge.out, params.seqSizeFile)

}

