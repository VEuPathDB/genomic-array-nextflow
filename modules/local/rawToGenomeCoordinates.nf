#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process rawToGenomeCoordinates {
    container 'docker.io/veupathdb/bioperl:latest'

    input:
    tuple val(meta), path(inputFile)

    output:
    tuple val(meta), path("genomeCoords.txt")

    script:
    """
    TransformRawDataToGenomeCoordinates  --inputFile ${params.input}/${inputFile} \\
                                         --outputFile genomeCoords.txt \\
                                         --probesBamFile '${params.platformBamFile}' \\
                                         --assayType ${params.assayType}
    """
}
