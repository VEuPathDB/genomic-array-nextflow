#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { rawToGenomeCoordinates } from '../modules/local/rawToGenomeCoordinates.nf'
include { wiggleResults } from '../modules/local/wiggleResults.nf'

workflow CHIPCHIP {
    take:
    samples


    main:
    rawToGenomeCoordinates(samples)

    peakFinderAndSmoother(rawToGenomeCoordinates.out)

    peakResults(peakFinderAndSmoother.out)

    adjustedTuple = peakFinderAndSmoother.out.map { meta, smoothed, peaks ->
        tuple (meta, smoothed)
    }

    wiggleResults(adjustedTuple, params.seqSizeFile)

    peakResults.out.studyConfig.collectFile(name: "insert_study_results", storeDir: params.outDir, keepHeader: true, skip: 1)
    
}


process peakFinderAndSmoother {

    // TODO:  Make a veupath version of this from Dockerfile in this repo
    container 'docker.io/jbrestel/genomicarray:latest'

    publishDir params.outDir, mode: 'copy', pattern: "*.peaks"

    input:
    tuple val(meta), path(genomeCoordinates)

    output:
    tuple val(meta), path("reformatted_smoothed_noheader.txt"), path("*.peaks")

    script:
    """
    java -Xmx2000m -classpath /app/ChIPChipPeakFinder.jar org.apidb.ggtools.array.ChIP_Chip_Peak_Finder $genomeCoordinates  peaks.txt smoothed.txt ${params.peakFinderArgs}

    reformatPeaks.pl  --inputFile peaks.txt --outputFile ${meta.id}.peaks

    reformatSmoothedProfiles.pl --inputFile smoothed.txt --outputFile reformatted_smoothed.txt
    tail -n +2 reformatted_smoothed.txt > reformatted_smoothed_noheader.txt
    """
}


process peakResults {
    container = 'biocontainers/tabix:v1.9-11-deb_cv1'

    publishDir params.outDir, mode: 'copy', pattern: "*.gz"
    publishDir params.outDir, mode: 'copy', pattern: "*.tbi"


    input:
    tuple val(meta), path(smoothed), path(peaks)

    output:
    tuple val(meta), path("*.bed.gz*")
    path("study.config"), emit: studyConfig

    script:
    """

    tail -n +2 $peaks | cut -f 1,2,3,5 | sort -k1,1 -k2,2n >${meta.id}_peaks.bed
    bgzip ${meta.id}_peaks.bed
    tabix -p bed ${meta.id}_peaks.bed.gz

    writeStudyConfig --file ${meta.id}_peaks.txt --outputFile study.config --name '${meta.id}_peaks (ChIP-chip)' --protocol 'chipChipPeaks' --sourceIdType segment --profileSetName '${params.profileSetName}'

    """

}

