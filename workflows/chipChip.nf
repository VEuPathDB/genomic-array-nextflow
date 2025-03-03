#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow CHIPCHIP {
    take:
    samples


    main:
    rawToGenomeCoordinates(samples)

    peakFinderAndSmoother(rawToGenomeCoordinates.out)

    peakResults(peakFinderAndSmoother.out)

    smoothedResults(peakFinderAndSmoother.out)

    peakResults.out.studyConfig.collectFile(name: "insert_study_results", storeDir: params.outDir, keepHeader: true, skip: 1)
    
}

process rawToGenomeCoordinates {
    container "jbrestel/bioperl"

    input:
    tuple val(meta), path(chipChipFile)

    output:
    tuple val(meta), path("genomeCoords.txt")

    script:
    """
    TransformRawDataToGenomeCoordinates  --inputFile ${params.input}/${chipChipFile} \\
                                         --outputFile genomeCoords.txt \\
                                         --probesBamFile '${params.platformBamFile}'
    """
}


process peakFinderAndSmoother {
    container "jbrestel/genomicarray"


    input:
    tuple val(meta), path(genomeCoordinates)

    output:
    tuple val(meta), path("reformatted_smoothed.txt"), path("reformatted_peaks.txt")

    script:
    """
    java -Xmx2000m -classpath /app/ChIPChipPeakFinder.jar org.apidb.ggtools.array.ChIP_Chip_Peak_Finder $genomeCoordinates  peaks.txt smoothed.txt ${params.peakFinderArgs}

    reformatPeaks.pl  --inputFile peaks.txt --outputFile reformatted_peaks.txt
    reformatSmoothedProfiles.pl --inputFile smoothed.txt --outputFile reformatted_smoothed.txt
    """
}


process peakResults {
    container = 'biocontainers/tabix:v1.9-11-deb_cv1'

    publishDir params.outDir, mode: 'copy', pattern: "*.gz"
    publishDir params.outDir, mode: 'copy', pattern: "*.tbi"
    publishDir params.outDir, mode: 'copy', pattern: "*.txt"

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

process smoothedResults {
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

    publishDir params.outDir, mode: 'copy', pattern: "*.bw"

    input:
    tuple val(meta), path(smoothed), path(peaks)


    output:
    tuple val(meta), path("*.bw")

    script:
    """

    cat $smoothed | tail -n +2 | sort -k1,1 -k2,2n >smoothed.bed
    bedGraphToBigWig smoothed.bed ${params.seqSizeFile} ${meta.id}_smoothed.bw
    """

}
