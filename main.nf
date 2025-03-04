#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHIPCHIP  } from './workflows/chipChip.nf'
include { CGH  } from './workflows/cgh.nf'

workflow {

    sampleRows = Channel.fromPath(params.input + "/" + params.samplesheetFileName)
        .splitCsv( skip:1)

    samples = sampleRows.map { row ->
        fileName = file(params.input + "/" + row[1]);
        return [ [id: row[0] ], fileName ]
    }

    if(params.assayType == "chipChip") {
        CHIPCHIP(samples);
    }

    if(params.assayType == "cghArray") {
        CGH(samples);
    }



}
