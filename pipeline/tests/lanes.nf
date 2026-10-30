nextflow.enable.dsl=2
include { VALIDATE_INPUTS } from '../modules/local/validate_inputs'
include { STAGE_FASTQ } from '../modules/local/stage_fastq'
include { fastqLaneInputs } from '../workflows/katana_core'

workflow {
    samples=Channel.value(file(params.input))
    raw=samples.splitCsv(header:true).flatMap { row -> [file(row.fastq_1), file(row.fastq_2)] }.collect()
    VALIDATE_INPUTS(samples, Channel.value(file(params.contrasts)), raw,
        Channel.value(file(params.fasta)), Channel.value(file(params.gtf)),
        Channel.value(file("${projectDir}/../bin/katana_validate_inputs.py")),
        file(params.input).parent.toString(), 4, '~ subject + condition', 'analysis')
    validated=VALIDATE_INPUTS.out.samples.splitCsv(header:true).map { row ->
        def inputs=fastqLaneInputs(row)
        tuple([id:row.sample, single_end:false, lane_count:inputs.lane_count], inputs.reads)
    }
    STAGE_FASTQ(validated, Channel.value(file("${projectDir}/../bin/kanata_stage_fastq.py")))
}
