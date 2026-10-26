nextflow.enable.dsl=2
include { POLICY_CHECK as QC; POLICY_CHECK as NETWORK_FIT; POLICY_CHECK as DEXSEQ_COUNT } from './resource_policy_task'
workflow RESOURCE_POLICY {
    QC(Channel.of(
        [id:'gene-body', workflow_id:'flow.qc.genebody', execution_group:'qc'],
        [id:'sequencing-summary', workflow_id:'flow.qc.sequencing', execution_group:'qc'],
        [id:'annotation', workflow_id:'flow.reference.bed12', execution_group:'reference']
    ))
    NETWORK_FIT(Channel.of([id:'parallel-fixture', execution_group:'network_fit',
        resources:[cpus:2, memory:'16 MB', time:'1m']]))
    DEXSEQ_COUNT(Channel.of([id:'serial-upstream-fixture', execution_group:'test',
        resources:[cpus:2, memory:'16 MB', time:'1m']]))
}
workflow { RESOURCE_POLICY() }
