#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include {run_seq_aln; run_psicoffee; run_struct_aln; run_comparison; run_evaluation; TCS_filtering; dssp_to_fasta} from "./modules/run_aln_comp_and_eval.nf"
include {run_alphafold2; split_multi_fasta} from "./modules/run_af2.nf"
include {run_seq2maps; run_dmpfold} from "./modules/run_dmpfold.nf"

params.db = "$HOME/db/uniref50.fasta"
params.list = "PF*/*.list"
params.input_fasta = "PF*/*selected_ref.fa"
params.template = "PF*/*selected_ref.template_list"
params.pdbs = "PF*/*.pdb"
params.AF2 = "PF*/AF2/*.pdb"
params.pdb_for_dssp = "PF*/*[0-9].pdb"

log.info """\
	Predictions MSA Analysis - version 0.1
	=====================================
	Input sequences (FASTA)			: ${params.input_fasta}
	Input lists of sequences		: ${params.list}
	Input template lists			: ${params.template}
	Input PDB structures			: ${params.pdbs}
	Input path to Database for PSI-Coffee	: ${params.db}
	"""
	.stripIndent()

Channel.fromPath(params.list).map { it -> [it.toString().split('\\/')[-2],it]}.set{lists}
Channel.fromPath(params.input_fasta).map { it -> [it.toString().split('\\/')[-2],it]}.set{seq_input}
Channel.fromPath(params.template).map { it -> [it.toString().split('\\/')[-2],it]}.set{templates}
Channel.fromPath(params.pdbs).map { it -> [it.toString().split('\\/')[-2],it]}.groupTuple().set{structures}
Channel.fromPath(params.AF2).map { it -> [it.toString().split('\\/')[-3],it]}.groupTuple().set{af2_models}
Channel.fromPath(params.pdb_for_dssp).map { it -> [it.toString().split('\\/')[-2],it]}.groupTuple().set{sec_struct}


workflow analysis {
	dssp_to_fasta(sec_struct)
	run_seq_aln(seq_input)
	run_psicoffee(params.db,seq_input)
	lists.combine(seq_input, by: 0).combine(templates, by: 0).combine(structures, by: 0).combine(af2_models, by: 0).set{struct_input}
	run_struct_aln(struct_input)
	run_seq_aln.out.seq_aln_output.combine(run_psicoffee.out.psicoffee_output, by: 0).combine(dssp_to_fasta.out.dssp_output, by: 0).combine(run_struct_aln.out.struct_output, by: 0).combine(templates, by: 0).combine(structures, by: 0).combine(af2_models, by: 0).set{comp_and_eval_input}
	run_comparison(comp_and_eval_input)
	run_evaluation(comp_and_eval_input)
	seq_input.combine(templates, by: 0).combine(run_evaluation.out.tcs_avg, by: 0).combine(dssp_to_fasta.out.dssp_output, by: 0).combine(structures, by: 0).combine(af2_models, by: 0).set{tcs_filter_input}
	TCS_filtering(params.db,tcs_filter_input)
}

workflow struct_pred_and_analysis {
	dssp_to_fasta(sec_struct)
	split_multi_fasta(seq_input)
	run_alphafold2(split_multi_fasta.out.transpose())
	run_seq2maps(split_multi_fasta.out.transpose())
	run_dmpfold(run_seq2maps.out.dmpfold_input)
	run_seq_aln(seq_input)
	run_psicoffee(params.db,seq_input)
	lists.combine(seq_input, by: 0).combine(templates, by: 0).combine(structures, by: 0).combine(run_alphafold2.out.af2_models.mix(run_dmpfold.out.dmpfold_pdbs).groupTuple(), by: 0).set{struct_input}
	run_struct_aln(struct_input)
	run_seq_aln.out.seq_aln_output.combine(run_psicoffee.out.psicoffee_output, by: 0).combine(dssp_to_fasta.out.dssp_output, by: 0).combine(run_struct_aln.out.struct_output, by: 0).combine(templates, by: 0).combine(structures, by: 0).combine(run_alphafold2.out.af2_models.mix(run_dmpfold.out.dmpfold_pdbs).groupTuple(), by: 0).set{comp_and_eval_input}
	run_comparison(comp_and_eval_input)
	run_evaluation(comp_and_eval_input)
	seq_input.combine(templates, by: 0).combine(run_evaluation.out.tcs_avg, by: 0).combine(dssp_to_fasta.out.dssp_output, by: 0).combine(structures, by: 0).combine(run_alphafold2.out.af2_models.mix(run_dmpfold.out.dmpfold_pdbs).groupTuple(), by: 0).set{tcs_filter_input}
	TCS_filtering(params.db,tcs_filter_input)
}

workflow {
	analysis()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
