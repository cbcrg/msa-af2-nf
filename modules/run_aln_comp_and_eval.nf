#!/usr/env nextflow

nextflow.enable.dsl=2

process run_seq_aln {
	tag "${fam}_seq_aln"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_famsa.fa_aln"), path("${fam}_selected_ref_ginsi.fa_aln"), path("${fam}_selected_ref_msaprobs.fa_aln"), path("${fam}_selected_ref_tcoffee.fa_aln"), path("${fam}_selected_ref_tcoffee.lib"), path("${fam}_selected_ref_ginsi.lib"), emit: seq_aln_output

	shell:
	'''
	famsa !{fasta} "!{fam}"_selected_ref_famsa.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_famsa.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_famsa.fa_aln	
	ginsi !{fasta} > "!{fam}"_selected_ref_ginsi.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_ginsi.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_ginsi.fa_aln
	msaprobs -o "!{fam}"_selected_ref_msaprobs.fa_aln !{fasta}
	sort_fasta.py "!{fam}"_selected_ref_msaprobs.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_msaprobs.fa_aln
	t_coffee !{fasta} -output fasta_aln -out_lib "!{fam}"_selected_ref_tcoffee.lib -outfile "!{fam}"_selected_ref_tcoffee.fa_aln -thread 4
	sort_fasta.py "!{fam}"_selected_ref_tcoffee.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_tcoffee.fa_aln
	t_coffee -in !{fasta} mafftginsi_pair -output fasta_aln -out_lib "!{fam}"_selected_ref_ginsi.lib -outfile "!{fam}"_selected_ref_ginsi_pair.fa_aln -thread 4
	'''
}

process run_psicoffee {
	tag "${fam}_psicoffee_aln"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4
	
	input:
    val db
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_psicoffee.fa_aln"), path("${fam}_selected_ref_psicoffee.lib"), emit: psicoffee_output

	shell:
    '''
	t_coffee -seq !{fasta} -special_mode psicoffee -blast_server LOCAL -protein_db !{db} -psitrim 100 -psiJ 3 -prot_min_cov 90 -prot_max_sim 100 -prot_min_sim 0 -output fasta_aln -out_lib "!{fam}"_selected_ref_psicoffee.lib -outfile "!{fam}"_selected_ref_psicoffee.fa_aln -thread 4
	sort_fasta.py "!{fam}"_selected_ref_psicoffee.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_psicoffee.fa_aln

	'''
}

process run_struct_aln {
	tag "${fam}_struct_aln"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(fasta), path(temp), path(pdbs), path(af2_models)

	output:
	tuple val(fam), path("*.fa_aln"), path("*.lib"), path("*.template_list"), emit: struct_output
	path("*_vs_ref_pdb_comparison_selected.tsv"), emit: pdb_comp

	shell:
	template 'struct_aln.sh'
}

process run_comparison {
	tag "${fam}_comparison"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

    input:
	tuple val(fam), path(famsa_aln), path(ginsi_aln), path(msaprobs_aln), path(tcoffee_aln), path(tcoffee_lib), path(ginsi_lib), path(psicoffee_aln), path(psicoffee_lib), path(dssp_fasta), path(struct_fasta_files), path(struct_libs), path(struct_templates), path(temps), path(pdbs), path(af2_models) //file(deepblast_aln), file(deepblast_lib),

	output:
	tuple path("*.sp*"), path("*.tc*"), path("*.avg*"), path("*.pair*"), emit: comp_output

	shell:
	template 'comparisons.sh'
}

process run_evaluation {
	tag "${fam}"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(famsa_aln), path(ginsi_aln), path(msaprobs_aln), path(tcoffee_aln), path(tcoffee_lib), path(ginsi_lib), path(psicoffee_aln), path(psicoffee_lib), path(dssp_fasta), path(struct_fasta_files), path(struct_libs), path(struct_templates), path(temps), path(pdbs), path(af2_models) //file(deepblast_aln), file(deepblast_lib),

	output:
	path("*"), emit: eval_output
	tuple val(fam), path("${fam}_selected.tcs.avg"), emit: tcs_avg

	shell:
	template 'evaluations.sh'
}

process TCS_filtering {
	tag "${fam}"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	val db
	tuple val(fam), path(fasta), path(template), path(avg_tcs), path(dssp_fasta), path(nat_pdbs), path(af2_models)

	output:
	path("*"), emit: tcs_filter_out

	shell:
	template 'filter_TCS.sh'
}

process dssp_to_fasta {
	tag "${fam}"
	publishDir "${fam}/results", mode: 'copy'

	input:
	tuple val(fam), path(sec_struct)

	output:
	path ("*"), emit: dssp_to_fasta_out
	tuple val(fam), path("${fam}.dssp"), emit: dssp_output

	shell:
	'''
	dssp2fasta.py
	cat *.dssp > "!{fam}".dssp
	'''
}