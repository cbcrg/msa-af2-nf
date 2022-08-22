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

process run_famsa {
	tag "${fam}_famsa"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_famsa.fa_aln"), emit: famsa_aln_output

	shell:
	'''
	famsa !{fasta} "!{fam}"_selected_ref_famsa.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_famsa.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_famsa.fa_aln
	'''
}

process run_ginsi {
	tag "${fam}_ginsi"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_ginsi.fa_aln"), emit: ginsi_aln_output

	shell:
	'''
	ginsi !{fasta} > "!{fam}"_selected_ref_ginsi.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_ginsi.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_ginsi.fa_aln
	'''
}

process run_msaprobs {
	tag "${fam}_msaprobs"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_msaprobs.fa_aln"), emit: msaprobs_aln_output

	shell:
	'''
	msaprobs -o "!{fam}"_selected_ref_msaprobs.fa_aln !{fasta}
	sort_fasta.py "!{fam}"_selected_ref_msaprobs.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_msaprobs.fa_aln
	'''
}

process run_tcoffee {
	tag "${fam}_tcoffee"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_tcoffee.fa_aln"), path("${fam}_selected_ref_tcoffee.lib"), emit: tcoffee_aln_output

	shell:
	'''
	t_coffee !{fasta} -output fasta_aln -out_lib "!{fam}"_selected_ref_tcoffee.lib -outfile "!{fam}"_selected_ref_tcoffee.fa_aln -thread 4
	sort_fasta.py "!{fam}"_selected_ref_tcoffee.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_tcoffee.fa_aln
	'''
}

process run_ginsi_pair {
	tag "${fam}_ginsi_pair"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta)

	output:
	tuple val(fam), path("${fam}_selected_ref_ginsi.lib"), emit: ginsi_pair_output

	shell:
	'''
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
	path "*_vs_ref_pdb_comparison_selected.tsv", emit: pdb_comp

	shell:
	template 'struct_aln.sh'
}

process run_mTMalign_PDB {
	tag "${fam}_mTMalign_PDB"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(fasta), path(temp), path(pdbs)

	output:
	tuple val(fam), path("${fam}_selected_ref_mtmalign.fa_aln"), emit: mtmalign_pdb_output

	shell:
    template 'run_mTMalign_PDB.sh'
}

process run_3DCoffee_PDB {
	tag "${fam}_3DCoffee_PDB"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta), path(temp), path(pdbs)

	output:
	tuple val(fam), path("${fam}_selected_ref_3dcoffee.fa_aln"), path("*.lib"), emit: output_3dcoffee_pdb

	shell:
	template 'run_3DCoffee_PDB.sh'
}

process run_3DCoffee_TMalign_PDB {
	tag "${fam}_3DCoffee_TMalign_PDB"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(fasta), path(temp), path(pdbs)

	output:
	tuple val(fam), path("${fam}_selected_ref_3dcoffee_TMalign.fa_aln"), emit: output_3dcoffee_tmalign_pdb

	shell:
	template 'run_3DCoffee_TMalign_PDB.sh'
}

process run_mTMalign_AF2 {
	tag "${fam}_mTMalign_AF2"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(fasta), path(af2_models)

	output:
	tuple val(fam), path("*.fa_aln"), emit: mtmalign_af2_output

	shell:
	template 'run_mTMalign_AF2.sh'
}

process run_3DCoffee_AF2 {
	tag "${fam}_3DCoffee_AF2"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(fasta), path(af2_models)

	output:
	tuple val(fam), path("*.fa_aln"), path("*.lib"), path("*.template_list"), emit: output_3dcoffee_af2

	shell:
	template 'run_3DCoffee_AF2.sh'
}

process run_3DCoffee_TMalign_AF2 {
	tag "${fam}_3DCoffee_TMalign_AF2"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(fasta), path(af2_models)

	output:
	tuple val(fam), path("*.fa_aln"), emit: output_3dcoffee_tmalign_af2

	shell:
	template 'run_3DCoffee_TMalign_AF2.sh'
}

process run_structure_comparison {
	tag "${fam}_structure_comparison"
    publishDir "${fam}/results", mode: 'copy'
	cpus 4

	input:
	tuple val(fam), path(list), path(temp), path(pdbs), path(af2_models)

	output:
	path "*_vs_ref_pdb_comparison_selected.tsv", emit: pdb_comp

	shell:
	template 'run_structure_comparison.sh'
}

process run_comparison {
	tag "${fam}_comparison"
	publishDir "${fam}/results", mode: 'copy'
	cpus 4

    input:
	tuple val(fam), path(famsa_aln),
        path(ginsi_aln),
        path(msaprobs_aln),
        path(tcoffee_aln),
        path(tcoffee_lib),
        path(ginsi_lib),
        path(psicoffee_aln),
        path(psicoffee_lib),
        path(dssp_fasta),
		path(mtmalign_pdb),
        path(aln_3dcoffee_pdb),
		path(aln_3dcoffee_pdb_libs),
        path(aln_3dcoffee_tmalign_pdb),
        path(mtmalign_af2),
        path(aln_3dcoffee_af2),
        path(aln_3dcoffee_af2_libs),
        path(aln_3dcoffee_af2_temps),
        path(aln_3dcoffee_tmalign_af2),
        path(temps),
        path(pdbs),
        path(af2_models)

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
	tuple val(fam), path(famsa_aln),
        path(ginsi_aln),
        path(msaprobs_aln),
        path(tcoffee_aln),
        path(tcoffee_lib),
        path(ginsi_lib),
        path(psicoffee_aln),
        path(psicoffee_lib),
        path(dssp_fasta),
		path(mtmalign_pdb),
        path(aln_3dcoffee_pdb),
		path(aln_3dcoffee_pdb_libs),
        path(aln_3dcoffee_tmalign_pdb),
        path(mtmalign_af2),
        path(aln_3dcoffee_af2),
        path(aln_3dcoffee_af2_libs),
        path(aln_3dcoffee_af2_temps),
        path(aln_3dcoffee_tmalign_af2),
        path(temps),
        path(pdbs),
        path(af2_models)

	output:
	path("*"), emit: eval_output
	tuple val(fam), path("${fam}_selected.tcs.avg"), emit: tcs_avg

	shell:
	template 'evaluations.sh'
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