#!/usr/env nextflow

nextflow.enable.dsl=2

process split_multi_fasta {
	tag "${fam_name}"
	publishDir "${fam_name}/results", mode: 'copy'

	input: 
	tuple val(fam_name), path(fasta_file)

	output:
	tuple val(fam_name), path ("*.fasta")

	"""
	for i in `grep ">" ${fasta_file} | tr -d ">"`; do
		seq_name=`echo \$i | tr '/' '_'`
		t_coffee -other_pg seq_reformat -in ${fasta_file} -action +extract_seq_list \$i > \$seq_name.fasta
	done
	"""
}

process run_alphafold2 {
	tag "${fasta.baseName}-${fam_name}"
	publishDir "${fam_name}/results", mode: 'copy' 

	input:
	tuple val(fam_name), path(fasta)

	output:
	path ("*"), emit: af2_output
	tuple val(fam_name), path ("${fasta.baseName}.alphafold.pdb"), emit: af2_models
 
	"""
	export CUDA_VISIBLE_DEVICES=0
	python3 /app/alphafold/run_alphafold.py --fasta_paths=${fasta} --max_template_date=2020-05-14 --preset=full_dbs --output_dir=\$PWD --model_names=model_1,model_2,model_3,model_4,model_5 --data_dir=/db/ --uniref90_database_path=/db/uniref90/uniref90.fasta --mgnify_database_path=/db/mgnify/mgy_clusters.fa --bfd_database_path=/db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniclust30_database_path=/db/uniclust30/uniclust30_2018_08/uniclust30_2018_08 --pdb70_database_path=/db/pdb70/pdb70 --template_mmcif_dir=/db/pdb_mmcif/mmcif_files --obsolete_pdbs_path=/db/pdb_mmcif/obsolete.dat
	cp "${fasta.baseName}"/ranked_0.pdb ./"${fasta.baseName}".alphafold.pdb

	"""
}
