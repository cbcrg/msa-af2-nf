#!/usr/env nextflow

//params.input_fasta = "*/*.fa"
params.input_fasta = "PF*/predictions/PF*/*.fa"
Channel.fromPath(params.input_fasta).map { it -> [it.toString().split('\\/')[-2],it.baseName,it]}.set{fasta_file}


process run_seq2maps {
	tag "${fasta.baseName}-${fam_name}"
	publishDir "${fam_name}/results", mode: 'copy' 

	input:
        tuple val(fam_name), path(fasta)

	output:
	path("*"), emit: run_seq2maps_output
	tuple val(fam_name), val(fasta.baseName), path("${fasta}"), path("${fasta.baseName}.21c"), path("${fasta.baseName}.map"), emit: dmpfold_input

	shell: 
	'''
	awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' !{fasta} > test
        mv test !{fasta}
	csh /DMPfold/seq2maps.csh !{fasta}
	'''
}
process run_dmpfold {
	tag "${fasta.baseName}-${fam_name}"
	publishDir "${fam_name}/results", mode: 'copy' 
        cpus 4

        input:
        tuple val(fam_name), val(seq_name), path(fasta), path(file_21c), path(file_map)

        output:
        path ("*"), emit: dmpfold_output
        tuple val(fam_name), path("${seq_name}.dmpfold.pdb"), emit: dmpfold_pdbs

        shell:

        '''
	bash /DMPfold/run_dmpfold.sh !{fasta} !{file_21c} !{file_map} ./!{seq_name}
	cp ./!{seq_name}/final_1.pdb ./!{seq_name}.dmpfold.pdb
        '''
}
