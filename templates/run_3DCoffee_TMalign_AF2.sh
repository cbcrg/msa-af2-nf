for x in alphafold; do
		for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
		for i in `cat !{list}`; do echo -e ">"$i "_P_" "${i//\//_}"."$x"'.pdb'; done > !{fam}_selected_ref_"$x".template_list
		for aln in TMalign; do
			t_coffee !{fasta} -method "$aln"_pair -template_file !{fam}_selected_ref_"$x".template_list -out_lib !{fam}_selected_ref_"$x"."$aln".lib
		done
		t_coffee !{fasta} -lib "!{fam}"_selected_ref_"$x".TMalign.lib -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln -thread 4
		sort_fasta.py "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln > tmp
		mv tmp "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln
	done