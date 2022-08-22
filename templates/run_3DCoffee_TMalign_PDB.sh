for aln in TMalign; do
	t_coffee !{fasta} -method "$aln"_pair -template_file !{fam}_selected_ref.template_list -out_lib !{fam}_selected_ref."$aln".lib
	done
	t_coffee !{fasta} -lib !{fam}_selected_ref.TMalign.lib -output fasta_aln -outfile !{fam}_selected_ref_3dcoffee_TMalign.fa_aln -thread 4
	sort_fasta.py !{fam}_selected_ref_3dcoffee_TMalign.fa_aln > tmp
	mv tmp !{fam}_selected_ref_3dcoffee_TMalign.fa_aln