for i in `cat !{list}`; do echo -e `grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`'.pdb'; done > "!{fam}"_selected_ref.mtmalign
#	for i in `cat !{list}`; do echo -e ">"$i _P_ "${i//\//_}"'.pdb'; done > "!{fam}"_selected_ref.template_list
	mTM-align -i "!{fam}"_selected_ref.mtmalign -o "!{fam}"_selected_ref_mtmalign
	for i in `cat !{list}`; do
	ref_pdb=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`
	ref_seq=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $1}'`
	sed -i "s|>$ref_pdb\.pdb|$ref_seq|g" result.fasta
	done
	mv result.fasta "!{fam}"_selected_ref_mtmalign.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_mtmalign.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_mtmalign.fa_aln
	for aln in sap TMalign; do
	t_coffee !{fasta} -method "$aln"_pair -template_file !{fam}_selected_ref.template_list -out_lib !{fam}_selected_ref."$aln".lib
	done
	t_coffee !{fasta} -lib !{fam}_selected_ref.sap.lib !{fam}_selected_ref.TMalign.lib -output fasta_aln -outfile !{fam}_selected_ref_3dcoffee.fa_aln -thread 4
	sort_fasta.py !{fam}_selected_ref_3dcoffee.fa_aln > tmp
	mv tmp !{fam}_selected_ref_3dcoffee.fa_aln
	t_coffee !{fasta} -lib !{fam}_selected_ref.TMalign.lib -output fasta_aln -outfile !{fam}_selected_ref_3dcoffee_TMalign.fa_aln -thread 4
	sort_fasta.py !{fam}_selected_ref_3dcoffee_TMalign.fa_aln > tmp
	mv tmp !{fam}_selected_ref_3dcoffee_TMalign.fa_aln

	for x in alphafold; do
		for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
		for i in `cat !{list}`; do echo -e ">"$i "_P_" "${i//\//_}"."$x"'.pdb'; done > !{fam}_selected_ref_"$x".template_list
		for i in `cat !{list}`; do echo -e "${i//\//_}"'.'$x'.pdb'; done > !{fam}_selected_ref_"$x".mtmalign
		mTM-align -i !{fam}_selected_ref_"$x".mtmalign -o !{fam}_selected_ref_mtmalign_"$x"
		for i in `cat !{list}`; do
			pdb_name="${i//\//_}.$x.pdb"
			sed -i "s|$pdb_name|$i|g" result.fasta
		done
		mv result.fasta "!{fam}"_selected_ref_mtmalign_"$x".fa_aln
		sort_fasta.py "!{fam}"_selected_ref_mtmalign_"$x".fa_aln > tmp
		mv tmp "!{fam}"_selected_ref_mtmalign_"$x".fa_aln
		for aln in sap TMalign; do
			t_coffee !{fasta} -method "$aln"_pair -template_file !{fam}_selected_ref_"$x".template_list -out_lib !{fam}_selected_ref_"$x"."$aln".lib
		done
		t_coffee !{fasta} -lib "!{fam}"_selected_ref_"$x".sap.lib "!{fam}"_selected_ref_"$x".TMalign.lib -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln -thread 4
		sort_fasta.py "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln > tmp
		mv tmp "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln
		t_coffee !{fasta} -lib "!{fam}"_selected_ref_"$x".TMalign.lib -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln -thread 4
		sort_fasta.py "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln > tmp
		mv tmp "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln

		for i in `cat !{list}`; do 
			ref_pdb=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`
			echo -e $i"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "RMSD of" | awk '{print $6}'`"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "^TM-score" | awk '{print $3}'`"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "GDT-TS" | awk '{print $2}'`"\t""!{fam}"; done >> "!{fam}"_"$x"_vs_ref_pdb_comparison_selected.tsv
	done
