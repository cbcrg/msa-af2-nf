for x in alphafold; do
		for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
		for i in `cat !{list}`; do 
			ref_pdb=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`
			echo -e $i"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "RMSD of" | awk '{print $6}'`"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "^TM-score" | awk '{print $3}'`"\t"`TMscore "${i//\//_}"."$x".pdb "${ref_pdb}".pdb | grep "GDT-TS" | awk '{print $2}'`"\t""!{fam}"; done >> "!{fam}"_"$x"_vs_ref_pdb_comparison_selected.tsv
	done