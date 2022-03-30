find PF* -maxdepth 0 -type d > list_of_fams
for fam in `cat list_of_fams`; do
	for mode in sp; do
		for ref_aln in mtmalign 3dcoffee 3dcoffee_TMalign; do
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode" >> selected_comparisons_ref_"$ref_aln"_"$mode".txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".avg >> selected_comparisons_ref_"$ref_aln"_avg_"$mode".txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".pair >> selected_comparisons_ref_"$ref_aln"_pair_"$mode".txt

			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".same_state >> selected_comparisons_ref_"$ref_aln"_"$mode".same_state.txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".avg.same_state >> selected_comparisons_ref_"$ref_aln"_avg_"$mode".same_state.txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".pair.same_state >> selected_comparisons_ref_"$ref_aln"_pair_"$mode".same_state.txt

			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".without_loops >> selected_comparisons_ref_"$ref_aln"_"$mode".without_loops.txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".avg.without_loops >> selected_comparisons_ref_"$ref_aln"_avg_"$mode".without_loops.txt
			cat "$fam"/results/"$fam"_selected_ref_"$ref_aln"."$mode".pair.without_loops >> selected_comparisons_ref_"$ref_aln"_pair_"$mode".without_loops.txt
		done
	done

	cat "$fam"/results/"$fam"_selected.tcs >> selected_comparisons_tcs.txt
	cat "$fam"/results/"$fam"_selected.tcs.avg >> selected_comparisons_tcs_avg.txt
	cat "$fam"/results/"$fam"_selected.irmsd >> selected_comparisons_irmsd.txt
	cat "$fam"/results/"$fam"_selected.nirmsd >> selected_comparisons_nirmsd.txt
	cat "$fam"/results/"$fam"_selected.nirmsd.avg >> selected_comparisons_nirmsd_avg.txt
	cat "$fam"/results/"$fam"_selected.nirmsd.pair >> selected_comparisons_nirmsd_pair.txt
	cat "$fam"/results/"$fam"_selected.eval >> selected_comparisons_eval.txt

	cat "$fam"/results/"$fam"_selected.pid >> selected_comparisons_pid.txt
	cat "$fam"/results/"$fam"_selected.pid.avg >> selected_comparisons_pid.avg.txt
	cat "$fam"/results/"$fam"_selected.pid.pair >> selected_comparisons_pid.pair.txt
	cat "$fam"/results/"$fam"_selected.len >> selected_comparisons_aln_length.txt

	cat "$fam"/results/"$fam"_alphafold_vs_ref_pdb_comparison_selected.tsv >> alphafold_vs_ref_pdb_comparison_selected.tsv
done

paste list_of_fams selected_comparisons_ref_3dcoffee_sp.txt > test
mv test selected_comparisons_ref_3dcoffee_sp.txt
#paste list_of_fams selected_comparisons_ref_3dcoffee_tc.txt > test
#mv test selected_comparisons_ref_3dcoffee_tc.txt
paste list_of_fams selected_comparisons_ref_3dcoffee_TMalign_sp.txt > test
mv test selected_comparisons_ref_3dcoffee_TMalign_sp.txt
#paste list_of_fams selected_comparisons_ref_3dcoffee_TMalign_tc.txt > test
#mv test selected_comparisons_ref_3dcoffee_TMalign_tc.txt
paste list_of_fams selected_comparisons_ref_mtmalign_sp.txt > test
mv test selected_comparisons_ref_mtmalign_sp.txt
#paste list_of_fams selected_comparisons_ref_mtmalign_tc.txt > test
#mv test selected_comparisons_ref_mtmalign_tc.txt

paste list_of_fams selected_comparisons_tcs.txt > test
mv test selected_comparisons_tcs.txt
paste list_of_fams selected_comparisons_irmsd.txt > test
mv test selected_comparisons_irmsd.txt
paste list_of_fams selected_comparisons_nirmsd.txt > test
mv test selected_comparisons_nirmsd.txt
paste list_of_fams selected_comparisons_eval.txt > test
mv test selected_comparisons_eval.txt