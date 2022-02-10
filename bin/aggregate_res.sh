find PF* -maxdepth 0 -type d > list_of_fams
for fam in `cat list_of_fams`; do
	for mode in sp tc; do
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

	cat "$fam"/results/"$fam"_dmpfold_vs_ref_pdb_comparison_selected.tsv >> dmpfold_vs_ref_pdb_comparison_selected.tsv
	cat "$fam"/results/"$fam"_alphafold_vs_ref_pdb_comparison_selected.tsv >> alphafold_vs_ref_pdb_comparison_selected.tsv
done

paste list_of_fams selected_comparisons_ref_3dcoffee_sp.txt > test
mv test selected_comparisons_ref_3dcoffee_sp.txt
paste list_of_fams selected_comparisons_ref_3dcoffee_tc.txt > test
mv test selected_comparisons_ref_3dcoffee_tc.txt
paste list_of_fams selected_comparisons_ref_3dcoffee_TMalign_sp.txt > test
mv test selected_comparisons_ref_3dcoffee_TMalign_sp.txt
paste list_of_fams selected_comparisons_ref_3dcoffee_TMalign_tc.txt > test
mv test selected_comparisons_ref_3dcoffee_TMalign_tc.txt
paste list_of_fams selected_comparisons_ref_mtmalign_sp.txt > test
mv test selected_comparisons_ref_mtmalign_sp.txt
paste list_of_fams selected_comparisons_ref_mtmalign_tc.txt > test
mv test selected_comparisons_ref_mtmalign_tc.txt

paste list_of_fams selected_comparisons_tcs.txt > test
mv test selected_comparisons_tcs.txt
paste list_of_fams selected_comparisons_irmsd.txt > test
mv test selected_comparisons_irmsd.txt
paste list_of_fams selected_comparisons_nirmsd.txt > test
mv test selected_comparisons_nirmsd.txt
paste list_of_fams selected_comparisons_eval.txt > test
mv test selected_comparisons_eval.txt

for fam in `cat list_of_fams`; do
        echo -e  $fam"\t"`wc -l "$fam"/"$fam".list | awk '{print $1}'` >> all_tcs_cutoffs_all_fams_selected.list
done
# 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do 
	for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
		for test_aln in famsa tcoffee msaprobs ginsi psicoffee 3dcoffee_dmpfold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_alphafold mtmalign_alphafold; do
			touch all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref"
			touch all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state
			touch all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops
		done
	done
	cp all_tcs_cutoffs_all_fams_selected.list all_tcs_cutoffs_all_fams_selected.list_ref_"$tcs_ref"
	cp all_tcs_cutoffs_all_fams_selected.list all_tcs_cutoffs_all_fams_selected.list_ref_"$tcs_ref".same_state
	cp all_tcs_cutoffs_all_fams_selected.list all_tcs_cutoffs_all_fams_selected.list_ref_"$tcs_ref".without_loops
done

for (( tcs=50; tcs<=95; tcs+=5 )); do
	for fam in `cat list_of_fams`; do
		cd "$fam"/results
		for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do # 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
			for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
				for test_aln in famsa tcoffee msaprobs ginsi psicoffee 3dcoffee_dmpfold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_alphafold mtmalign_alphafold; do
					awk '{ if ($1 == "NA") print $1; else print $2}' "$fam"_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref" >> ../../all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref"
					awk '{ if ($1 == "NA") print $1; else print $2}' "$fam"_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".same_state >> ../../all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".same_state
					awk '{ if ($1 == "NA") print $1; else print $2}' "$fam"_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops >> ../../all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
				done
			done
			cat "$fam"_selected.list_tcs_above_"$tcs"_ref_"$tcs_ref" >> ../../all_fams_selected.list_tcs_above_"$tcs"_ref_"$tcs_ref"
			cat "$fam".list_tcs_above_"$tcs"_ref_"$tcs_ref" >> ../../all_fams.list_tcs_above_"$tcs"_ref_"$tcs_ref"
		done
		cd ../..
	done
	for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do # 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
		paste all_tcs_cutoffs_all_fams_selected.list_ref_"$tcs_ref" all_fams_selected.list_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp2 && mv tmp2 all_tcs_cutoffs_all_fams_selected.list_ref_"$tcs_ref"
	done
	for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do # 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
		for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
			for test_aln in famsa tcoffee msaprobs ginsi psicoffee 3dcoffee_dmpfold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_alphafold mtmalign_alphafold; do
				paste all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp3 && mv tmp3 all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref"
				paste all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".same_state > tmp3 && mv tmp3 all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state
				paste all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops > tmp3 && mv tmp3 all_tcs_cutoffs_all_fams_selected_ref_"$test_aln"_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops
			done
		done
	done	
done


paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_dmpfold_vs_ref_3dcoffee.sp
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_dmpfold_vs_ref_mtmalign.sp

paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.same_state all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.same_state > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_dmpfold_vs_ref_3dcoffee.sp.same_state
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.same_state all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.same_state all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.same_state > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_dmpfold_vs_ref_mtmalign.sp.same_state

paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.without_loops all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_3dcoffee.sp_ref_3dcoffee_dmpfold.without_loops > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_dmpfold_vs_ref_3dcoffee.sp.without_loops
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.without_loops all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.without_loops all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_mtmalign.sp_ref_mtmalign_dmpfold.without_loops > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_dmpfold_vs_ref_mtmalign.sp.without_loops

paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_alphafold_vs_ref_3dcoffee.sp
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_alphafold all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_alphafold all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_mtmalign.sp_ref_mtmalign_alphafold > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_alphafold_vs_ref_mtmalign.sp

paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.same_state all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.same_state > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_alphafold_vs_ref_3dcoffee.sp.same_state
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.same_state all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.same_state all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.same_state > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_alphafold_vs_ref_mtmalign.sp.same_state

paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.without_loops all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_3dcoffee.sp_ref_3dcoffee_alphafold.without_loops > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_alphafold_vs_ref_3dcoffee.sp.without_loops
paste all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.without_loops all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.without_loops all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_mtmalign.sp_ref_mtmalign_alphafold.without_loops > all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_alphafold_vs_ref_mtmalign.sp.without_loops

for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do # 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
	for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
		paste all_tcs_cutoffs_all_fams_selected_ref_famsa_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_msaprobs_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_tcoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref" > all_tcs_cutoffs_all_fams_selected_all_vs_ref_"$ref_aln".sp_ref_"$tcs_ref"
		paste all_tcs_cutoffs_all_fams_selected_ref_famsa_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_msaprobs_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_tcoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state > all_tcs_cutoffs_all_fams_selected_all_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".same_state
		paste all_tcs_cutoffs_all_fams_selected_ref_famsa_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_ginsi_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_msaprobs_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_tcoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops all_tcs_cutoffs_all_fams_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops > all_tcs_cutoffs_all_fams_selected_all_vs_ref_"$ref_aln".sp_ref_"$tcs_ref".without_loops
		for threshold in 75 80 85 90 95; do
			for fam in `cat list_of_fams`; do 
				paste "$fam"/results/"$fam"_selected_ref_famsa_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_ginsi_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_msaprobs_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_tcoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_psicoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" "$fam"/results/"$fam"_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref" > tmp
				cat tmp >> selected_comparisons_ref_"$ref_aln"_pair_sp.txt_tcs_above_"$threshold"_ref_"$tcs_ref"
				rm tmp
			
				paste "$fam"/results/"$fam"_selected_ref_famsa_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_ginsi_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_msaprobs_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_tcoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_psicoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state "$fam"/results/"$fam"_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".same_state > tmp
				cat tmp >> selected_comparisons_ref_"$ref_aln"_pair_sp.txt_tcs_above_"$threshold"_ref_"$tcs_ref".same_state
				rm tmp

				paste "$fam"/results/"$fam"_selected_ref_famsa_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_ginsi_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_msaprobs_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_tcoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_psicoffee_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_3dcoffee_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_mtmalign_dmpfold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_3dcoffee_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_3dcoffee_TMalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops "$fam"/results/"$fam"_selected_ref_mtmalign_alphafold_vs_ref_"$ref_aln".sp.pair_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops > tmp
				cat tmp >> selected_comparisons_ref_"$ref_aln"_pair_sp.txt_tcs_above_"$threshold"_ref_"$tcs_ref".without_loops
				rm tmp
			done
		done
	done
done

#for fam in `cat list_of_fams`; do paste "$fam"/"$fam"_selected_ref_ginsi_vs_ref_3dcoffee.sp.pair_tcs_above_75_ref_3dcoffee_dmpfold "$fam"/"$fam"_selected_ref_psicoffee_vs_ref_3dcoffee.sp.pair_tcs_above_75_ref_3dcoffee_dmpfold "$fam"/"$fam"_selected_ref_3dcoffee_dmpfold_vs_ref_3dcoffee.sp.pair_tcs_above_75_ref_3dcoffee_dmpfold > tmp; cat tmp >> selected_comparisons_ref_3dcoffee_pair_sp.txt_tcs_above_75_ref_3dcoffee_dmpfold; rm tmp; done
#for fam in `cat list_of_fams`; do paste "$fam"/"$fam"_selected_ref_ginsi_vs_ref_mtmalign.sp.pair_tcs_above_75_ref_mtmalign_dmpfold "$fam"/"$fam"_selected_ref_psicoffee_vs_ref_mtmalign.sp.pair_tcs_above_75_ref_mtmalign_dmpfold "$fam"/"$fam"_selected_ref_mtmalign_dmpfold_vs_ref_mtmalign.sp.pair_tcs_above_75_ref_mtmalign_dmpfold > tmp; cat tmp >> selected_comparisons_ref_mtmalign_pair_sp.txt_tcs_above_75_ref_mtmalign_dmpfold; rm tmp; done
