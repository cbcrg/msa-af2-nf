#find PF* -maxdepth 0 -type d > list_of_fams
for fam in `cat list_of_fams`; do
	cd "$fam"/results/
	for aln in mtmalign_alphafold ginsi psicoffee mtmalign; do
		t_coffee -other_pg seq_reformat -in "$fam"_selected_ref_"$aln".fa_aln -output sim | grep TOT | awk -F '\t' -v var=$fam '{print var"\t"$4}' > "$fam"_selected_ref_"$aln".pid
		t_coffee -other_pg seq_reformat -in "$fam"_selected_ref_"$aln".fa_aln -output sim | grep AVG | head -n -1 | awk -F '\t' -v var=$fam '{print $3"\t"$5"\t"var}' > "$fam"_selected_ref_"$aln".pid.avg
		t_coffee -other_pg seq_reformat -in "$fam"_selected_ref_"$aln".fa_aln -output sim | grep BOT | awk -v var=$fam '{print $5"\t"$6"\t"$7"\t"var}' > "$fam"_selected_ref_"$aln".pid.pair
	done
	paste "$fam"_selected_ref_ginsi.pid "$fam"_selected_ref_psicoffee.pid "$fam"_selected_ref_mtmalign.pid "$fam"_selected_ref_mtmalign_alphafold.pid > "$fam"_selected.pid
	paste "$fam"_selected_ref_ginsi.pid.avg "$fam"_selected_ref_psicoffee.pid.avg "$fam"_selected_ref_mtmalign.pid.avg "$fam"_selected_ref_mtmalign_alphafold.pid.avg > "$fam"_selected.pid.avg
	paste "$fam"_selected_ref_ginsi.pid.pair "$fam"_selected_ref_psicoffee.pid.pair "$fam"_selected_ref_mtmalign.pid.pair "$fam"_selected_ref_mtmalign_alphafold.pid.pair > "$fam"_selected.pid.pair
	cd ../../
done

cat PF*/results/PF*_selected.pid > selected_comparisons_pid.txt
cat PF*/results/PF*_selected.pid.avg > selected_comparisons_pid.avg.txt
cat PF*/results/PF*_selected.pid.pair > selected_comparisons_pid.pair.txt
