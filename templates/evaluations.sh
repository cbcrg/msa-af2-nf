t_coffee -infile "!{fam}"_selected_ref_tcoffee.fa_aln -evaluate -lib "!{fam}"_selected_ref_tcoffee.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_tcoffee.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_tcoffee.tcs
grep -P "/.*:" "!{fam}"_selected_ref_tcoffee.score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_tcoffee.tcs.avg
t_coffee -infile "!{fam}"_selected_ref_psicoffee.fa_aln -evaluate -lib "!{fam}"_selected_ref_psicoffee.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_psicoffee.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_psicoffee.tcs
grep -P "/.*:" "!{fam}"_selected_ref_psicoffee.score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_psicoffee.tcs.avg

t_coffee -infile "!{fam}"_selected_ref_ginsi.fa_aln -evaluate -lib "!{fam}"_selected_ref_ginsi.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_ginsi.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_ginsi.tcs
grep -P "/.*:" "!{fam}"_selected_ref_ginsi.score_ascii | awk -v var=!{fam} '{print $1"\t"var"\t"$3}' > "!{fam}"_selected_ref_ginsi.tcs.avg

t_coffee -infile "!{fam}"_selected_ref_3dcoffee.fa_aln -evaluate -lib "!{fam}"_selected_ref.sap.lib "!{fam}"_selected_ref.TMalign.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_3dcoffee.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_3dcoffee.tcs
grep -P "/.*:" "!{fam}"_selected_ref_3dcoffee.score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_3dcoffee.tcs.avg
t_coffee -infile "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln -evaluate -lib "!{fam}"_selected_ref.TMalign.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_3dcoffee_TMalign.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_3dcoffee_TMalign.tcs
grep -P "/.*:" "!{fam}"_selected_ref_3dcoffee_TMalign.score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_3dcoffee_TMalign.tcs.avg
t_coffee -infile "!{fam}"_selected_ref_mtmalign.fa_aln -evaluate -lib "!{fam}"_selected_ref.TMalign.lib -output score_ascii
grep SCORE "!{fam}"_selected_ref_mtmalign.score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_mtmalign.tcs
grep -P "/.*:" "!{fam}"_selected_ref_mtmalign.score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_mtmalign.tcs.avg
       	
paste "!{fam}"_selected_ref_ginsi.tcs "!{fam}"_selected_ref_tcoffee.tcs "!{fam}"_selected_ref_psicoffee.tcs "!{fam}"_selected_ref_3dcoffee.tcs "!{fam}"_selected_ref_3dcoffee_TMalign.tcs "!{fam}"_selected_ref_mtmalign.tcs > "!{fam}"_selected.tcs
paste "!{fam}"_selected_ref_ginsi.tcs.avg "!{fam}"_selected_ref_tcoffee.tcs.avg "!{fam}"_selected_ref_psicoffee.tcs.avg "!{fam}"_selected_ref_3dcoffee.tcs.avg "!{fam}"_selected_ref_3dcoffee_TMalign.tcs.avg "!{fam}"_selected_ref_mtmalign.tcs.avg > "!{fam}"_selected.tcs.avg

touch "!{fam}"_selected.nirmsd.avg
touch "!{fam}"_selected.nirmsd.pair
touch "!{fam}"_selected.pid.avg
touch "!{fam}"_selected.pid.pair
for pref in "!{fam}"_selected_ref_famsa "!{fam}"_selected_ref_ginsi "!{fam}"_selected_ref_msaprobs "!{fam}"_selected_ref_tcoffee "!{fam}"_selected_ref_psicoffee "!{fam}"_selected_ref_3dcoffee "!{fam}"_selected_ref_3dcoffee_TMalign "!{fam}"_selected_ref_mtmalign; do
	t_coffee -other_pg irmsd "$pref".fa_aln -template_file "!{fam}"_selected_ref.template_list > "$pref".all_irmsd
	grep -P "TOTAL\s+iRMSD" "$pref".all_irmsd | awk '{print $3}' > "$pref".irmsd
	grep -P "TOTAL\s+NiRMSD" "$pref".all_irmsd | awk '{print $3}' > "$pref".nirmsd
	grep -P "\sAVERAGE\s+NiRMS" "$pref".all_irmsd | awk -v var=!{fam} '{print $5"\t"var"\t"$3}' | tr -d "[]" > "$pref".nirmsd.avg
	grep -P "\sPAIRWISE\s+NiRMSD" "$pref".all_irmsd | awk -v var=!{fam} '{print $5"\t"$7"\t"var"\t"$3}' | tr -d "[]" > "$pref".nirmsd.pair
	grep -P "TOTAL\s+EVAL" "$pref".all_irmsd | awk '{print $3}' > "$pref".eval

	cat "$pref".irmsd >> test.irmsd
	cat "$pref".nirmsd >> test.nirmsd
	cat "$pref".eval >> test.eval

	paste "!{fam}"_selected.nirmsd.avg "$pref".nirmsd.avg > new_selected.nirmsd.avg
    mv new_selected.nirmsd.avg "!{fam}"_selected.nirmsd.avg
	paste "!{fam}"_selected.nirmsd.pair "$pref".nirmsd.pair > new_selected.nirmsd.pair
    mv new_selected.nirmsd.pair "!{fam}"_selected.nirmsd.pair

	t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep TOT | awk -F '\t' -v var=!{fam} '{print var"\t"$4}' > "$pref".pid
	t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep AVG | head -n -1 | awk -F '\t' -v var=!{fam} '{print $3"\t"$5"\t"var}' > "$pref".pid.avg
	t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep BOT | awk -v var=!{fam} '{print $5"\t"$6"\t"$7"\t"var}' > "$pref".pid.pair

	cat "$pref".pid >> test.pid

	paste "!{fam}"_selected.pid.avg "$pref".pid.avg > new_selected.pid.avg
    mv new_selected.pid.avg "!{fam}"_selected.pid.avg
	paste "!{fam}"_selected.pid.pair "$pref".pid.pair > new_selected.pid.pair
    mv new_selected.pid.pair "!{fam}"_selected.pid.pair

	echo -e !{fam}"\t"`esl-alistat "$pref".fa_aln | grep "Alignment length" | awk '{print $3}'` >> test.len

done
	
paste -s test.irmsd > "!{fam}"_selected.irmsd
paste -s test.nirmsd > "!{fam}"_selected.nirmsd
paste -s test.eval > "!{fam}"_selected.eval
paste -s test.pid > "!{fam}"_selected.pid
paste -s test.len > "!{fam}"_selected.len
rm test.irmsd test.nirmsd test.eval test.pid test.len

for x in alphafold; do
	for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
	t_coffee -infile "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln -evaluate -lib "!{fam}"_selected_ref_"$x".sap.lib "!{fam}"_selected_ref_"$x".TMalign.lib -output score_ascii
	grep SCORE "!{fam}"_selected_ref_3dcoffee_"$x".score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_3dcoffee_"$x".tcs
	grep -P "/.*:" "!{fam}"_selected_ref_3dcoffee_"$x".score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_3dcoffee_"$x".tcs.avg
	t_coffee -infile "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln -evaluate -lib "!{fam}"_selected_ref_"$x".TMalign.lib -output score_ascii
	grep SCORE "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".tcs
	grep -P "/.*:" "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".tcs.avg
	t_coffee -infile "!{fam}"_selected_ref_mtmalign_"$x".fa_aln -evaluate -lib "!{fam}"_selected_ref_"$x".TMalign.lib -output score_ascii
	grep SCORE "!{fam}"_selected_ref_mtmalign_"$x".score_ascii | tr "=" "\t" | cut -f2 > "!{fam}"_selected_ref_mtmalign_"$x".tcs
	grep -P "/.*:" "!{fam}"_selected_ref_mtmalign_"$x".score_ascii | awk -v var=!{fam} '{print $3}' > "!{fam}"_selected_ref_mtmalign_"$x".tcs.avg


	paste "!{fam}"_selected.tcs "!{fam}"_selected_ref_3dcoffee_"$x".tcs "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".tcs "!{fam}"_selected_ref_mtmalign_"$x".tcs > new_selected.tcs
	mv new_selected.tcs "!{fam}"_selected.tcs
		
	paste "!{fam}"_selected.tcs.avg "!{fam}"_selected_ref_3dcoffee_"$x".tcs.avg "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".tcs.avg "!{fam}"_selected_ref_mtmalign_"$x".tcs.avg > new_selected.tcs.avg
	mv new_selected.tcs.avg "!{fam}"_selected.tcs.avg

	for pref in "!{fam}"_selected_ref_3dcoffee_"$x" "!{fam}"_selected_ref_3dcoffee_TMalign_"$x" "!{fam}"_selected_ref_mtmalign_"$x"; do
		t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep TOT | awk -F '\t' -v var=!{fam} '{print var"\t"$4}' > "$pref".pid
		t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep AVG | head -n -1 | awk -F '\t' -v var=!{fam} '{print $3"\t"$5"\t"var}' > "$pref".pid.avg
		t_coffee -other_pg seq_reformat -in "$pref".fa_aln -output sim | grep BOT | awk -v var=!{fam} '{print $5"\t"$6"\t"$7"\t"var}' > "$pref".pid.pair

		paste "!{fam}"_selected.pid "$pref".pid > new_selected.pid
        mv new_selected.pid "!{fam}"_selected.pid
		paste "!{fam}"_selected.pid.avg "$pref".pid.avg > new_selected.pid.avg
        mv new_selected.pid.avg "!{fam}"_selected.pid.avg
		paste "!{fam}"_selected.pid.pair "$pref".pid.pair > new_selected.pid.pair
        mv new_selected.pid.pair "!{fam}"_selected.pid.pair

		echo -e !{fam}"\t"`esl-alistat "$pref".fa_aln | grep "Alignment length" | awk '{print $3}'` > "$pref".len
		paste "!{fam}"_selected.len "$pref".len > new_selected.len
		mv new_selected.len "!{fam}"_selected.len
	done


	for ref_struct in "" _"$x"; do	
		for pref in "!{fam}"_selected_ref_3dcoffee_"$x" "!{fam}"_selected_ref_3dcoffee_TMalign_"$x" "!{fam}"_selected_ref_mtmalign_"$x"; do
			t_coffee -other_pg irmsd "$pref".fa_aln -template_file "!{fam}"_selected_ref"$ref_struct".template_list > "$pref".all_irmsd"$ref_struct"
			grep -P "TOTAL\s+iRMSD" "$pref".all_irmsd"$ref_struct" | awk '{print $3}' > "$pref".irmsd"$ref_struct"
			grep -P "TOTAL\s+NiRMSD" "$pref".all_irmsd"$ref_struct" | awk '{print $3}' > "$pref".nirmsd"$ref_struct"
			grep -P "\sAVERAGE\s+NiRMS" "$pref".all_irmsd"$ref_struct" | awk -v var=!{fam} '{print $5"\t"var"\t"$3}' | tr -d "[]" > "$pref".nirmsd.avg"$ref_struct"
			grep -P "\sPAIRWISE\s+NiRMSD" "$pref".all_irmsd"$ref_struct" | awk -v var=!{fam} '{print $5"\t"$7"\t"var"\t"$3}' | tr -d "[]" > "$pref".nirmsd.pair"$ref_struct"
			grep -P "TOTAL\s+EVAL" "$pref".all_irmsd"$ref_struct" | awk '{print $3}' > "$pref".eval"$ref_struct"

			paste "!{fam}"_selected.irmsd "$pref".irmsd"$ref_struct" > new_selected.irmsd
			mv new_selected.irmsd "!{fam}"_selected.irmsd
			paste "!{fam}"_selected.nirmsd "$pref".nirmsd"$ref_struct" > new_selected.nirmsd
            mv new_selected.nirmsd "!{fam}"_selected.nirmsd
			paste "!{fam}"_selected.nirmsd.avg "$pref".nirmsd.avg"$ref_struct" > new_selected.nirmsd.avg
            mv new_selected.nirmsd.avg "!{fam}"_selected.nirmsd.avg
			paste "!{fam}"_selected.nirmsd.pair "$pref".nirmsd.pair"$ref_struct" > new_selected.nirmsd.pair
            mv new_selected.nirmsd.pair "!{fam}"_selected.nirmsd.pair
			paste "!{fam}"_selected.eval "$pref".eval"$ref_struct" > new_selected.eval
			mv new_selected.eval "!{fam}"_selected.eval
		done
	done
done
