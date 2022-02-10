for x in dmpfold alphafold; do
        for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
done
tcs_prev=0
for (( tcs=50; tcs<=95; tcs+=5 )); do
        awk -v var=$tcs '{ if ($9 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_3dcoffee_dmpfold
        awk -v var=$tcs '{ if ($10 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_3dcoffee_TMalign_dmpfold
        awk -v var=$tcs '{ if ($11 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_mtmalign_dmpfold
        awk -v var=$tcs '{ if ($12 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_3dcoffee_alphafold
        awk -v var=$tcs '{ if ($13 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_3dcoffee_TMalign_alphafold
        awk -v var=$tcs '{ if ($14 >= var) print $1}' !{avg_tcs} > "!{fam}".list_tcs_above_"$tcs"_ref_mtmalign_alphafold
        for tcs_ref in mtmalign_dmpfold mtmalign_alphafold 3dcoffee_dmpfold 3dcoffee_alphafold; do # 3dcoffee_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold
                if [ `wc -l "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref" | awk '{print $1}'` -lt 2 ]; then
                                for mode in sp tc; do
		                        for test_aln in famsa tcoffee msaprobs ginsi psicoffee 3dcoffee_dmpfold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_alphafold mtmalign_alphafold; do #deepblast
                                                for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
                                                        echo "NA" > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref"
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref"
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref"

                                                        echo "NA" > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref".same_state
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref".same_state
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref".same_state

                                                        echo "NA" > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
                                                        echo > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
                                                done
                                        done
		                done
                else
                        if [ "$tcs_prev" -gt 45 ] && cmp --silent -- "!{fam}".list_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"; then 
                                cp "!{fam}"_selected_ref.template_list_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref.template_list_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref.mtmalign_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref.mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref.fa_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_famsa.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_famsa.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_tcoffee.lib_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_tcoffee.lib_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_tcoffee.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_tcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_msaprobs.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_msaprobs.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_ginsi.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_ginsi.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_psicoffee.lib_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_psicoffee.lib_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_psicoffee.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_psicoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_mtmalign_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_mtmalign.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_mtmalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}".dssp_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for aln in sap TMalign; do
                                         cp "!{fam}"_selected_ref."$aln".lib_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref."$aln".lib_tcs_above_"$tcs"_ref_"$tcs_ref"
                                done
                                cp "!{fam}"_selected_ref_3dcoffee.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_3dcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                cp "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for x in dmpfold alphafold; do
		                        for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do echo -e ">$i" "_P_" "${i//\//_}"."$x"; done > "!{fam}"_selected_ref_"$x".template_list_tcs_above_"$tcs"_ref_"$tcs_ref"
                                        for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do echo -e "${i//\//_}"."$x".pdb; done > "!{fam}"_selected_ref_"$x".mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                        cp "!{fam}"_selected_ref_mtmalign_"$x"_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_mtmalign_"$x"_tcs_above_"$tcs"_ref_"$tcs_ref"
                                        cp "!{fam}"_selected_ref_mtmalign_"$x".fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_mtmalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                        for aln in sap TMalign; do
                                               cp "!{fam}"_selected_ref_"$x"."$aln".lib_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_"$x"."$aln".lib_tcs_above_"$tcs"_ref_"$tcs_ref" 
                                        done
                                        cp "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                        cp "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln_tcs_above_"$tcs_prev"_ref_"$tcs_ref" "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                done
                        else
                                
                                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do grep "$i" !{template}; done > "!{fam}"_selected_ref.template_list_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do echo -e `grep "$i" !{template} | awk '{print $3}'`'.pdb'; done > "!{fam}"_selected_ref.mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                xargs samtools faidx !{fasta} < "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref" > "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref"
                                xargs samtools faidx !{dssp_fasta} < "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref" > "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref"

                                famsa "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" "!{fam}"_selected_ref_famsa.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
	                        sort_fasta.py "!{fam}"_selected_ref_famsa.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
	                        mv tmp "!{fam}"_selected_ref_famsa.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                        
                                t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -output fasta_aln -out_lib "!{fam}"_selected_ref_tcoffee.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -outfile "!{fam}"_selected_ref_tcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
	                        sort_fasta.py "!{fam}"_selected_ref_tcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
	                        mv tmp "!{fam}"_selected_ref_tcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                msaprobs -o "!{fam}"_selected_ref_msaprobs.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref"
			        sort_fasta.py "!{fam}"_selected_ref_msaprobs.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_msaprobs.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
			        ginsi "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" > "!{fam}"_selected_ref_ginsi.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
			        sort_fasta.py "!{fam}"_selected_ref_ginsi.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_ginsi.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                t_coffee -seq "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -special_mode psicoffee -blast_server LOCAL -protein_db !{db} -psitrim 100 -psiJ 3 -prot_min_cov 90 -prot_max_sim 100 -prot_min_sim 0 -output fasta_aln -out_lib "!{fam}"_selected_ref_psicoffee.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -outfile "!{fam}"_selected_ref_psicoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
			        sort_fasta.py "!{fam}"_selected_ref_psicoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_psicoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

			        #t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs" -method /tcoffee/t_coffee/src/deepblast.tc_method -output fasta_aln -out_lib "!{fam}"_selected_ref_deepblast.lib_tcs_above_"$tcs" -outfile "!{fam}"_selected_ref_deepblast.fa_aln_tcs_above_"$tcs" -thread 4
	                        #sort_fasta.py "!{fam}"_selected_ref_deepblast.fa_aln_tcs_above_"$tcs" > tmp
        	                #mv tmp "!{fam}"_selected_ref_deepblast.fa_aln_tcs_above_"$tcs"	

                                mTM-align -i "!{fam}"_selected_ref.mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref" -o "!{fam}"_selected_ref_mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do
                                        ref_pdb=`grep "$i" !{template} | awk '{print $3}'`
                                        ref_seq=`grep "$i" !{template} | awk '{print $1}'`
                                        sed -i "s|>$ref_pdb\.pdb|$ref_seq|g" result.fasta
                                done
                                mv result.fasta "!{fam}"_selected_ref_mtmalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
		                sort_fasta.py "!{fam}"_selected_ref_mtmalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
		                mv tmp "!{fam}"_selected_ref_mtmalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                for aln in sap TMalign; do
                                        t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -method "$aln"_pair -template_file "!{fam}"_selected_ref.template_list_tcs_above_"$tcs"_ref_"$tcs_ref" -out_lib "!{fam}"_selected_ref."$aln".lib_tcs_above_"$tcs"_ref_"$tcs_ref"
                                done
                                t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -lib "!{fam}"_selected_ref.sap.lib_tcs_above_"$tcs"_ref_"$tcs_ref" "!{fam}"_selected_ref.TMalign.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
			        sort_fasta.py "!{fam}"_selected_ref_3dcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_3dcoffee.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -lib "!{fam}"_selected_ref.TMalign.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
			        sort_fasta.py "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_3dcoffee_TMalign.fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                for x in dmpfold alphafold; do
		                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do echo -e ">$i" "_P_" "${i//\//_}"."$x"; done > "!{fam}"_selected_ref_"$x".template_list_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do echo -e "${i//\//_}"."$x".pdb; done > "!{fam}"_selected_ref_"$x".mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref"
                                mTM-align -i "!{fam}"_selected_ref_"$x".mtmalign_tcs_above_"$tcs"_ref_"$tcs_ref" -o "!{fam}"_selected_ref_mtmalign_"$x"_tcs_above_"$tcs"_ref_"$tcs_ref"
                                for i in `cat "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref"`; do
                                        pdb_name="${i//\//_}.$x.pdb"
                                        sed -i "s|$pdb_name|$i|g" result.fasta
                                done
                                mv result.fasta "!{fam}"_selected_ref_mtmalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
		                sort_fasta.py "!{fam}"_selected_ref_mtmalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
		                mv tmp "!{fam}"_selected_ref_mtmalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                for aln in sap TMalign; do
                                        t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -method "$aln"_pair -template_file "!{fam}"_selected_ref_"$x".template_list_tcs_above_"$tcs"_ref_"$tcs_ref" -out_lib "!{fam}"_selected_ref_"$x"."$aln".lib_tcs_above_"$tcs"_ref_"$tcs_ref"
                                done
                                t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -lib "!{fam}"_selected_ref_"$x".sap.lib_tcs_above_"$tcs"_ref_"$tcs_ref" "!{fam}"_selected_ref_"$x".TMalign.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
			        sort_fasta.py "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_3dcoffee_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"

                                t_coffee "!{fam}"_selected_ref.fa_tcs_above_"$tcs"_ref_"$tcs_ref" -lib "!{fam}"_selected_ref_"$x".TMalign.lib_tcs_above_"$tcs"_ref_"$tcs_ref" -output fasta_aln -outfile "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -thread 4
			        sort_fasta.py "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" > tmp
			        mv tmp "!{fam}"_selected_ref_3dcoffee_TMalign_"$x".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref"
                                done

                        fi
                        for mode in sp tc; do
                                for ref_aln in 3dcoffee 3dcoffee_TMalign mtmalign; do
			                for test_aln in famsa tcoffee msaprobs ginsi psicoffee 3dcoffee_dmpfold 3dcoffee_TMalign_dmpfold mtmalign_dmpfold 3dcoffee_alphafold 3dcoffee_TMalign_alphafold mtmalign_alphafold; do #deepblast
                                                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" | grep -v "seq1" | grep -v "*" | awk -v var=!{fam} '{print var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref"
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -io_format s | awk -v var=!{fam} '{print $1"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref"
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -io_format p | awk -v var=!{fam} '{print $1"\t"$2"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref"

                                                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]=[Struc]' | grep -v "seq1" | grep -v "*" | awk -v var=!{fam} '{print var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref".same_state
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]=[Struc]' -io_format s | awk -v var=!{fam} '{print $1"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref".same_state
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]=[Struc]' -io_format p | awk -v var=!{fam} '{print $1"\t"$2"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref".same_state

                                                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]+[H][E]=[Struc]' | grep -v "seq1" | grep -v "*" | awk -v var=!{fam} '{print var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode"_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]+[H][E]=[Struc]' -io_format s | awk -v var=!{fam} '{print $1"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".avg_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops
                        	                t_coffee -other_pg aln_compare -al1 "!{fam}"_selected_ref_"$ref_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -al2 "!{fam}"_selected_ref_"$test_aln".fa_aln_tcs_above_"$tcs"_ref_"$tcs_ref" -compare_mode "$mode" -st "!{fam}".dssp_tcs_above_"$tcs"_ref_"$tcs_ref" pep -io_cat '[H][H]+[E][E]+[H][E]=[Struc]' -io_format p | awk -v var=!{fam} '{print $1"\t"$2"\t"var"\t"$4}' > "!{fam}"_selected_ref_"$test_aln"_vs_ref_"$ref_aln"."$mode".pair_tcs_above_"$tcs"_ref_"$tcs_ref".without_loops

                                        done
                                done
                        done
                fi
                echo `wc -l "!{fam}".list_tcs_above_"$tcs"_ref_"$tcs_ref" | awk '{print $1}'` >> "!{fam}"_selected.list_tcs_above_"$tcs"_ref_"$tcs_ref"
        done
        tcs_prev=$tcs
done

