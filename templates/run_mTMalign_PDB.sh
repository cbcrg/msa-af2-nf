for i in `cat !{list}`; do echo -e `grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`'.pdb'; done > "!{fam}"_selected_ref.mtmalign
	mTM-align -i "!{fam}"_selected_ref.mtmalign -o "!{fam}"_selected_ref_mtmalign
	for i in `cat !{list}`; do
	ref_pdb=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`
	ref_seq=`grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $1}'`
	sed -i "s|>$ref_pdb\.pdb|$ref_seq|g" result.fasta
	done
	mv result.fasta "!{fam}"_selected_ref_mtmalign.fa_aln
	sort_fasta.py "!{fam}"_selected_ref_mtmalign.fa_aln > tmp
	mv tmp "!{fam}"_selected_ref_mtmalign.fa_aln