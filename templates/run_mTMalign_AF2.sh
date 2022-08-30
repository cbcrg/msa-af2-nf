for x in alphafold; do
		for i in `find *."$x".pdb`; do extract_from_pdb -force -infile $i > test.pdb; mv test.pdb $i; done
		for i in `cat !{list}`; do echo -e "${i//\//_}"'.'$x'.pdb'; done > !{fam}_selected_ref_"$x".mtmalign
		mTM-align -i !{fam}_selected_ref_"$x".mtmalign -o !{fam}_selected_ref_mtmalign_"$x"
		for i in `cat !{list}`; do
			pdb_name="${i//\//_}.$x.pdb"
			sed -i "s|$pdb_name|$i|g" result.fasta
		done
		mv result.fasta "!{fam}"_selected_ref_mtmalign_"$x".fa_aln
		sort_fasta.py "!{fam}"_selected_ref_mtmalign_"$x".fa_aln > tmp
		mv tmp "!{fam}"_selected_ref_mtmalign_"$x".fa_aln
	done