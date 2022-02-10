for fam in `cat list_of_fams`; do
	cd $fam
	for seq in `cat "$fam".list`; do
		dir_name=${seq//\//_}
		echo -e $seq"\t"`extract_plddts.py $dir_name/ranking_debug.json`"\t"$fam
	done
	cd ../
done
