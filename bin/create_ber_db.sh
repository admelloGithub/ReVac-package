#!/bin/sh

print_usage()
{
	progname=`basename $0`
	cat << END
usage: $progname -i <btab_hits> -o <nucdb_fasta_out>
	-d <db> -f <db_format> -b <bin_dir> -y <fetch_bin_dir> [-p] [-m <prot_nuc_id_map>]
END
	exit 1
}

# Default to nucleotide seqs
protein=F

while getopts "i:o:d:f:m:b:c:a:x:p" opt
do
	case $opt in
		i) hits=$OPTARG;;
		o) out=$OPTARG;;
		d) db=$OPTARG;;
		f) db_format=$OPTARG;;
		m) id_map=$OPTARG;;
		b) bin_dir=$OPTARG;;
	        c) cdbfasta_path=$OPTARG;;
		a) formatdb_path=$OPTARG;;
		x) xdformat_path=$OPTARG;;
		p) protein='T';;
	esac
done

test -z $hits && echo "No btab hits provided" && print_usage
test -z $out &&	echo "No output nucleotide fasta provided" && print_usage
test -z $db && echo "No db provided" && print_usage
test -z $db_format && echo "No db format provided" && print_usage
test -z $bin_dir && echo "No bin directory provided" && print_usage
test $protein == 'F' -a -z "$id_map" && echo "No protein/nucleotide map provided" && print_usage

if [ -f $hits ] && [ ! -s $hits ]
then
	touch $out
	exit 0
fi

temp_file_ext='.id_list'
id_temp_file=$out$temp_file_ext

id_to_fetch=""
# Collect IDs, separate by space, and store to id_to_fetch variable
if [ $protein == 'F' ]
then
	for prot_id in `cut -f1 $hits | sort -u`; do
		if [ $prot_id ]
		then
			hit_id=`grep "${prot_id//\./\\\.}	" $id_map | cut -f2 | sort -u`
			id_to_fetch="$hit_id $id_to_fetch"
		fi
	done
else
	id_to_fetch=`cut -f6 $hits | sort -u | perl -pe 's/\n/ /g'`
fi

# Pass IDs to a file
echo $id_to_fetch 1>$id_temp_file

# Call fetch_fasta_from_db, using a list file of IDs
if [ "$id_temp_file" ]
then
	$bin_dir/fetch_fasta_from_db -I "$id_temp_file" -d $db -p $protein -f $db_format -o "$out" --cdbfasta_path "$cdbfasta_path" --formatdb_path "$formatdb_path" --xdformat_path "$xdformat_path"
	ec=$?
	if [ $ec -ne 0 ]
	then
		echo "Error fetching sequences"
		exit $ec
	fi
else
	exit 1
fi
