while getopts ":i:f:d:n:c:t:h" o; do
    case "${o}" in

	i) i=${OPTARG}
            ;;
	f) f=${OPTARG}
            ;;
        d) d=${OPTARG}
            ;;
        n) n=${OPTARG}
            ;;
	c) c=${OPTARG}
            ;;
	t) t=${OPTARG}
            ;;
	h) echo "
                        This script will concatenate all the aligned one-liner fasta formatted .aln files in the input folder using Phyutility and write a partition file.

                        List of non-optional arguments:

                        -i path to the input folder containing aligned .aln files (write ./ to launch the script in the current folder).
                        -f path to the output folder.
			-d data type (DNA, AA, ..) or evolutionary model (i.e. WAG for amminocaids)
                        -n output suffix.
                        -c if TRUE also a codon-formatted partition file is generated.
                        -t if NEXUS the concatenation file will be written as .nxs nesxus formatted file, if PHYLIP also a .phy pylip concatenation file will be generated.
			-h help page.
"
               exit
           ;;
         esac
     done

if [ -z "$i" ] || [ -z "$f" ] || [ -z "$d" ] || [ -z "$n" ] || [ -z "$c" ] || [ -z "$t" ]
then
 echo " WARNING! non-optional argument/s is missing "
exit
fi

#####################################################################################################################

mkdir $f

cp $i/*.aln $f

cd $f

for i in *.aln; 

	do echo -n ' ' $i >> tmp.prt; 

done; 

java -jar ~/phyutility.jar -concat -in $(cat tmp.prt) -out $n"_concatenation.nxs"

	rm tmp.prt

export length=$(ls *.aln | wc -l)

	head -3 $n"_concatenation.nxs" | tail -1 | awk '{ for (i=3;i<=NF;i+=2) $i="" } 2' >> a.tmp;

	awk -F " " '{$1 = ""; print $0}' a.tmp > a_cor.tmp
	
	head -3 $n"_concatenation.nxs" | tail -1 | awk '{ for (i=2;i<=NF;i+=2) $i="" } 2' >> b.tmp;

	sed -s 's/[][]//g' b.tmp > b_cor.tmp

for j in $(eval echo {1..$length}); 
	
	do export a=$(awk '{print $'$j'}' a_cor.tmp); 

	export b=$(awk '{print $'$j'}' b_cor.tmp);
 
	echo $d"," $b "=" $a >> $n"_genes_partitions.txt"; 
	
	done;

rm a.tmp b.tmp a_cor.tmp b_cor.tmp *.aln phyutility.log

#####################################################################################################################


if [ "$c" == "TRUE" ];

	then

#
	sed -s 's/gene1/gene/' $n"_genes_partitions.txt" >> tmp_codon_partitions.txt

	while read i; 

	do export r=$(echo $i | awk -F " " '{print $4}' | awk -F "-" '{print $1}'); export m=$((r+1)); export l=$((r+2)); export q=$(echo $i | awk -F " " '{print $2}');

	echo $i"/3" | sed 's/'$q'/'$q'_st/'>> $n"_codon_partitions.txt";
	echo $i"/3" | sed 's/'$q'/'$q'_nd/' | sed 's/'$r'/'$m'/' >> $n"_codon_partitions.txt"; 
	echo $i"/3" | sed 's/'$q'/'$q'_rd/' | sed 's/'$r'/'$l'/' >> $n"_codon_partitions.txt";

	done < tmp_codon_partitions.txt;

	rm tmp_codon_partitions.txt
fi

#####################################################################################################################

if [ "$t" == "PHYLIP" ];

	then
	
	ntax=$(grep "DIMENSIONS" $n"_concatenation.nxs" | awk -F "=" '{print $2}' | awk '{print $1}');

	nchar=$(grep "DIMENSIONS" $n"_concatenation.nxs" | awk -F "=" '{print $3}' | awk -F ";" '{print $1}');

	head -$(( ntax + 6 )) $n"_concatenation.nxs" | tail -$ntax >> $n"_concatenation.phy.tmp"

	echo -e "$ntax $nchar \n$(cat $n"_concatenation.phy.tmp")" > $n"_concatenation.phy"

	rm $n"_concatenation.phy.tmp" $n"_concatenation.nxs"

fi
