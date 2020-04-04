while getopts  ":i:f:t:h" o; do

    case "${o}" in

        i) i=${OPTARG}
            ;;
	f) f=${OPTARG}
            ;;
	t) t=${OPTARG}
            ;;
	h) echo "
			This script is used to align nucleotide MSA based on the relative protein MSA by a combination of transeq & mafft & pal2nal (conda insatallation tested).

			List of non-optional arguments:

                        -i	input folder, containing interleaved or onliner .fa / .fas /.fasta files (both interleaved and oneliner entires are allowed in the same file).
			-f	output folder, just specify a name to have it inside the current folder.
			-t	number of threads.
                        -h	help page.
"
               exit
          ;;
       \?) echo "WARNING! -$OPTARG isn't a valid option"
           exit
          ;;
       :) echo "WARNING! missing -$OPTARG value"
          exit
          ;;
       esac
 done
if [ -z "$i" ] || [ -z "$t" ] || [ -z "$f" ]
then
 echo " WARNING! non-optional argument/s is missing "
exit
fi
  
################################################################################################################
	
	mkdir ./$f

	cp $i/*.fa* $f

	cd $f

	for j in *.fa*; do export n=$(echo $j | awk -F "." '{print $1}');
		
		awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $j > "$n".n.fasta;

		rm $j;
		
		transeq -sequence $n.n.fasta -outseq $n.p.fasta;

		mafft $n.p.fasta > $n.mafft.p.aln;

		pal2nal.pl $n.mafft.p.aln $n.n.fasta -output fasta >> tmp1.aln;

		awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < tmp1.aln > tmp2.aln;

		for i in tmp2.aln;

			 do sed -s 's/_1//g' tmp2.aln > $n".mafft.n.aln"; done

		rm tmp*.aln;

	done;