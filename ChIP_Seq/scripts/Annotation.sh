# Annotation

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')

tput setaf 3; tput bold; echo " "
tput setaf 3; tput bold; echo "    Annotating Regions in The Genome with Homer"
tput setaf 2; tput bold; echo " "

annotation=$(grep "ANNOTATION" config | cut -d ":" -f2)

case $annotation in

TRUE)

fasta=$(grep "PATH TO FASTA" config | cut -d ":" -f2)
gtf=$(grep "PATH TO GTF" config | cut -d ":" -f2)

pwd=$(pwd)

if [ -d "$pwd/peak_calling/merged" ]

	then

		ln -s peak_calling/merged/* .

		for peaks in $(ls *.tsv)

			do

			sample=$(echo $peaks | cut -d "_" -f2)
			out=$(echo "annotation_$sample")

			annotatePeaks.pl $peaks $fasta -gtf $gtf > $out

			done

			tput setaf 2; tput bold; echo " "
			tput setaf 3; tput bold; echo "    Annotation Complete!"
			tput setaf 2; tput bold; echo " "

	else

		ln -s peak_calling/*.broadPeak .

		for peaks in $(ls *.broadPeak)

                        do

                        sample=$(echo $peaks | cut -d "_" -f2)
                        out=$(echo "annotation_$sample")

                        annotatePeaks.pl $peaks $fasta -gtf $gtf > $out

                        done

                        tput setaf 2; tput bold; echo " "
                        tput setaf 3; tput bold; echo "    Annotation Complete!"
                        tput setaf 2; tput bold; echo " "

fi

;;

FALSE)

;;

*)

esac

mkdir annotated
mv annotation* annotated/
rm *.broadPeak
