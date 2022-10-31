#!/bin/bash

peak_c=$(grep "PEAK_CALLING_MACS" config | cut -d ":" -f2)

case $peak_c in

TRUE)

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Peak Calling Analysis with MACS2"
	tput setaf 2; tput bold; echo " "

ln -s downsampled/*.bam .
ln -s downsampled/*.bam.bai .

	while [[ $treatment != "EXIT" ]]

        do

        tput setaf 6; tput bold; echo "Set treatment bam file. To abord peak call analysis type (EXIT)"
        tput setaf 2; tput bold; echo " "
        read treatment

		if [[ $treatment = "EXIT" ]]

			then
                       	break

			else
			tput setaf 6; tput bold; echo "Set input bam file. If there is not input file type (no_input) "
                        tput setaf 2; tput bold; echo " "
                        read input

			tput setaf 6; tput bold; echo "Set output file, e.g sample name"
                        tput setaf 2; tput bold; echo " "
                        read out

			g=$(grep "GENOME_SP" config | cut -d ":" -f2)
			broad=$(grep "BROAD" config | cut -d ":" -f2)


			                        case $broad in

		                                TRUE)
						if [[ $input != "no_input" ]]

							then
                	                		# broad (with input)
                        		        	macs2 callpeak --treatment $treatment --control $input --nomodel --broad --format BAM --gsize $g -n $out
							# Filtering (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
							grep -v "#" "${out}_peaks.xls" | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' > "${out}_peaks.filtered.xls"
							else
							macs2 callpeak --treatment $treatment --nomodel --broad --format BAM --gsize $g -n $out
						fi

		                                ;;

                		                FALSE)
						if [[ $input != "no_input" ]]

                                                        then
	      	                        		# narrow (with input)
			                                macs2 callpeak --treatment $treatment --control $input --nomodel --format BAM --gsize $g -n $out
							# Filtering (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
                                                        grep -v "#" "${out}_peaks.xls" | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' > "${out}_peaks.filtered.xls"
							else
							macs2 callpeak --treatment $treatment --nomodel --format BAM --gsize $g -n $out
						fi

                		                ;;

        	                	        *)

                        	       		 ;;

			                        esac
                fi

done

mkdir peak_calling
# narrow
mv *.narrowPeak peak_calling
mv *.xls peak_calling
# broad
mv *.broadPeak peak_calling
mv *.gappedPeak peak_calling


;;

FALSE)

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    Peak Calling Analysis with SICER"
tput setaf 2; tput bold; echo " "

ln -s downsampled/*.bam .
ln -s downsampled/*.bam.bai .
mkdir peak_calling

	while [[ $treatment != "EXIT" ]]

        do

        tput setaf 6; tput bold; echo "Set treatment bam file. To abord peak call analysis type (EXIT)"
        tput setaf 2; tput bold; echo " "
        read treatment

		if [[ $treatment = "EXIT" ]]

			then
                       	break

			else
			tput setaf 6; tput bold; echo "Set input bam file. If there is not input file type (no_input) "
                        tput setaf 2; tput bold; echo " "
                        read input

			tput setaf 6; tput bold; echo "Set output file, e.g sample name"
                        tput setaf 2; tput bold; echo " "
                        read out

				if [[ $input != "no_input" ]]

					then
                	                # (with input)
					sicer -t $treatment -c $input -s mm10 -w 200 -g 600 --false_discovery_rate 0.001 --effective_genome_fraction 0.88 -o $out -cpu 7
					mv $out peak_calling
					else
					sicer -t $treatment -s mm10 -w 200 -g 600 --false_discovery_rate 0.001 --effective_genome_fraction 0.88 -o $out -cpu 7
					mv $out peak_calling
				fi

                fi

done

;;

*)

esac
rm *.bam
rm *.bai
