#determining if coral primers match ctenophore sequences

mamba activate bioinformatics #activate environment
cd ~/Desktop/AW/data #set current directory
ls #list directory contents

#identifies only the ctenophore sequences containing F and R primer with 2 or fewer mismatches
cutadapt -g CGTGAAACCGYTRRAAGGG...CGTCTTGAAACACGGACCAA -e 2 -O 10 --action=retain --trimmed-only -o retained_cteno.fasta 28S_Ctenos_with_Outgroups.fasta > retained_cteno.txt
mv retained_cteno.fasta ~/Desktop/AW/results
mv retained_cteno.txt ~/Desktop/AW/results

#identifies and trims only the ctenophore sequences containing F and R primer with 2 or fewer mismatches
cutadapt -g CGTGAAACCGYTRRAAGGG...CGTCTTGAAACACGGACCAA -e 2 -O 10 --trimmed-only -o trimmed_cteno.fasta 28S_Ctenos_with_Outgroups.fasta > trimmed_cteno.txt
mv trimmed_cteno.fasta ~/Desktop/AW/results
mv trimmed_cteno.txt ~/Desktop/AW/results

#generating match list of ctenophore ASVs to curate with LULU

#converting the coral seqs in a .tsv format to a .fasta file
awk -F '\t' '{printf ">%s\n%s\n",$1,$2}' Cteno_preclassifier.tsv > Cteno_preclassifier.fasta 

#creating a matchlist by first creating a BLAST database...
makeblastdb -in Cteno_preclassifier.fasta -parse_seqids -dbtype nucl 

