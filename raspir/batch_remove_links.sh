#/bin/bash

#unlink certain files found by ls command

for i in `ls *.bam`

	do
	unlink $i
	unlink $i.bai

done







