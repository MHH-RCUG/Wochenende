#/bin/bash

#link certain files found by ls command

for i in `ls ../*mq30.bam`

	do
	# softlink bai and bam to current dir
	ln -s $i .
	ln -s $i.bai .

done







