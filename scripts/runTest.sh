#!/bin/sh

for N in 100 1000 1000000 100000000
do
	for seed in 0 1 2 3 4 5
	do
		./benchmark $N $seed
		if [ $? != 0 ]
		then
			echo "Test failed"
			exit 1
		fi
	done
done

