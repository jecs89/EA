	for iter in 500
		do
		echo $iter iterations
		for f in 1 2 3 4 5
			do
			for i in 100 250 500 750 1000 
				do
				./qiear  $f $i 0 1 $iter 4 0.00000000000000001
			done
		done
	done