for f in ./data/size_scaling/*.in
do
	./simulation $f ${f%.*}.out
	echo ${f%.*}
done
