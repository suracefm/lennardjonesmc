for f in ./data/phasediag2/*.in
do
	./simulation $f ${f%.*}.out
	echo ${f%.*}
done
