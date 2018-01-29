for f in ./data/*.in
do
	./simulation $f ${f%.*}.out
done
