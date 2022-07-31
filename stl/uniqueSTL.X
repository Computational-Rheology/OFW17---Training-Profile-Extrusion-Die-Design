# changest the stl to give the name as for the OF patch,
# then produces an unique file for the grid production

for N in `cat mylist`; do
	cat ${N}.stl | sed -e '1d;$d' > ${N}_.tmp
	echo "solid ${N}" >> ${N}_1.stl
	cat ${N}_.tmp >> ${N}_1.stl
	echo "endsolid ${N}" >> ${N}_1.stl
        touch total.stl
        surfaceAdd ${N}_1.stl total.stl total.stl
	rm ${N}_.tmp
done

surfaceFeatureEdges total.stl total.fms

