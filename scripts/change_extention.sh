# Rename all *.m to *.vtk
for f in *.m; do
	vtk=${f%.m}.vtk
    echo "m2vtk $f $vtk"
	m2vtk $f $vtk
done
