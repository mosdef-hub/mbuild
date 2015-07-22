# --------------------------------------------------------------------
# Convert iPython notebooks to htmls.
# If they contain visualization code, generate visualization code too.
# --------------------------------------------------------------------

cd ../../docs/tutorials
for FILE in *.ipynb
do
	FILENAME="${FILE%.*}"
#	cp $FILE ${FILENAME}_tmp.ipynb

	sed -e 's/.visualize()/.visualize(export_topology=True)/g' $FILE > ${FILENAME}_vis.ipynb
	sed -i 's/.visualize(show_ports=True)/.visualize(show_ports=True, export_topology=True)/g' ${FILENAME}_vis.ipynb
	


done