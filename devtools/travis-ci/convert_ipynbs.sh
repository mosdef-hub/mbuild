# --------------------------------------------------------------------
# Convert iPython notebooks to htmls.
# If they contain visualization code, generate visualization code too.
# --------------------------------------------------------------------

cd ../../docs/tutorials
TPL_PATH=`cd ../assets/tpl; pwd`

# create build folder



for FILE in *.ipynb
do
	FILENAME="${FILE%.*}"
#	cp $FILE ${FILENAME}_tmp.ipynb

	sed -e 's/.visualize()/.visualize(export_topology=True)/g' $FILE > ${FILENAME}_tmp.ipynb
	sed -e 's/.visualize(show_ports=True)/.visualize(show_ports=True, export_topology=True)/g' ${FILENAME}_tmp.ipynb > ${FILENAME}_vis.ipynb
	rm ${FILENAME}_tmp.ipynb

	echo "Running ${FILENAME} notebook non-interactively..."
	ipython nbconvert --to=html --ExecutePreprocessor.enabled=True ${FILENAME}_vis.ipynb > /dev/null
	rm ${FILENAME}_vis.html
	rm ${FILENAME}_vis.ipynb

	mv topology_*.json ../assets/json/

	#pwd

	ipython nbconvert --TemplateExporter.template_path="['.', '../assets/tpl/']" --to html --template mbuild_ipynb_template.tpl $FILE

done